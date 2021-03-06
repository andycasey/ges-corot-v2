#!/usr/bin/python

"""
Propagate relevant flag information from one node to others.
"""

import logging
import numpy as np
import yaml

from code import GESDatabase
from code.gesdb import UnknownNodeError

# Connect to database.
db_filename = "db.yaml"
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)

# Create a database object.
database = GESDatabase(**credentials)

logger = logging.getLogger("ges")


with open("flags.yaml", "r") as fp:
    qc_flags = yaml.load(fp)


# Clear any previous propagations before starting.
logger.info("Clearing previous propagations and setting all to have passed_quality_control = True")
database.update(
    """ UPDATE  results
           SET  propagated_tech_from_result_id = null,
                propagated_peculi_from_result_id = null,
                propagated_remark_from_result_id = null,
                propagated_tech = '',
                propagated_peculi = '',
                propagated_remark = '',
                passed_quality_control = true
         WHERE  passed_quality_control = false;""")

database.connection.commit()

# Identify spurious spectra and mark them as such.
N_peculiar_spectra = {}
peculiar_spectra_kwds = dict(
    and_or="and", sigma_discrepant=3, teff_discrepant=250, feh_discrepant=1.0)

for wg in (1, ):

    logger.info("Querying for peculiar spectra in WG{}".format(wg))

    kwds = dict(wg=wg)
    kwds.update(peculiar_spectra_kwds)

    peculiar_spectrum_query = """
        with t4 as (
        select id, cname, filename, avg_filename_teff, avg_cname_teff, stddev_cname_teff, abs((avg_cname_teff - avg_filename_teff)/(0.00001 + stddev_cname_teff)) as abs_sigma_teff_discrepant, avg_filename_feh, avg_cname_feh, stddev_cname_feh, abs((avg_cname_feh - avg_filename_feh)/(0.00001 + stddev_cname_feh)) as abs_sigma_feh_discrepant FROM (with ar as (select distinct on (filename) id, cname, trim(filename) as filename, avg(teff) over w as avg_filename_teff, avg(feh) over w as avg_filename_feh from (with n as (select id from nodes where wg = {wg}) select distinct on (r.filename, r.node_id) r.id, r.cname, trim(r.filename) as filename, r.node_id, teff, feh from n, results as r where r.node_id = n.id and teff <> 'NaN' or feh <> 'NaN') t window w as (partition by filename)) select ar.id, ar.cname, ar.filename, ar.avg_filename_teff, avg(avg_filename_teff) over w2 as avg_cname_teff, stddev(avg_filename_teff) over w2 as stddev_cname_teff, ar.avg_filename_feh, avg(avg_filename_feh) over w2 as avg_cname_feh, stddev(avg_filename_feh) over w2 as stddev_cname_feh FROM ar window w2 as (partition by cname)) t3)
        select * from t4 where (t4.abs_sigma_teff_discrepant > {sigma_discrepant} and t4.abs_sigma_teff_discrepant <> 'NaN' and abs(t4.avg_cname_teff - avg_filename_teff) >= {teff_discrepant}) {and_or} (t4.abs_sigma_feh_discrepant > {sigma_discrepant} and abs(t4.avg_cname_feh - t4.avg_filename_feh) >= {feh_discrepant} and t4.abs_sigma_feh_discrepant <> 'NaN') order by cname asc;""".format(
            **kwds)

    peculiar_spectra = database.retrieve_table(peculiar_spectrum_query)
    if peculiar_spectra is None: continue

    N_peculiar_spectra[wg] = len(peculiar_spectra)

    for row in peculiar_spectra:

        filenames = row["filename"].strip().split("|")
        logger.info("Propagating {}/{}/{}".format(
            row["id"], row["cname"], row["filename"]))

        n = 0
        for filename in filenames:
            n += database.update(
                """ UPDATE results
                       SET propagated_tech = '10106-{}-00-00-A',
                           propagated_tech_from_result_id = '{}',
                           passed_quality_control = false
                     WHERE filename LIKE '%{}%';""".format(
                        wg, int(row["id"]), filename))

        if n > 0:
            logger.info("--> affected {} results".format(n))

database.connection.commit()


# PROPAGATE BY SPECTRUM
def propagate_by_spectrum(flag, constraint=None, commit=False):
    
    constraint_str = "" if constraint is None else " AND {}".format(constraint)
    affected = database.retrieve_table(
        """ SELECT  id, TRIM(filename) AS filename, TRIM(tech) as tech
            FROM    results
            WHERE   (tech LIKE '%{0}-%' OR TRIM(TECH) = '{0}') {1}
        """.format(flag, constraint_str))
    
    if affected is None:
        return 0

    N = 0
    for row in affected:
    
        # Each row can have multiple TECH flags, so first identify the TECH flag
        # that PostgreSQL matched on.
        tech_flags = row["tech"].strip().split("|")
        for tech_flag in tech_flags:
            if "{}-".format(flag) in tech_flag:
                matched_tech_flag = tech_flag
                break

        else:
            raise ValueError(
                "cannot identify tech flag {} from the SQL match: {}".format(
                    flag, row["tech"].strip()))


        # Update other results using the same filename(s).
        filenames = row["filename"].strip().split("|")
        for j, filename in enumerate(filenames):
            if not filename: continue

            n = database.update(
                """ UPDATE  results
                    SET     propagated_tech_from_result_id = '{}',
                            propagated_tech = '{}',
                            passed_quality_control = false
                    WHERE   filename LIKE '%{}%'
                      AND   passed_quality_control = true
                """.format(
                    int(row["id"]), matched_tech_flag, filename))
            
            N += n

            if n > 0:
                logger.info("Propagated ({}/{}/{}) to {} other entries".format(
                    row["id"], matched_tech_flag, filename, n))    

    if commit:
        database.connection.commit()

    return N


# PROPAGATE BY SPECTRUM
N_propagations = {}

for key, value in qc_flags.get("propagate_flags_by_spectrum", {}).items():
    if key == "no_constraint":
        for flag in value:
            N_propagations.setdefault(flag, 0)
            N_propagations[flag] += propagate_by_spectrum(flag, None)
    
    else:
        flag = key
        constraint = value.get("constraint", None)
        N_propagations.setdefault(flag, 0)
        N_propagations[flag] += propagate_by_spectrum(flag, constraint)

database.connection.commit()


# PROPAGATE BY CNAME
def propagate_by_cname(flag, constraint=None, commit=False):

    constraint_str = "" if constraint is None else " AND {}".format(constraint)
    affected = database.retrieve_table(
        """ SELECT  id, cname, TRIM(tech) as tech
              FROM  results
             WHERE  (tech LIKE '%{0}-%' OR TRIM(tech) = '{0}') {1}
        """.format(flag, constraint_str))

    if affected is None:
        return 0

    N = 0
    for row in affected:
        # Each row can have multiple TECH flags, so first identify the TECH flag
        # that PostgreSQL matched on.
        tech_flags = row["tech"].strip().split("|")
        for tech_flag in tech_flags:
            if "{}-".format(flag) in tech_flag:
                matched_tech_flag = tech_flag
                break

        else:
            raise ValueError(
                "cannot identify tech flag {} from the SQL match: {}".format(
                    flag, row["tech"].strip()))

        # Update other results matching this CNAME.
        n = database.update(
            """ UPDATE  results
                   SET  propagated_tech_from_result_id = '{}',
                        propagated_tech = '{}',
                        passed_quality_control = false
                 WHERE  cname = '{}'
                   AND  passed_quality_control = true;
            """.format(int(row["id"]), matched_tech_flag, row["cname"]))
        N += n

        if n > 0:
            logger.info("Propagated ({}/{}/{}) to {} other entries".format(
                row["id"], matched_tech_flag, row["cname"], n))

    if commit:
        database.connection.commit()

    return N


# PROPAGATE BY CNAME
for key, value in qc_flags.get("propagate_flags_by_cname", {}).items():

    if key == "no_constraint":
        for flag in value:
            N_propagations.setdefault(flag, 0)
            N_propagations[flag] += propagate_by_cname(flag, None)

    else:
        flag = key
        constraint = value.get("constraint", None)
        N_propagations.setdefault(flag, 0)
        N_propagations[flag] += propagate_by_cname(flag, constraint)

database.connection.commit()


# NODE-SPECIFIC FLAGS
N_marked_as_poor_quality = {}
for key, value in qc_flags.get("node_specific_flags", {}).items():

    if key == "no_constraint":
        for flag in value:
            N_marked_as_poor_quality.setdefault(flag, 0)
            N = database.update(
                """ UPDATE  results
                       SET  passed_quality_control = false
                     WHERE  (tech LIKE '%{0}-%' OR trim(tech) = '{0}')
                       AND  passed_quality_control = true;
                """.format(flag))
            N_marked_as_poor_quality[flag] += N

            if N > 0:
                logger.info(
                    "Marked {} results as poor quality due to matching flag {}"\
                    .format(N, flag))
    else:
        flag = key
        constraint = value.get("constraint", None)
        constraint_str = "" if constraint is None else " AND {}".format(constraint)
        N_marked_as_poor_quality.setdefault(flag, 0)

        N = database.update(
            """ UPDATE  results
                   SET  passed_quality_control = false
                 WHERE  (tech LIKE '%{0}-%' OR trim(tech) = '{0}')
                   AND  passed_quality_control = true
                   {1};
            """.format(flag, constraint_str))
        N_marked_as_poor_quality[flag] += N

        if N > 0:
            logger.info(
                "Marked {} results as poor quality due to matching flag {} "\
                "and constraint {}".format(N, flag, constraint))


database.connection.commit()


# Removing all results for UVES/Damiani because the benchmark values are exactly
# the same as the accepted values.
try:
    node_id = database.retrieve_node_id(1, "UVES-Damiani")

except UnknownNodeError:
    None

else:
    database.execute(
        "UPDATE results SET passed_quality_control = false WHERE node_id = %s", 
        (node_id, ))

"""
# Remove all Elena/Carmela results because they are +/- (50*n) K multiples offset
# from the accepted values.
try:
    node_id = database.retrieve_node_id(1, "UVES-Carmela-Elena")

except UnknownNodeError:
    None

else:
    database.execute(
        "UPDATE results SET passed_quality_control = false WHERE node_id = %s",
        (node_id, ))
"""

# Remove superfluous nodes.
contributing_node_ids = database.retrieve_table(
    """ SELECT DISTINCT ON (node_id) node_id
        FROM results
        WHERE (teff <> 'NaN' AND passed_quality_control) 
    """)["node_id"]

all_node_ids = database.retrieve_table("SELECT id FROM nodes")["id"]
superfluous_node_ids = list(set(all_node_ids).difference(contributing_node_ids))
for node_id in superfluous_node_ids:
    database.execute("DELETE FROM nodes WHERE id = %s", (node_id, ))

logger.info("Removed {} superfluous nodes".format(len(superfluous_node_ids)))


# The Photometry node have 65 instances where the SETUP entry does not match the
# FILENAME entry.
uves_node_id = database.retrieve_node_id(1, "UVES-Photometry")
giraffe_node_id = database.retrieve_node_id(1, "GIRAFFE-Photometry")
database.execute(
    """ UPDATE  results
        SET     setup = 'GIRAFFE',
                node_id = '{giraffe_node_id}'
        WHERE   filename LIKE 'gir_%'
          AND   setup LIKE 'UVES%'
          AND   node_id = '{uves_node_id}';
    """.format(giraffe_node_id=giraffe_node_id, uves_node_id=uves_node_id))

database.execute(
    """ UPDATE  results
        SET     setup = 'UVES',
                node_id = '{uves_node_id}'
        WHERE   filename LIKE 'u%'
          AND   setup LIKE 'GIRAFFE%'
          AND   node_id = '{giraffe_node_id}';
    """.format(giraffe_node_id=giraffe_node_id, uves_node_id=uves_node_id))


# Identify any spurious results from the same spectrum.
either, both = (3, 2.5) # sigma thresholds

for wg in (1, ):

    results = database.retrieve_table(
        """ SELECT id, cname, node_id, teff, feh
            FROM   results 
            WHERE  passed_quality_control;""")
    if results is None: continue

    results = results.group_by("cname")
    
    spurious_results = []
    for group in results.groups:

        teff_sigmas = np.abs((np.nanmedian(group["teff"]) - group["teff"])/np.nanstd(group["teff"]))
        feh_sigmas = np.abs((np.nanmedian(group["feh"]) - group["feh"])/np.nanstd(group["feh"]))

        both_match = ((teff_sigmas > both) * (feh_sigmas > both))
        either_match = ((teff_sigmas > either) + (feh_sigmas > either))
        joint_match = both_match + either_match

        if any(joint_match):
            spurious_results.extend(group["id"][joint_match])

    if any(spurious_results):
        database.execute(
            """ UPDATE  results
                SET     passed_quality_control = false
                WHERE   id IN %s
            """, (tuple(spurious_results), ))

database.connection.commit()
