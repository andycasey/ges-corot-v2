#!/usr/bin/python

""" Create database for the GES/CoRoT homogenisation. """

import logging
import numpy as np
import psycopg2 as pg
import yaml
from glob import glob


# For fake data generation
import os
from astropy.io import fits

from code import GESDatabase


db_filename = "db.yaml"
nodes_filename = "nodes.yaml"
schema_filename = "code/schema.sql"
masterlist_filename = "fits-templates/masterlist/MasterGES_CoRoT_14Oct2016.fits"

# Connect to database.
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)
    connection = pg.connect(**credentials)

logger = logging.getLogger("ges")
logger.info("Connected to database.")

# Create the tables.
cursor = connection.cursor()
logger.info("Creating tables from {}...".format(schema_filename))
with open(schema_filename, "r") as fp:
    cursor.execute(fp.read())
cursor.close()
logger.info("Tables created.")

connection.commit()
connection.close()

# Create a database object.
database = GESDatabase(**credentials)

# Create nodes.
with open(nodes_filename, "r") as fp:
    all_nodes = yaml.load(fp)

for wg, node_names in all_nodes.items():
    for node_name in node_names:
        node_id = database.create_or_retrieve_node_id(wg, node_name)

# Ingest the masterlist of spectra.
N_ingested = database.ingest_spectra_masterlist(masterlist_filename)

# Ingest results from the nodes.
for filename in glob("node-results/*/*.fits"):
    N = database.ingest_node_results(filename, extension=1)
    logger.info("Ingested {} results from {}".format(N, filename))

database.connection.commit()

# Fix the MH/FEH issue:
database.execute(
    """ UPDATE results
           SET feh = mh
         WHERE feh = 'NaN'
           AND mh <> 'NaN';""")

# We don't have GES_FLD for these stars, and we need them.
# Set object as GES_FLD
database.execute("UPDATE spectra SET ges_fld = object")

# Note that there is an issue with the CNAMEs of UVES benchmark spectra. There
# are two CNAME entries for alf_Cen_A, two for alf_Cet, and two for GJ880.
database.execute(
    """ UPDATE spectra
           SET cname = '14392972-6049560'
         WHERE ges_fld like 'alf_Cen_A%'""")

database.execute(
    """ UPDATE spectra
           SET cname = '03021676+0405219'
         WHERE ges_fld like 'alf_Cet%'""")

database.execute(
    """ UPDATE spectra
           SET cname = '22563384+1633085'
         WHERE ges_fld like 'GJ880%'""")

# Some of the nodes did not keep the logg values from the masterlist file.
masterlist = fits.open(masterlist_filename)
for i in range(len(masterlist[1].data)):
    cname = masterlist[1].data["CNAME"][i]
    logg = masterlist[1].data["LOGG"][i]
    e_logg = masterlist[1].data["E_LOGG"][i]
    database.execute(
        """ UPDATE results 
            SET logg = %s, e_logg = %s 
            WHERE cname = %s """, (float(logg), float(e_logg), cname))

database.connection.commit()

# Remove superfluous fake nodes.
contributing_node_ids = database.retrieve_table(
    """ SELECT DISTINCT ON (node_id) node_id
        FROM results
        WHERE teff <> 'NaN'
    """)["node_id"]

all_node_ids = database.retrieve_table("SELECT id FROM nodes")["id"]
for node_id in set(all_node_ids).difference(contributing_node_ids):
    database.execute("DELETE FROM nodes WHERE id = %s", (node_id, ))

database.connection.commit()

logger.info("Ingestion complete.")

