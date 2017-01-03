#!/usr/bin/python


import yaml
import logging
import os
import matplotlib.pyplot as plt
from glob import glob

from code import (GESDatabase, plot, summary)
from astropy.table import Table

# Initialize logging.
logger = logging.getLogger("ges.corot")
logger.setLevel(logging.DEBUG)

handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
    "%(asctime)s [%(levelname)-8s] %(message)s"))
logger.addHandler(handler)


# Create a database object.
db_filename = "db.yaml"
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)
database = GESDatabase(**credentials)


prefix = "ges-corot"
debug = False
wgs = (1, )
parameters = ("teff", "logg", "feh", "xi")
extents = dict(teff=(4000, 6500), logg=(0, 5), feh=(-3, 0.5))
isochrones = glob("isochrones/*.dat")

# Create node-level folders.
nodes = database.retrieve_table("SELECT id, wg, TRIM(name) as name FROM nodes")
for node in nodes:
    folder = "figures/qc/wg{}/{}".format(node["wg"], node["name"])
    if not os.path.exists(folder):
        os.mkdir(folder)



# Node-to-node comparisons within a WG
for wg in wgs:
    for parameter in parameters:

        print("Plotting node-to-node comparison of {} for wg {}".format(
            parameter, wg))
        try:
            fig = plot.compare_nodes_within_wg(database, wg, parameter)

        except:
            if debug: raise
            logger.exception(
                "Could not create compare_nodes_within_wg figure for {}/{}:"\
                .format(wg, parameter))
            continue

        else:
            basename = "figures/qc/wg{wg}/{prefix}-{wg}-{parameter}"\
                .format(prefix=prefix, wg=wg, parameter=parameter)

            fig.savefig("{}.pdf".format(basename))
            fig.savefig("{}.png".format(basename))

plt.close("all")


# Node-to-node comparisons within a WG
for wg in wgs:
    for parameter in parameters:

        print("Plotting zoomed node-to-node comparison of {} for wg {}".format(
            parameter, wg))
        try:
            fig = plot.compare_nodes_within_wg(
                    database, wg, parameter, 
                    extent=extents.get(parameter, None))

        except:
            if debug: raise
            logger.exception(
                "Could not create compare_nodes_within_wg figure for {}/{}:"\
                .format(wg, parameter))
            continue

        else:
            basename = "figures/qc/wg{wg}/{prefix}-{wg}-{parameter}-zoom"\
                .format(prefix=prefix, wg=wg, parameter=parameter)

            fig.savefig("{}.pdf".format(basename))
            fig.savefig("{}.png".format(basename))

plt.close("all")

# Node-level figures.
for node in nodes:

    print("Plotting benchmark performance for {}".format(node))
    try:
        fig = plot.node_benchmark_performance(database, node["wg"], node["name"],
            ylims=dict(teff=500, logg=0.5, mh=0.5))

    except:
        if debug: raise
        logger.exception(
            "Could not create node_benchmark_performance figure for {}/{}:"\
            .format(node["wg"], node["name"]))
        continue

    else:
        if fig is not None:
            basename = "figures/qc/wg{wg}/{name}/{prefix}-{wg}-{name}-benchmarks-zoom".format(
                wg=node["wg"], name=node["name"], prefix=prefix)

            fig.savefig("{}.pdf".format(basename))
            fig.savefig("{}.png".format(basename))

        else:
            logger.warn("No benchmark information for {}/{}".format(
                node["wg"], node["name"]))

plt.close("all")



for node in nodes:

    print("Plotting benchmark performance for {}".format(node))
    try:
        fig = plot.node_benchmark_performance(database, node["wg"], node["name"])

    except:
        if debug: raise
        logger.exception(
            "Could not create node_benchmark_performance figure for {}/{}:"\
            .format(node["wg"], node["name"]))
        continue

    else:
        if fig is not None:
            basename = "figures/qc/wg{wg}/{name}/{prefix}-{wg}-{name}-benchmarks".format(
                wg=node["wg"], name=node["name"], prefix=prefix)

            fig.savefig("{}.pdf".format(basename))
            fig.savefig("{}.png".format(basename))

        else:
            logger.warn("No benchmark information for {}/{}".format(
                node["wg"], node["name"]))

plt.close("all")



# H-R diagrams by setup.
for node in nodes:

    print("Plotting HR diagrams for node {}".format(node))

    try:
        fig = plot.hrd_by_setup(database, node["wg"], node["name"])

    except:
        if debug: raise
        logger.exception(
            "Could not create hrd_by_setup figure for {}/{}:"\
            .format(node["wg"], node["name"]))
        continue

    else:
        if fig is None: continue
        basename = "figures/qc/wg{wg}/{name}/{prefix}-{wg}-{name}-hrd-by-setup".format(
            wg=node["wg"], name=node["name"], prefix=prefix)

        fig.savefig("{}.pdf".format(basename))
        fig.savefig("{}.png".format(basename))

plt.close("all")







# Total HRD.
for node in nodes:

    print("Plotting total HRD for {}".format(node))
    try:
        fig = plot.hrd(database, node["wg"], node["name"])

    except:
        if debug: raise
        logger.exception(
            "Could not create hrd figure for {}/{}:".format(
                node["wg"], node["name"]))
        continue

    else:
        basename = "figures/qc/wg{wg}/{name}/{prefix}-{wg}-{name}-hrd".format(
            wg=node["wg"], name=node["name"], prefix=prefix)

        fig.savefig("{}.pdf".format(basename))
        fig.savefig("{}.png".format(basename))

plt.close("all")


# Plot results for the Sun.
for node in nodes:

    print("PLotting solar results for {}".format(node))

    try:
        fig = plot.hrd(database, node["wg"], node["name"], mark=(5777, 4.4),
            where="CNAME = 'ssssssss-sssssss'")

    except:
        if debug: raise
        logger.exception(
            "Could not create hrd (SUN) figure for {}/{}:".format(
                node["wg"], node["name"]))
        continue

    else:
        basename = "figures/qc/wg{wg}/{name}/{prefix}-{wg}-{name}-solar".format(
            wg=node["wg"], name=node["name"], prefix=prefix)

        fig.savefig("{}.pdf".format(basename))
        fig.savefig("{}.png".format(basename))


plt.close("all")



# Plot distributions of stellar parameters.
for node in nodes:

    print("Plotting stellar parameter distributions for {}".format(node))

    try:
        fig = plot.stellar_parameter_histograms(database, node["wg"], node["name"])

    except:
        if debug: raise
        logger.exception(
            "Could not create stellar_parameter_histograms for {}/{}:".format(
                node["wg"], node["name"]))
        continue

    else:
        basename = "figures/qc/wg{wg}/{name}/{prefix}-{wg}-{name}-histogram"\
            .format(wg=node["wg"], name=node["name"], prefix=prefix)

        fig.savefig("{}.pdf".format(basename))
        fig.savefig("{}.png".format(basename))


plt.close("all")

# Plot distributions of errors in stellar parameters.
for node in nodes:

    print("Plotting stellar parameter error distributions for {}".format(node))
    
    try:
        fig = plot.stellar_parameter_error_histograms(database, node["wg"], node["name"])

    except:
        if debug: raise
        logger.exception(
            "Could not create stellar_parameter_error_histograms for {}/{}:".format(
                node["wg"], node["name"]))
        continue

    else:
        basename = "figures/qc/wg{wg}/{name}/{prefix}-{wg}-{name}-error-histogram".format(
            wg=node["wg"], name=node["name"], prefix=prefix)

        fig.savefig("{}.pdf".format(basename))
        fig.savefig("{}.png".format(basename))

plt.close("all")



# Create summary tables.
parameter_summary = summary.stellar_parameter_summary(database)
parameter_summary.write("figures/qc/parameter-summary.txt",
    format="ascii")

parameter_ranges = summary.stellar_parameter_range(database)
parameter_ranges.write("figures/qc/parameter-range-summary.txt",
    format="ascii")

for node in nodes:

    tech = summary.tech_flags(database, node["wg"], node["name"])
    if tech is not None:
        tech.write(
            "figures/{prefix}-{wg}-{name}-tech-summary.txt".format(
                prefix=prefix, wg=node["wg"], name=node["name"]),
            format="ascii")

