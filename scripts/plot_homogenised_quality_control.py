
""" Make some quality control plots for the homogenised results. """

import yaml
import logging
import numpy as np
import os
from astropy.table import Table
from collections import OrderedDict

from code import GESDatabase, plot
from code.model.ensemble import EnsembleModel


PLOT_DIR = "figures/homogenisation/"

# Initialize logging.
logger = logging.getLogger("ges_corot")

# Create a database object.
db_filename = "db.yaml"
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)
database = GESDatabase(**credentials)

prefix, wg = ("ges-corot", 1)

# Load the "benchmarks"
benchmarks = Table.read("fits-templates/benchmarks/GES_iDR5_FGKM_Benchmarks_ARC_29092016.fits")
# Only use "benchmarks" with TEFF < 8000 K
benchmarks = benchmarks[benchmarks["TEFF"] < 8000]


# Plot the benchmarks.
basename = "{prefix}-wg{wg}-benchmarks".format(prefix=prefix, wg=wg)
fig = plot.wg_benchmark_performance(database, wg, benchmarks, 
    skip_missing=False, show_recommended=True)
fig.savefig(os.path.join(PLOT_DIR, "{}.png".format(basename)))
fig.savefig(os.path.join(PLOT_DIR, "{}.pdf".format(basename)))


# Make a H-R diagram.
