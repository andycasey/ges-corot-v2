
""" Make diagnostic plots from the UVES homogenisation models. """

import yaml
import logging
import numpy as np
import os
from astropy.table import Table
from collections import OrderedDict


from code import GESDatabase, plot
from code.model.ensemble import EnsembleModel

PLOT_DIR = "figures/homogenisation/uves/"

# Initialize logging.
logger = logging.getLogger("ges_corot")

# Create a database object.
db_filename = "db.yaml"
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)
database = GESDatabase(**credentials)

wg = 1
parameters = ("feh", "teff", "feh")

model_paths = "homogenisation-uves-wg{wg}-{parameter}.model"

plot_kwds = dict(format_node_name=lambda name: name[5:])

minimum_node_uncertainty_kwds = {
    "teff": dict(vmin=100, vmax=1000),
    "feh": dict(vmin=None, vmax=None)
}

path_kwds = dict(prefix="ges-corot-uves", wg=1)
for parameter in parameters:

    path_kwds["parameter"] = parameter

    model_path = model_paths.format(wg=wg, parameter=parameter)

    model = EnsembleModel.read(model_path, database)

    # Plot the correlation coefficients between different nodes.
    basename = "{prefix}-{wg}-{parameter}-correlations".format(**path_kwds)
    fig_correlations = model.plot_node_correlations(**plot_kwds)
    fig_correlations.savefig(os.path.join(PLOT_DIR, "{}.png".format(basename)))
    fig_correlations.savefig(os.path.join(PLOT_DIR, "{}.pdf".format(basename)))

    # Plot the biases of each node.
    basename = "{prefix}-{wg}-{parameter}-biases".format(**path_kwds)
    fig_biases = model.plot_biases(**plot_kwds)
    fig_biases.savefig(os.path.join(PLOT_DIR, "{}.png".format(basename)))
    fig_biases.savefig(os.path.join(PLOT_DIR, "{}.pdf".format(basename)))


    # Plot the systematic error across parameter space.
    kwds = dict(
        slicer=np.mean, colorbar_label="{node_name} mean uncertainty in {parameter}")
    kwds.update(plot_kwds)
    kwds.update(minimum_node_uncertainty_kwds.get(parameter, {}))

    for node_name in model._metadata["node_names"]:
        path_kwds["node"] = node_name
        basename = "{prefix}-{wg}-{parameter}-{node}-uncertainty".format(**path_kwds)
        
        fig_uncertainty = model.plot_minimum_node_uncertainty(node_name, **kwds)
        fig_uncertainty.savefig(os.path.join(PLOT_DIR, "{}.png".format(basename)))
        fig_uncertainty.savefig(os.path.join(PLOT_DIR, "{}.pdf".format(basename)))
