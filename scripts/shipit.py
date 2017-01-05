
""" Ship a file containing homogenised GES/CoRoT results. """

import yaml
import logging

from code import GESDatabase, ship

# Initialize logging.
logger = logging.getLogger("ges")

# Create a database object.
db_filename = "db.yaml"
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)
database = GESDatabase(**credentials)

# Produce a homogenised file.
ship.homogenised_catalog(database,
    "fits-templates/masterlist/MasterGES_CoRoT_14Oct2016.fits",
    "outputs/ges-corot-homogenised.fits",
    wg=1, overwrite=True)
