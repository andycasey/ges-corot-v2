
"""
Homogenisation models.
"""

import yaml
import logging
import numpy as np
import os
from astropy.table import Table
from collections import OrderedDict


from code import GESDatabase
from code.model.ensemble import EnsembleModel, MedianModel

# Initialize logging.
logger = logging.getLogger("ges_corot")

# Create a database object.
db_filename = "db.yaml"
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)
database = GESDatabase(**credentials)

# Load the "benchmarks"
# Only use "benchmarks" with TEFF < 8000 K
benchmarks = Table.read("fits-templates/benchmarks/GES_iDR5_FGKM_Benchmarks_ARC_29092016.fits")
benchmarks = benchmarks[benchmarks["TEFF"] < 8000]
benchmarks["E_FEH"] = 0.05

model_paths = "homogenisation-uves-wg{wg}-{parameter}.model"

wgs = (1, )
parameter_scales = OrderedDict([
    ("teff", 250.0),
    ("feh", 0.10)
])

lower_sigma = dict(teff=10, feh=0.01)

sample_kwds = dict(chains=4, iter=20000)

finite = np.isfinite(benchmarks["TEFF"] * benchmarks["LOGG"] * benchmarks["FEH"])
benchmarks = benchmarks[finite]

models = {}
for wg in wgs:
    models[wg] = {}
    for parameter, scale in parameter_scales.items():

        model_path = model_paths.format(wg=wg, parameter=parameter)
        
        if os.path.exists(model_path): 
            model = EnsembleModel.read(model_path, database)
            
        else:
            model = EnsembleModel(database, wg, parameter, benchmarks,
                model_path="code/model/ensemble-model-5node.stan")
            data, metadata = model._prepare_data(
                default_sigma_calibrator=scale, 
                sql_constraint="n.name like 'UVES-%'")
            assert all([n.startswith("UVES-") for n in metadata["node_names"]])
            
            data["lower_sigma"] = lower_sigma[parameter]
            init = {
                "truths": data["mu_calibrator"],
                "biases": np.zeros(data["N"]),
                "missing_estimates": np.random.uniform(
                    data["lower_bound"], data["upper_bound"], size=data["TM"]),
                "sigma_sys_constant": scale * np.ones(data["N"]),
                "L_corr": np.eye(data["N"])
            }

            # Polynomial coefficients
            for i in range(1, 9):
                for tlf in "tlf":
                    for j, ba in enumerate("ba", start=1):
                        init["vs_{tlf}{ba}{i}".format(i=i, tlf=tlf, ba=ba)] = j

            for i in range(7, 10):
                init["vs_tc{}".format(i)] = np.ones(data["N"])

            # Optimize and sample.
            op_params = model.optimize(data, init=init, iter=100000)
            fit = model.sample(data, init=op_params, **sample_kwds)

            model.write(model_path, overwrite=True)


        model.homogenise_stars_matching_query(
            "SELECT DISTINCT ON (cname) cname FROM results WHERE setup LIKE 'UVES%'",
            sql_constraint="setup like 'UVES%'")

