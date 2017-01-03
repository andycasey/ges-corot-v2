
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

model_paths = "homogenisation-uves-{parameter}.model"

wgs = (1, )
parameter_scales = OrderedDict([
    ("teff", 250),
    ("feh", 0.10)
])

sample_kwds = dict(chains=4, iter=10000)

finite = np.isfinite(benchmarks["TEFF"] * benchmarks["FEH"])
benchmarks = benchmarks[finite]

models = {}
for wg in wgs:
    models[wg] = {}
    for parameter, scale in parameter_scales.items():

        model_path = model_paths.format(wg=wg, parameter=parameter)
        
        if os.path.exists(model_path): 
            model = EnsembleModel.read(model_path, database)
            
        else:
            model = EnsembleModel(database, wg, parameter, benchmarks)
            data, metadata = model._prepare_data(
                default_sigma_calibrator=scale, sql_constraint="n.id not in (10, 11, 13, 14, 15)")
            assert all([n.startswith("UVES-") for n in metadata["node_names"]])


            init = {
                "truths": data["mu_calibrator"],
                "biases": np.zeros(data["N"]),
                "missing_estimates": np.random.uniform(
                    data["lower_bound"], data["upper_bound"], size=data["TM"]),
                "sigma_sys_constant": scale * np.ones(data["N"]),
                "L_corr": np.eye(data["N"])
            }

            init.update({"vs_tb{}".format(i): 1 for i in range(1, 9)})
            init.update({"vs_lb{}".format(i): 1 for i in range(1, 9)})
            init.update({"vs_fb{}".format(i): 1 for i in range(1, 9)})

            init.update({"vs_ta{}".format(i): 2 for i in range(1, 9)})
            init.update({"vs_la{}".format(i): 2 for i in range(1, 9)})
            init.update({"vs_fa{}".format(i): 2 for i in range(1, 9)})
            
            init["vs_tc7"] = 1 * np.ones(data["N"])
            init["vs_tc8"] = 1 * np.ones(data["N"])
            init["vs_tc9"] = 1 * np.ones(data["N"])

            print("Number of model parameters for {}: {}".format(parameter,
                sum([np.array(v).size for v in init.values()])))

            op_params = model.optimize(data, init=init, iter=100000)
            fit = model.sample(data, init=op_params, **sample_kwds)

            model.write(model_path, overwrite=True, 
                __ignore_model_pars=("Sigma", "full_rank_estimates"))



        #model.homogenise_all_stars(update_database=True)
        #model.homogenise_stars_matching_query(
        #    """SELECT distinct on (r.cname) r.cname FROM results as r, spectra as s, nodes as n where n.id = r.node_id and s.cname = r.cname and (
        #        s.ges_type like '%_OC%' or s.ges_type like '%_GC%' or s.ges_type like '%_CL%') and n.wg = 11""")
        #model.homogenise_stars_matching_query(
        #    """SELECT distinct on (r.cname) r.cname FROM results as r, spectra as s, nodes as n where n.wg = 11 and s.cname = r.cname and s.ges_fld like 'Rup134%'""")
        
