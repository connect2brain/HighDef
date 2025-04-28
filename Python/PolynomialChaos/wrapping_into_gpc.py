# Follows what_to_do_with_simnibs_gpc_out.py
# Goal now: Wrap the per-trial GPC-regressions (that provide uncertainty of conductivity -> uncertainty of E-mag in all elements of ROI) into a greater gpc-model

import argparse, os
from datetime import datetime
import numpy as np
import pandas as pd

import h5py
import pynibs
from pynibs.util.simnibs import calc_e_in_midlayer_roi
from pynibs.regression import regress_data

from collections import OrderedDict
from functools import partial

import simnibs
from simnibs import opt_struct, mesh_io
from simnibs.optimization import opt_struct
from simnibs.simulation import fem, sim_struct, save_matlab_sim_struct
from simnibs.simulation import gpc, gPC_regression

import pygpc
# from pygpc.io import read_gpc_pkl, read_data_hdf5
import logging

from uq_models import WrapperModel, SubModel



def add_normalized_beta_parameter(shape):
    return pygpc.Beta(pdf_shape=shape, pdf_limits=[-1, 1])








parser = argparse.ArgumentParser(description='Runs generalized polynomial chaos expansion for Uncertainty quantification')
parser.add_argument('-f', '--fn_subject', help='path of subject', required=True, type=str)
parser.add_argument('-e', '--exp_id', help='id of the experiment', required=True, type=str)
parser.add_argument('-n', '--n_cpu', help='How many cpus to use', type=int, default=23)
args = parser.parse_args()



logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
ROOT = "/mnt/d/HighDef-operate/HighDef"
file_handler = logging.FileHandler(f"{ROOT}/wrapped_gpc_run-{datetime.now().strftime('%Y-%m-%dT%H%M%S')}.log", encoding='utf-8')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(name)s: %(message)s')
file_handler.setFormatter(formatter)
file_handler.setLevel(logging.DEBUG)
logger.addHandler(file_handler)



fn_subject = os.path.abspath(args.fn_subject)
subject_id = os.path.split(fn_subject)[1].split(".")[0]
exp_id = args.exp_id
mesh_id = "mesh0"
hemisphere = exp_id.split('-')[1][0] # map-R2L/map-R -> ['map', 'R...'] -> 'R...' -> 'R'
roi_id = f"midlayer_{hemisphere.lower()}"
n_cpu = args.n_cpu

logger.info(f'Generalized polynomial chaos expansion for:')
logger.info(f'subject = {subject_id}, exp_id = {exp_id}, mesh_id = {mesh_id}, roi_id = {roi_id}')
logger.info(f'running on {n_cpu} CPUs')

config = {"fn_subject": fn_subject, "exp_id": exp_id, "mesh_id": mesh_id, "hemisphere": hemisphere, "roi_id": roi_id, "n_cpu": n_cpu, "timestamp": datetime.now().strftime("%Y%m%dT%H%M%S")}


subject = pynibs.load_subject(fn_subject)
fn_raw_data_table = subject.exp[exp_id]['response_csv']
raw_data = pd.read_csv(os.path.join(subject.subject_folder, fn_raw_data_table))
raw_data = raw_data.iloc[:200,]

mesh_simn = mesh_io.read_msh(subject.mesh[mesh_id]['fn_mesh_msh'])
roi = pynibs.load_roi_surface_obj_from_hdf5(subject.mesh[mesh_id]['fn_mesh_hdf5'])[roi_id]
cropped = mesh_simn.crop_mesh(tags=[1, 2])







model = WrapperModel(subject, cropped, roi, raw_data, config)

parameters = OrderedDict()
parameters["c_wm"]  = add_normalized_beta_parameter([3, 3])
parameters["c_gm"]  = add_normalized_beta_parameter([3, 3])
parameters["c_csf"] = add_normalized_beta_parameter([3, 3])
parameters["c_cb"]  = add_normalized_beta_parameter([3, 3])

problem = pygpc.Problem(model, parameters)





options = {}
options["order_start"] = 2
options["order_end"] = 10
options["solver"] = "LarsLasso"
options["interaction_order"] = 2
options["order_max_norm"] = 0.7
options["n_cpu"] = 0
options["adaptive_sampling"] = True
options["gradient_enhanced"] = False
options["fn_results"] = f"/mnt/d/HighDef-operate/TMS_UQ/wrapped/{subject_id}/exp_{exp_id}/gpc" 
options["error_type"] = "loocv"
options["error_norm"] = "absolute"
options["eps"] = 1.0
options["GPU"] = True

algorithm = pygpc.RegAdaptive(problem=problem, options=options)

session = pygpc.Session(algorithm)

session, coeffs, results = session.run()

print(f"How often was the model evaluated: {model.how_often_evaluated} times")

mean = session.gpc[0].get_mean(coeffs)
print(f"mean = {mean}")
std  = session.gpc[0].get_std(coeffs)
print(f"std = {std}")