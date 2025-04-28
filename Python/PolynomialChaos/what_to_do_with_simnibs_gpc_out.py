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
from simnibs.simulation import gpc

import logging

parser = argparse.ArgumentParser(description='Exploratory; how to put the output of SimNIBS-builtin UQ to use for a specific question')
parser.add_argument('-f', '--fn_subject', help='path of subject', required=True, type=str)
parser.add_argument('-e', '--exp_id', help='id of the experiment', required=True, type=str)
parser.add_argument('-n', '--n_cpu', help='How many cpus to use', type=int, default=21)
args = parser.parse_args()

fn_subject = os.path.abspath(args.fn_subject)
subject_id = os.path.split(fn_subject)[1].split(".")[0]
exp_id = args.exp_id
mesh_id = "mesh0"
hemisphere = exp_id.split('-')[1][0] # map-R2L/map-R -> ['map', 'R...'] -> 'R...' -> 'R'
roi_id = f"midlayer_{hemisphere.lower()}"

subject = pynibs.load_subject(fn_subject)
fn_raw_data_table = subject.exp[exp_id]['response_csv']
raw_data = pd.read_csv(os.path.join(subject.subject_folder, fn_raw_data_table))

mesh_simn = mesh_io.read_msh(subject.mesh[mesh_id]['fn_mesh_msh'])
roi = pynibs.load_roi_surface_obj_from_hdf5(subject.mesh[mesh_id]['fn_mesh_hdf5'])[roi_id]

from simnibs import Msh

print(f"n tris in mesh = {np.sum(mesh_simn.elm.elm_type == 2)}")
print(f"n tets in mesh = {np.sum(mesh_simn.elm.elm_type == 4)}")
print(f"n nodes in mesh = {mesh_simn.nodes.nr}")
print(f"tags of elements = {mesh_simn.elm.tag1}")
print(f"# elements with tag == 1 = {np.sum(mesh_simn.elm.tag1 == 1)}")
print(f"# elements with tag == 2 = {np.sum(mesh_simn.elm.tag1 == 2)}")
print(f"                together = {np.sum(mesh_simn.elm.tag1 == 1) + np.sum(mesh_simn.elm.tag1 == 2)}")

cropped = mesh_simn.crop_mesh(tags=[1, 2])


i_trial = 1



fn_hdf5 = f'tms_gpc/TMS_UQ/{subject_id}/exp_{exp_id}/e_{i_trial:04.0f}_gpc.hdf5'
# Read the regression object from the HDF5 file
regression = gpc.gPC_regression.read_hdf5(fn_hdf5)


# Example from doc:
def percentile_99(Es):
    # The function will receive the electric field in a format
    # N_sims x N_elm x 3
    # for each simulation, we calculate the 99th percentile
    prc = np.zeros(Es.shape[0])
    for i, E in enumerate(Es):
        prc[i] = simnibs.ElementData(E, mesh=mesh_simn).get_percentiles(99)[0]

    return prc


# This works just fine!
def just_select_a_few(Es):
    result = np.zeros((Es.shape[0], 5))
    for i, E in enumerate(Es): # E is N_elm x 3
        result[i,:] = np.sqrt(np.sum(np.power(E[:5, :], 2), axis=1))
    return result
# gpc_coeffs = regression.expand_quantity(just_select_a_few)


def get_midlayer_E_mag(Es):
    E_mag = np.zeros((Es.shape[0], roi.n_tris))
    for i, E in enumerate(Es):
        res = calc_e_in_midlayer_roi(E, roi=roi, mesh=cropped, qoi=["mag"], layer_gm_wm_info=None)
        assert(len(res.shape) < 2 or res.shape[1] < 2)
        E_mag[i,:] = res.flatten()
    return E_mag


f_in_question = get_midlayer_E_mag
# Calculate the gPC coefficients
gpc_coeffs = regression.expand_quantity(f_in_question)  # I.e. i can obtain the uncertainty of the electric field



# Could i do:
# Electric field across whole mesh in one simulation (of multiple for each trial)
#  -> Compute E_mag in the region of interest in this one simulation (of multiple for each trial)
#  (maybe: separately do this for each element of ROI:)
#  -> Obtain Distribution of E_mag in each element of ROI in this trial: Can sample from distribution of E_mag in each trial in each element
#  -> Set up enveloping gPC-problem that does the regression and samples from these distributions: Only needs to do the regression step, not E-field simulation!

# Question (once i'm done with this): Can this actually be done also for the uncertainty of the coil position in each trial?


print(f"Function: {f_in_question}")
print("Mean Value: ", regression.mean(gpc_coeffs).shape)
print("Standard Deviation: ", regression.std(gpc_coeffs).shape)

# Draw 10 samples for the specified function
# samples_in will be (n_samples, n_uncertain_variables) -- here, n_uncertain_variables = 4 (WM, GM, CSF, compact Bone conductivities)
samples_in, samples_out = regression.MC_sampling(gpc_coeffs, 50)
print(f"Samples in:  shape={samples_in.shape}, min={np.min(samples_in, axis=0)}, max={np.max(samples_in, axis=0)}")
print(f"Samples out: {samples_out.shape}")

