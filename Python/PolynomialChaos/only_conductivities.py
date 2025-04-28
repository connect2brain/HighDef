# Since confuctivities AND coil location is not computable (too heavy), opt for the established: only conductivity uncertainty

import argparse, os, shutil
from datetime import datetime
import numpy as np
import pandas as pd

import h5py
import pynibs
from pynibs.util.simnibs import calc_e_in_midlayer_roi
from pynibs.regression import regress_data

from collections import OrderedDict
from functools import partial

from simnibs import opt_struct, mesh_io
from simnibs.optimization import opt_struct
from simnibs.simulation import fem, sim_struct, save_matlab_sim_struct
from simnibs.simulation import gpc

import logging



parser = argparse.ArgumentParser(description='Runs generalized polynomial chaos expansion for Uncertainty quantification')
parser.add_argument('-f', '--fn_subject', help='path of subject', required=True, type=str)
parser.add_argument('-e', '--exp_id', help='id of the experiment', required=True, type=str)
parser.add_argument('-s', '--start_trial', help='Which trial to start from; if given, data may be overwritten', type=int, default=None)
parser.add_argument('-p', '--stop_trial', help='Which trial to stop (exclusive; start:stop); if given, data may be overwritten -- only when start_trial and stop_trial are both unset (default) will any preexisting data be moved to a folder labeled as outdated', type=int, default=None)
args = parser.parse_args()



logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
ROOT = "/mnt/d/HighDef-operate"
file_handler = logging.FileHandler(f"{ROOT}/TMS_UQ/gpc_only_cond_run-{datetime.now().strftime('%Y-%m-%dT%H%M%S')}.log", encoding='utf-8')
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

logger.info(f'Generalized polynomial chaos expansion for:')
logger.info(f'subject = {subject_id}, exp_id = {exp_id}, mesh_id = {mesh_id}, roi_id = {roi_id}')
logger.info(f'from trial {args.start_trial if args.start_trial is not None else 0} to trial {args.stop_trial if args.stop_trial is not None else "last"}')



config = {"fn_subject": fn_subject, "exp_id": exp_id, "mesh_id": mesh_id, "hemisphere": hemisphere, "roi_id": roi_id, "n_cpu": 1, "timestamp": datetime.now().strftime("%Y%m%dT%H%M%S")}

subject = pynibs.load_subject(fn_subject)
fn_raw_data_table = subject.exp[exp_id]['response_csv']
raw_data = pd.read_csv(os.path.join(subject.subject_folder, fn_raw_data_table))


# LIMIT TRIALS FOR NOW!
if args.start_trial is not None:
    if args.stop_trial is not None:
        raw_data = raw_data.iloc[args.start_trial:args.stop_trial,]
    else:
        raw_data = raw_data.iloc[args.start_trial:,]
elif args.stop_trial is not None: # i.e. start_trial is unset, but stop_trial is set
    raw_data = raw_data.iloc[0:args.stop_trial,]
# raw_data = raw_data.iloc[:10,]




coil_transforms = np.tile(np.reshape(np.eye(4), (4,4,1)), (1,1,len(raw_data)))
for j,axis in enumerate(['x', 'y', 'z', 'p']):
    for i,component in enumerate([f'{axis}{u}' for u in range(1,4)]):
        print(f'Retrieving {component} into {i},{j}')
        coil_transforms[i,j,:] = raw_data[component]

orientation = np.unique(raw_data['Coordinate_space'])
if len(orientation) > 1:
    NotImplementedError('The coordinate system switched during the recording!')
orientation = orientation[0] # Unpack the singular element
print(f'\nCoordinate system = {orientation}')

if orientation == "LPS":
    print('Transforming from LPS to RAS coordinates!')
    # flip first two axes
    coil_transforms[0,:,:] *= -1 # L -> R
    coil_transforms[1,:,:] *= -1 # P -> A
    orientation = "RAS"

subject.exp[exp_id]["fn_exp_hdf5"] = [os.path.join(subject.subject_folder, "exp", "UQ", exp_id, mesh_id, "experiment.hdf5")]
temp_dir = os.path.join(os.path.split(subject.exp[exp_id]['fn_exp_hdf5'][0])[0], "nnav2simnibs", mesh_id)
fn_nifti = subject.mri[0]["fn_mri_T1"]
fn_nifti_conform = os.path.join(subject.mesh[mesh_id]["mesh_folder"], subject.mesh[mesh_id]["fn_mri_conform"])
coil_transforms_simnibs = pynibs.expio.nnav2simnibs(fn_exp_nii=fn_nifti, 
                                                    fn_conform_nii=fn_nifti_conform, 
                                                    m_nnav=coil_transforms, 
                                                    nnav_system='Localite', 
                                                    orientation=orientation, 
                                                    target='simnibs', 
                                                    temp_dir=temp_dir)

pos_matrices = [coil_transforms_simnibs[:, :, i] for i in range(coil_transforms_simnibs.shape[-1])]
n_trials = len(pos_matrices)

mesh_simn = mesh_io.read_msh(subject.mesh[mesh_id]['fn_mesh_msh'])
roi = pynibs.load_roi_surface_obj_from_hdf5(subject.mesh[mesh_id]['fn_mesh_hdf5'])[roi_id]

postpro = partial(
        calc_e_in_midlayer_roi,
        roi=roi,
        mesh=mesh_simn,
        qoi=["E", "mag"],
        layer_gm_wm_info=None,
)

didt_list = [0.96*i for i in raw_data['Intensity_percentMSO']]



## Define the TMS simulation
tms = sim_struct.TMSLIST()
tms.fnamecoil = "/mnt/c/Users/bnplab-admin/SimNIBS-4.0/simnibs_env/Lib/site-packages/simnibs/resources/coil_models/Drakaki_BrainStim_2022/MagVenture_Cool-B35.ccd"
tms.mesh = mesh_simn
for i, M in enumerate(pos_matrices):
    pos = tms.add_position()
    pos.matsimnibs = M
    pos.didt = didt_list[i]

#tms.postprocess = postpro tms.postprocess has to be a string or list of strings from ["E", "e", "J", "j"]

# Define the uncertain conductivities
tms.cond[0].distribution_type = 'beta'
tms.cond[0].distribution_parameters = [3, 3, .1, .4] # WM
tms.cond[1].distribution_type = 'beta'
tms.cond[1].distribution_parameters = [3, 3, .1, .6] # GM
tms.cond[2].distribution_type = 'beta'
tms.cond[2].distribution_parameters = [3, 3, 1.2, 1.8] # CSF
tms.cond[6].distribution_type = 'beta'
tms.cond[6].distribution_parameters = [3, 3, 0.003, 0.012] # compact bone


# Run the UQ calling with White and Gray matter as an ROI and tolerance of 1e-2
pathname_gpc = f'{ROOT}/TMS_UQ/raw/{subject_id}/exp_{exp_id}/'
if args.start_trial is not None or args.stop_trial is not None:
    pathname_gpc = f'{pathname_gpc}trials_{args.start_trial if args.start_trial is not None else "first"}-{args.stop_trial if args.stop_trial is not None else "last"}/'

if not os.path.exists(pathname_gpc):
    os.makedirs(pathname_gpc)
else:
    if args.start_trial is None and args.stop_trial is None:
        shutil.move(pathname_gpc, f'{ROOT}/TMS_UQ/raw/outdated/{datetime.now().strftime("%Y%m%dT%H%M%S")}/{subject_id}/exp_{exp_id}/')
        os.makedirs(pathname_gpc)
    else:
        logger.info(f"Since starttrial andor stoptrial were given (from {args.start_trial if args.start_trial is not None else 'first'} to {args.stop_trial if args.stop_trial is not None else 'last'}), will NOT move current data --- this may overwrite some data!")


# Add the following changes manually:
# ll. 512 -- 513
# def run_tms_gpc(poslist, fn_simu, cpus=1, tissues=[2], eps=1e-2,
#                max_iter=1000, min_iter=2, data_poly_ratio=2, offset_trial=1): # <-- offset_trial added myself
# AND:
# l. 556: fn_hdf5 = fn_simu+'_{0:0=4d}_gpc.hdf5'.format(i + offset_trial)
gpc.run_tms_gpc(tms, f"{pathname_gpc}e", eps=0.5, tissues=[1, 2])


"""


#Yes, postprocessing. I specifically have to run a regression to two external responses (A, B) for each element in my region of interest.
for i_trial in range(len(raw_data['Intensity_percentMSO'])):
    fn_hdf5 = f'tms_gpc/TMS_UQ/{subject_id}_{exp_id}_{i_trial:04.0f}_gpc.hdf5'
    # Read the regression object from the HDF5 file
    regression = gpc.gPC_regression.read_hdf5(fn_hdf5)
    regression.MC_sampling"
"""