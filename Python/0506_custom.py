
from functools import partial
import pandas as pd
import numpy as np
import pynibs
from pynibs.util.simnibs import calc_e_in_midlayer_roi
import os
import h5py
from simnibs import opt_struct, mesh_io
from simnibs.optimization import opt_struct
from simnibs.simulation import fem, sim_struct, save_matlab_sim_struct
from time import sleep
import datetime
import argparse



parser = argparse.ArgumentParser(description='Creates a new subject in path')
parser.add_argument('-f', '--fn_subject', help='path of new subject', required=True, type=str)
parser.add_argument('-e', '--exp_id', help='id of the experiment', required=True, type=str)
parser.add_argument('-m', '--mesh_id', help='Mesh ID', required=True, type=str)
parser.add_argument('-r', '--roi_id', help='ROI ID', required=True, type=str)
parser.add_argument('-n', '--n_cpu', help='How many cpus to use', type=int, default=21)
args = parser.parse_args()

fn_subject = os.path.abspath(args.fn_subject)
subject_id = os.path.split(fn_subject)[1]
exp_id = args.exp_id
mesh_id = args.mesh_id
roi_id = args.roi_id
n_cpu = args.n_cpu

# Before:
#fn_subject = 'TMS_localization/TMS_loc_results/sub-001/sub-001.hdf5'
#exp_id = 'main'
#mesh_id = 'mesh0'
#roi_id = "midlayer_larger"
#fn_raw_data_table = "sub-001_raw.csv"
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##


subject = pynibs.load_subject(fn_subject)
fn_raw_data_table = subject.exp[exp_id]['response_csv']

print(f'\nE-field simulation for {subject_id}\n   subject-file: {fn_subject}\n  response-file: {fn_raw_data_table}\n     experiment: {exp_id}\n           mesh: {mesh_id}\n            ROI: {roi_id}\n        Columns: {", ".join(subject.exp[exp_id]["response_columns"])}\n\n')

raw_data = pd.read_csv(os.path.join(subject.subject_folder, fn_raw_data_table))

# This is the "result file"
subject.exp[exp_id]["fn_exp_hdf5"] = [os.path.join(subject.subject_folder, "exp", exp_id, mesh_id, "experiment.hdf5")]
# Create path if necessary:
if not os.path.exists(os.path.split(subject.exp[exp_id]["fn_exp_hdf5"][0])[0]):
    os.makedirs(os.path.split(subject.exp[exp_id]["fn_exp_hdf5"][0])[0])




## In original script 05:
# As i understand, this makes modifications to the subject-object
"""
if subject.exp[exp_id]['nnav_system'].lower() == "localite":
    pynibs.expio.localite.merge_exp_data_localite(subject=subject,                                                      <- The pynibs subject object (~hdf5 file)
                                                  exp_idx=exp_id,                                                       <- (see above)
                                                  mesh_idx=mesh_idx,                                                    <- (see above)
                                                  coil_outlier_corr_cond=coil_outlier_corr,                             <- Boolean flag
                                                  remove_coil_skin_distance_outlier=remove_coil_skin_distance_outlier,  <- Boolean flag
                                                  coil_distance_corr=coil_distance_corr,                                <- Boolean flag
                                                  verbose=verbose,                                                      <- Boolean flag
                                                  drop_mep_idx=drop_mep_idx,                                            <- None, or a list of indices?
                                                  mep_onsets=mep_onsets,                                                <- None, or a list of indices/times?
                                                  channels=channels,                                                    <- EMG channels
                                                  cfs_data_column=cfs_data_column,                                      <- Sth related to this file-type CFS (where the data is stored)
                                                  plot=plot,                                                            <- Boolean flag
                                                  start_mep=start_mep,                                                  <- Start of time frame after TMS pulse where p2p value is evaluated (in ms); [default=18 ms]
                                                  end_mep=end_mep)                                                      <- End of time frame after TMS pulse where p2p value is evaluated (in ms); [default=35 ms]
INSIDE THIS FUNCTION:
    mep_paths_lst = subject.exp[exp_idx]['fn_data'] <- e.g. [[subject_folder + "/exp/m1/mep/mep.mat"]
    tms_paths_lst = subject.exp[exp_idx]['fn_tms_nav']
    ...
    calls:      pynibs.combine_nnav_mep
    and later:  pynibs.coil_outlier_correction_cond
                pynibs.coil_distance_correction
            
Where does the data go?
dict_lst, then into results_dct
coil positions are put into this as "coil_0", "coil_1" etc; intensity as "current"
Thence:
df_stim_data -- which does not include information about the meps

df_phys_data_postproc_emg is put into the experiment hdf5 (fn_exp_hdf5) under phys_data/postproc/EMG
Yes: In the experiment.hdf5, we find phys_data/postproc/EMG/block0_values to contain the p2p-amplitudes!
The raw EMG is saved into the same experiment file under phys_data/raw/EMG
Meanwhile, the "stim_data" is saved in the same file under stim_data
    -> "stim_data"~>"condition" seems to be just the pulse index (?)
    -> "stim_data"~>"coil_0" etc are the coil positions

Anyways, this all seems to get stored in the exp/m1/mesh_0/experiment.hdf5 file (subject.exp[exp_id]["fn_exp_hdf5"] above)

Localite coordinates are getting transformed by call to nnav2simnibs (exp.py line 1062) during call 
to combine_nnav_mep (called in localite.py line 672)

From exp.py (line 1061 ff):
    m_simnibs = np.moveaxis(coil_array[idx, :, :, :], 0, 2)
    m_simnibs = nnav2simnibs(fn_exp_nii=nii_exp_path[0],
                                fn_conform_nii=nii_conform_path,
                                m_nnav=m_simnibs,
                                nnav_system=nnav_system,
                                mesh_approach=mesh_approach,
                                temp_dir=temp_dir)

What does the axis-moving do? coil_array has shape (3={coil0,coil1,mean?}, n_trials, 4, 4)
The call np.moveaxis(coil_array[idx, :, :, :], 0, 2) -~-> shape=(4,4,n_trials), which is needed for nnav2simnibs

Then nnav2simnibs is called (see calls below), then, the result is reshaped to be (N, 4, 4) again, and later turned 
into entries coil0_00, coil0_01 etc.
Then, combine_nnav_mep returns a list of dictionaries, one per trial, with keys coil0_00 etc.

"""

fn_nifti = subject.mri[0]["fn_mri_T1"]
fn_nifti_conform = os.path.join(subject.mesh[mesh_id]["mesh_folder"], subject.mesh[mesh_id]["fn_mri_conform"])

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

temp_dir = os.path.join(os.path.split(subject.exp[exp_id]['fn_exp_hdf5'][0])[0], "nnav2simnibs", mesh_id)

coil_transforms_simnibs = pynibs.expio.nnav2simnibs(fn_exp_nii=fn_nifti, 
                                                    fn_conform_nii=fn_nifti_conform, 
                                                    m_nnav=coil_transforms, 
                                                    nnav_system='Localite', 
                                                    orientation=orientation, 
                                                    target='simnibs', 
                                                    temp_dir=temp_dir)
# coil_transforms_simnibs is now in shape (4,4,num_trials)


## Again from script_05, adapted:
fn_coil_pos = "plot_coil_pos.hdf5"
fn_coil_pos = os.path.join(os.path.split(subject.exp[exp_id]["fn_exp_hdf5"][0])[0], fn_coil_pos)

# Plot coil positions:
pynibs.create_stimsite_from_matsimnibs(fn_hdf=fn_coil_pos, matsimnibs=coil_transforms_simnibs, datanames=subject.exp[exp_id]['response_columns'], data=np.reshape(raw_data[subject.exp[exp_id]['response_columns']], (len(raw_data), len(subject.exp[exp_id]['response_columns']))), overwrite=True)

fn_matsimnibs = "matsimnibs.hdf5"
fn_matsimnibs = os.path.join(os.path.split(subject.exp[exp_id]["fn_exp_hdf5"][0])[0], fn_matsimnibs)

print(f"Writing: {fn_matsimnibs}")
with h5py.File(fn_matsimnibs, "w") as f:
    f.create_dataset("matsimnibs", data=coil_transforms_simnibs)


"""
Dissecting script 06




"""
# output folder
fn_out = os.path.join(os.path.split(fn_subject)[0], 'results',
                      f"exp_{exp_id}", "electric_field",
                      f"mesh_{mesh_id}",
                      f"roi_{roi_id}")
if not os.path.exists(fn_out):
    os.makedirs(fn_out)
# file containing the head model (.msh format)
fn_mesh_msh = subject.mesh[mesh_id]["fn_mesh_msh"]
# file containing the FEM options for SimNIBS (will be created)
fn_fem_options = os.path.join(fn_out, 'FEM_config.mat')
fn_result_suffix = ""
fn_e = os.path.join(fn_out, f"e{fn_result_suffix}.hdf5")


anisotropy_type = "scalar"

tms_opt = opt_struct.TMSoptimize()
tms_opt.fnamehead = fn_mesh_msh
tms_opt.pathfem   = fn_out
tms_opt.fnamecoil = '/mnt/c/Users/bnplab-admin/SimNIBS-4.0/simnibs_env/Lib/site-packages/simnibs/resources/coil_models/Drakaki_BrainStim_2022/MagVenture_Cool-B35.ccd'
tms_opt.target = None
tms_opt.open_in_gmsh = False
tms_opt.anisotropy_type = anisotropy_type

if os.path.exists(fn_fem_options):
    print(f"Found old FEM config file {fn_fem_options} -- overwriting it!")
    sleep(1.5)
    save_matlab_sim_struct(tms_opt, fn_fem_options)

print(f" Loading ROI surface from: {subject.mesh[mesh_id]['fn_mesh_hdf5']}")
print(f" yields: {pynibs.load_roi_surface_obj_from_hdf5(subject.mesh[mesh_id]['fn_mesh_hdf5'])}")

roi = pynibs.load_roi_surface_obj_from_hdf5(subject.mesh[mesh_id]['fn_mesh_hdf5'])[roi_id]
n_tris_roi = roi.n_tris


"""
in 06_calc_e.py, the script run_simnibs is called.
cmd = f"{sys.executable} " + os.path.join(scripts_folder, "run_simnibs.py") + \
              " --folder " + fn_out + " --fn_subject " + fn_subject + " --fn_coilpos " + fn_coilpos_hdf5 + \
              " --cpus " + str(n_cpu) + " --anisotropy_type " + anisotropy_type + " -m " + mesh_idx + \
              " -r " + "'" + str(roi_idx) + "'" + " --hdf5_fn " + "'" + fn_e + "'" + \
              (" --calc_theta_grad" if args_dict["calc_theta_grad"] else "")
print(f"Running {cmd}")
os.system(cmd)
"""
# Adapted from Script "run_simnibs"
"""# import after setting environment solver options
os.environ['SIMNIBS_SOLVER_OPTIONS'] = args.solver_opt  -- the default value here is pardiso, and in the above call, it is 
not explicitly set to anything else (? so always pardiso?)
"""

qois = ['E', 'mag', 'norm', 'tan']
fields = 'E'
dataset = '/tmp'

print(f"Reading head model (msh): {subject.mesh[mesh_id]['fn_mesh_msh']}")
mesh_simn = mesh_io.read_msh(subject.mesh[mesh_id]['fn_mesh_msh'])

print(f"Reading head model (hdf5): {subject.mesh[mesh_id]['fn_mesh_hdf5']}")
mesh_pyf = pynibs.load_mesh_hdf5(subject.mesh[mesh_id]['fn_mesh_hdf5'])

print(f"Reading ROI #{roi_id}: {subject.mesh[mesh_id]['fn_mesh_hdf5']}")
roi = pynibs.load_roi_surface_obj_from_hdf5(subject.mesh[mesh_id]['fn_mesh_hdf5'])[roi_id]
fn_hdf5 = fn_e #os.path.join(fn_out, fn_e)
print(f"fn_hdf5 = {fn_hdf5}")
if os.path.exists(fn_hdf5):
    os.rename(fn_hdf5, f'{os.path.splitext(fn_hdf5)[0]}_outdated_at_{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}.hdf5')


print(f"Preparing conductivity tensors")
conductivities = sim_struct.SimuList.cond2elmdata(tms_opt, mesh=mesh_simn)

print(f"Preparing simulations for {coil_transforms_simnibs.shape[-1]} coil positions...")
pos_matrices = [coil_transforms_simnibs[:, :, i] for i in range(coil_transforms_simnibs.shape[-1])]

print(f" DOING AD-HOC conversion of MSO to DIDT! (96 A/us * % MSO)")
didt_list = [0.96*i for i in raw_data['Intensity_percentMSO']]

# From run_simnibs:
layer_gm_wm_info = None
postpro = partial(
        calc_e_in_midlayer_roi,
        roi=roi,
        mesh=mesh_simn,
        qoi=qois,
        layer_gm_wm_info=layer_gm_wm_info,
)

t = datetime.datetime.now()
print(f"\n\nStarting: fem.tms_many_simulations at {t.strftime('%Y-%m-%d %H:%M:%S')}")

fem.tms_many_simulations(
        mesh=mesh_simn,
        cond=conductivities,
        fn_coil=tms_opt.fnamecoil,
        matsimnibs_list=pos_matrices,
        didt_list=didt_list,
        fn_hdf5=fn_hdf5,
        dataset=dataset,
        post_pro=postpro,
        solver_options=tms_opt.solver_options,
        n_workers=n_cpu,
        field=fields)

print(f"Done. FEM solving took {datetime.datetime.now() - t}")



print(f"Reshaping .hdf5")

with h5py.File(fn_hdf5, 'r') as h5:
    # undo zero padding by storing just as many values as elements on layer
    num_roi_elmts = roi.node_number_list.shape[0]
    E = h5[dataset][:, :num_roi_elmts, 0:3]
    E_mag = h5[dataset][:, :num_roi_elmts, 3]
    E_norm = h5[dataset][:, :num_roi_elmts, 4]
    E_tan = h5[dataset][:, :num_roi_elmts, 5]

    qoi_idx = 6

with h5py.File(fn_hdf5.replace('.hdf5','reshaped.hdf5'), 'w') as h5:
    h5.create_dataset('/E', data=E)
    h5.create_dataset('/E_mag', data=E_mag)
    h5.create_dataset('/E_norm', data=E_norm)
    h5.create_dataset('/E_tan', data=E_tan)






# Adapted from Script 6:
with h5py.File(fn_e, "a") as f:
    print(f"Reshaping fields.")
    E = f['tmp'][:, :, 0:3]
    E_mag = f['tmp'][:, :, 3]
    E_norm = f['tmp'][:, :, 4]
    E_tan = f['tmp'][:, :, 5]
    f.create_dataset('/E', data=E)
    f.create_dataset('/E_mag', data=E_mag)
    f.create_dataset('/E_norm', data=E_norm)
    f.create_dataset('/E_tan', data=E_tan)

print(f"{' '.join(20*['^'])}\n      Done.")




