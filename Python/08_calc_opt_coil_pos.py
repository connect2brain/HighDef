#!/usr/bin/env python

"""
Determine optimal coil position/orientation from localization experiment.
"""

import os
import re
import pandas as pd
import h5py
import pynibs
import simnibs
import warnings
import argparse
import numpy as np
from simnibs import opt_struct, mesh_io

if int(simnibs.__version__[0]) < 4:
    from simnibs.utils.nnav import write_tms_navigator_im
else:
    from simnibs.utils.nnav import localite
    loc = localite()
    write_tms_navigator_im = loc.write

# set up command line argument parser
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--fn_subject', help='Path to *.hdf5-file of subject', required=True)
parser.add_argument('-m', '--mesh_idx', help='Mesh ID', required=True, type=str)
parser.add_argument('-e', '--exp_idx', help='Experiment ID', required=True, type=str)
parser.add_argument('-n', '--n_cpu', help='How many cpus to use', type=int, default=4)
parser.add_argument('-a', '--anisotropy_type', help='Anisotropy "vn" or "scalar"', type=str, default="vn")
parser.add_argument('-p', '--pardiso', action='store_true', help="Use pardiso solver")
parser.add_argument('-t', '--target', nargs='+', required=True, help='R2 results file or x/y/z coordinates')
parser.add_argument('-l', '--label', required=False, help='label for the target optimization to be used for the storage folder (e.g. CsE_FDI_in_uV)', default='target')
parser.add_argument('-q', '--qoi', help='Electric field component ("E_mag", "E_norm", "E_tan")', required=False,
                    default='E_mag')
parser.add_argument('-c', '--fn_coil', help='Path to TMS coil', required=True)
parser.add_argument('--opt_search_radius', help='Optimization search radius', default=15, type=float) # previous default: 20 mm
parser.add_argument('--opt_search_angle', help='Optimization search angle', default=180, type=float)
parser.add_argument('--opt_angle_resolution', help='Optimization angle resolution', default=15, type=float) # previous default: 7.5°
parser.add_argument('--opt_spatial_resolution', help='Optimization spatial resolution', default=2, type=float)
parser.add_argument('--opt_smooth', help='Surface smoothing', default=100, type=float)
parser.add_argument('--opt_distance', help='Skin-coil distance', default=1, type=float)
parser.add_argument('--recompute', help='Force recompute everything', default=False, type=bool, required=False)
parser.description = 'Determine optimal coil position/orientation for cortical target.\n'
args = parser.parse_args()

# print simulation info
# ================================================================================
print("-" * 64)
print(f"{parser.description}")
args_dict = vars(args)

print('Parameters:')
for key in args_dict.keys():
    print(f"{key: >15}: {args_dict[key]}")
print("-" * 64)

fn_subject = os.path.abspath(args.fn_subject)
mesh_id = args.mesh_idx
exp_id = args.exp_idx
target = args.target
n_cpu = args.n_cpu
anisotropy_type = args.anisotropy_type
pardiso = args.pardiso
qoi = args.qoi
fn_coil = os.path.abspath(args.fn_coil)
opt_search_radius = args.opt_search_radius
opt_search_angle = args.opt_search_angle
opt_angle_resolution = args.opt_angle_resolution
opt_spatial_resolution = args.opt_spatial_resolution
opt_smooth = args.opt_smooth
opt_distance = args.opt_distance
label = args.label

fn_config = '/home/bnplab-admin/TMS_localization/config.json'
config = pd.read_json(fn_config)


if len(target) > 1:
    target = np.array(target).astype(float)
    print(f"Defining user defined target: {target}")
else:
    target = target[0]

# determine target coordinates from r2 map
if type(target) is str:
    print(f"Reading target from file: {target}")
    # read r2 data
    fn_r2_data = target
    fn_r2_geo = os.path.splitext(fn_r2_data)[0][:-4] + "geo.hdf5"

    # find element with maximum goodness-of-fit
    with h5py.File(fn_r2_data, 'r') as gof_res:
        qof = gof_res[f'data/tris/c_{qoi}'][:]
        if np.max(qof) < config.R2_rejection.only_apply_if_best_R2_is_lower_than:
            lqof = np.log(gof_res[f'data/tris/c_{qoi}'][:])
            Q3 = np.quantile(lqof, 0.75)
            iqr = Q3 - np.quantile(lqof, 0.25)
            cutoff = config.R2_rejection.Cutoff_IQR_factor*iqr + Q3

            banned_elements = lqof > cutoff
            print(f"Rejecting {np.sum(banned_elements)} outliers (Q3 = {Q3}, IQR = {iqr} => Cutoff = {cutoff})!")
            # print(np.where(banned_elements)[0]+1)
            np.savetxt(f"{os.path.splitext(fn_r2_data)[0][:-4]}rejected.csv", np.where(banned_elements)[0]+1, fmt="%i")
            qof[banned_elements] = 0
            best_elm_id = np.argmax(qof)
        else:
            print("No outlier rejection!")
            best_elm_id = np.argmax(qof)

    # get location of best element on cortex
    with h5py.File(fn_r2_geo, 'r') as gof_geo:
        nodes_ids = gof_geo['/mesh/elm/triangle_number_list'][best_elm_id]
        target = np.mean(gof_geo['/mesh/nodes/node_coord'][:][nodes_ids], axis=0)

# Subject information
# ========================================================================================================
subject_dir = os.path.split(fn_subject)[0]
subject = pynibs.load_subject(fn_subject)

# file containing the head model (.msh format)
fn_mesh_msh = subject.mesh[mesh_id]["fn_mesh_msh"]

assert os.path.exists(fn_coil), f"Coil file '{fn_coil}' is missing."
assert os.path.exists(fn_mesh_msh), f"Mesh file '{fn_mesh_msh}' is missing."

optim_dir = os.path.join(subject_dir,
                         f"opt",
                         exp_id,
                         f"{label}_{np.array2string(target, formatter={'float_kind': '{0:+0.2f}'.format})}",
                         mesh_id,
                         f"{qoi}")
if not os.path.exists(optim_dir):
    os.makedirs(optim_dir)

print(f'\ttarget:   {target}')
print(f'\theadmesh: {fn_mesh_msh}')
print(f'\toptim_dir: {optim_dir}')

fn_conform_nii = os.path.join(subject.mesh[mesh_id]["mesh_folder"], subject.mesh[mesh_id]["fn_mri_conform"])
try:
    fn_exp_nii = subject.exp[exp_id]["fn_mri_nii"][0][0]
except KeyError:
    print(f"Experiment {exp_id} not found.")
    fn_exp_nii = fn_conform_nii

# Initialize structure
# ========================================================================================================
tms_opt = opt_struct.TMSoptimize()
tms_opt.fnamehead = fn_mesh_msh
tms_opt.pathfem = optim_dir
tms_opt.fnamecoil = fn_coil
tms_opt.target = target
tms_opt.anisotropy_type = anisotropy_type

if subject.mesh[mesh_id]['fn_tensor_vn'] is not None:
    tms_opt.fname_tensor = os.path.join(subject.mesh[mesh_id]["mesh_folder"],
                                        subject.mesh[mesh_id]['fn_tensor_vn'])
else:
    tms_opt.fname_tensor = None

if pardiso:
    tms_opt.solver_options = 'pardiso'  # 'pardiso'  # faster, lots of RAM
else:
    tms_opt.solver_options = None
tms_opt.open_in_gmsh = False

tms_opt.search_radius = opt_search_radius
tms_opt.search_angle = opt_search_angle
tms_opt.angle_resolution = opt_angle_resolution
tms_opt.spatial_resolution = opt_spatial_resolution
tms_opt.smooth = opt_smooth
tms_opt.distance = opt_distance

tms_opt.method = 'ADM'

# optimization target area in mm
tms_opt.target_size = .5
mesh = mesh_io.read_msh(fn_mesh_msh)

# check if the target is within GM
# ========================================================================================================
bar = mesh.elements_baricenters()[:]
dist = np.linalg.norm(bar - target, axis=1)
elm = mesh.elm.elm_number[
    (dist < tms_opt.target_size) *
    np.isin(mesh.elm.tag1, [2]) *
    np.isin(mesh.elm.elm_type, [4])
    ]
if not len(elm):
    warnings.warn(f"\tNo elements found at {target}. Changing target.")
    sub_elms = bar[np.isin(mesh.elm.tag1, [2]) * np.isin(mesh.elm.elm_type, [4])]
    new_sub_coords = sub_elms[np.where(
        np.linalg.norm(sub_elms - target, axis=1)
        ==
        np.min(np.linalg.norm(sub_elms - target, axis=1)))[0][0]]

    tms_opt.target = new_sub_coords
    print(f'\tcoords: {new_sub_coords}')

# Run optimization
# ========================================================================================================
if not (os.path.exists(f'{optim_dir}{os.sep}opt_coil_pos.xml') and
        os.path.exists(f'{optim_dir}{os.sep}opt_coil_pos.hdf5')) or args.recompute:
    tms_opt._prepare()
    pos_matrices = tms_opt._get_coil_positions()
    tms_opt.pos_ydir = []
    tms_opt.centre = []
    pynibs.write_coil_pos_hdf5(fn_hdf=os.path.join(optim_dir, "search_positions.hdf5"),
                               centers=np.vstack([p[0:3, 3] for p in pos_matrices]),
                               m0=np.vstack([p[0:3, 0] for p in pos_matrices]),
                               m1=np.vstack([p[0:3, 1] for p in pos_matrices]),
                               m2=np.vstack([p[0:3, 2] for p in pos_matrices]),
                               datanames=None, data=None, overwrite=True)

    opt_matsim = np.squeeze(tms_opt.run(allow_multiple_runs=True, cpus=n_cpu))

    if len(opt_matsim.shape) == 2:
        opt_matsim = opt_matsim[:, :, np.newaxis]

    matlocalite = pynibs.nnav2simnibs(fn_exp_nii=fn_exp_nii,
                                      fn_conform_nii=fn_conform_nii,
                                      m_nnav=opt_matsim,
                                      target='nnav',
                                      nnav_system='localite')
    
    # write instrument marker file
    write_tms_navigator_im(opt_matsim, os.path.join(optim_dir, 'opt_coil_pos.xml'), names='opt', overwrite=True)
    # There was a mistake here: write_tms_navigator_im already flips the axes as needed to go from simnibs to localite, so
    # the call to matlocalite is superfluous and leads to a net non-change
    pynibs.create_stimsite_from_matsimnibs(os.path.join(optim_dir, 'opt_coil_pos.hdf5'), opt_matsim, overwrite=True)

else:
    print(f"Optimal coordinates already exist in:")
    print(f" > {os.path.join(optim_dir, 'opt_coil_pos.xml')}")
    print(f" > {os.path.join(optim_dir, 'opt_coil_pos.hdf5')}")
    print(f"Skipping optimization and reading optimal coordinates...")

    with h5py.File(os.path.join(optim_dir, 'opt_coil_pos.hdf5'), "r") as f:
        opt_matsim = f["matsimnibs"][:]

# save all coil positions together with electric field values in .xmdf file for plotting
fn_search_positions = os.path.join(optim_dir, "search_positions.hdf5")
fn_e = os.path.join(optim_dir, subject.id + "_TMS_optimize_" + os.path.split(fn_coil)[1].replace('.nii.gz', '_nii') +
                    ".hdf5")

with h5py.File(fn_search_positions, "r") as f:
    centers = f["centers"][:]
    m0 = f["m0"][:]
    m1 = f["m1"][:]
    m2 = f["m2"][:]

os.unlink(fn_search_positions)

# get e-fields at target for each coil position
if os.path.exists(fn_e) and not args.recompute:
    with h5py.File(fn_e, "r") as f:
        try:
            e = f["tms_optimization/E_norm"][:]  # SimNIBS 3
        except KeyError:
            e = f["tms_optimization/E_magn"][:]  # SimNIBS 4
else:
    # read e-field from .pos file
    fn_e = os.path.join(optim_dir, "coil_positions.geo")
    # line starting with 'SP'(float, float, float) {e_mag};
    exp = r"SP\([+-]?[0-9]*[.]?[0-9]+, [+-]?[0-9]*[.]?[0-9]+, [+-]?[0-9]*[.]?[0-9]+\){([+-]?[0-9]*[.]?[0-9]+)};"
    try:
        with open(fn_e, 'r') as f:
            e = []
            for line in f:
                try:
                    e.append(float(re.findall(exp, line)[0]))
                except IndexError:
                    pass
        e = np.array(e)
    except Exception as exc:
        print("Cannot read electrical fields in target for all search positions:")
        print(f"{exc}")
        e = None

pynibs.write_coil_pos_hdf5(fn_hdf=fn_search_positions,
                           centers=centers,
                           m0=m0,
                           m1=m1,
                           m2=m2,
                           datanames=["E_mag"],
                           data=e,
                           overwrite=True)

# run final simulation at optimal coil position
# ========================================================================================================

if not (os.path.exists(f'{optim_dir}{os.sep}opt_coil_pos.xml') and
        os.path.exists(f'{optim_dir}{os.sep}opt_coil_pos.hdf5') and
        os.path.exists(fn_e)) or args.recompute:
    # Setting up SimNIBS SESSION
    # ========================================================================================================
    print("Setting up SimNIBS SESSION")
    S = simnibs.sim_struct.SESSION()
    S.fnamehead = fn_mesh_msh
    S.pathfem = os.path.join(optim_dir, "e_opt")
    S.fields = "vDeE"
    S.open_in_gmsh = False
    S.fname_tensor = None
    S.map_to_surf = True
    S.map_to_fsavg = False
    S.map_to_MNI = False
    S.subpath = os.path.join(subject.mesh[mesh_id]["mesh_folder"], f"m2m_{subject.id}")

    if not os.path.exists(S.pathfem):
        os.makedirs(S.pathfem)

    # Define the TMS simulation and setting up conductivities
    # ========================================================================================================
    tms = [S.add_tmslist()]
    tms[0].fnamecoil = fn_coil
    tms[0].cond[0].value = 0.126  # WM
    tms[0].cond[1].value = 0.275  # GM
    tms[0].cond[2].value = 1.654  # CSF
    tms[0].cond[3].value = 0.01   # Skull
    tms[0].cond[4].value = 0.465  # Scalp
    tms[0].anisotropy_type = anisotropy_type

    if subject.mesh[mesh_id]['fn_tensor_vn'] is not None:
        tms[0].fn_tensor_nifti = os.path.join(subject.mesh[mesh_id]["mesh_folder"],
                                            subject.mesh[mesh_id]['fn_tensor_vn'])
    else:
        tms[0].fn_tensor_nifti = None

    tms[0].excentricity_scale = None

    # Define the optimal coil position
    # ========================================================================================================
    pos = tms[0].add_position()
    pos.matsimnibs = opt_matsim[:, :, 0]
    pos.distance = 0.1
    pos.didt = 1e6  # coil current rate of change in A/s (1e6 A/s = 1 A/us)
    pos.postprocess = S.fields

    # Running simulations
    # ========================================================================================================
    print("Running SimNIBS")
    filenames = simnibs.run_simnibs(S, cpus=1)

    print("Saving results in .hdf5 format.")
    if type(filenames) is not list:
        filenames = [filenames]

    pynibs.simnibs_results_msh2hdf5(fn_msh=filenames,
                                    fn_hdf5=[os.path.join(S.pathfem, "e")],
                                    S=S,
                                    pos_tms_idx=[0],    # indices of different "coils"
                                    pos_local_idx=[0],  # indices of associated positions
                                    subject=subject,
                                    mesh_idx=mesh_id,
                                    overwrite=True,
                                    mid2roi=True,
                                    verbose=True)
print("=== DONE ===")
