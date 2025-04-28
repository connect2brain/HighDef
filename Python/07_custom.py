import os
import pandas as pd
import numpy as np
import pynibs
import argparse
from pynibs.regression import regress_data
import h5py
from scipy.optimize import curve_fit
import multiprocessing
import time
import warnings


# # Custom functions for regressing 3D E-field
# Variants: 
# (1) Multivariate Regression: using scipy.optimize.fitcurve
# (2) Support Vector Regression: sklearn.svm.SVR

def sigmoid3D(X, slope, L, w0, w1, w2, x0):
    w = np.array([w0, w1, w2])
    return L/(1+np.exp(-slope*(w.T@X - x0)))

def abssigmoid3D(X, slope, L, w0, w1, w2, x0):
    w = np.array([w0, w1, w2])
    return L/(1+np.exp(-slope*(np.abs(w.T@X) - x0)))

def regress3D(X, y, p0=None):
    """
    Performs sigmoid regression from (3, n) predictors to 
    (n,) dependent variable
    Returns:
        qof: R² quality of fit
    """
    # Initial guess!
    # Initial guess for direction vector: forward (y=1)
    if p0 is None:
        p0 = [1, np.quantile(y, 0.8), 0, 1, 0, np.mean(np.sqrt(np.sum(X**2,axis=0)))]

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            popt, pcov = curve_fit(abssigmoid3D, X, y, p0=p0, nan_policy='omit')
            predicted = abssigmoid3D(X, *popt)
            R_squared = 1-(np.mean(np.power(predicted - y, 2)) / np.var(y))
    except RuntimeError:
        R_squared = np.nan
        popt = np.nan*np.zeros((6,))
    return R_squared, popt

class Triangle:
    def __init__(self, X, y, p0=None):
        self.X = X
        self.y = y
        self.p0 = p0

def fit(element):
    return regress3D(element.X, element.y, element.p0)

def regress_all_triangles_parallel(E, response, connectivity, n_cpu):
    """
    params:
        E: (n_trials, n_triangles, 3) matrix describing E-field at each triangle in each trial
        response: (n_trials,) labels
        connectivity: (n_triangles, 3): Matrix with three node-indices per triangle
        n_cpu: number of CPU cores to use for parallel computation
    """
    E = np.transpose(E) # reverses order of axes to: (3, n_triangles, n_trials)
    start = time.time()
    n_triangles = E.shape[1]
    triangles = [Triangle(np.squeeze(E[:,i_triangle,:]), response) for i_triangle in range(n_triangles)]


    n_cpu = min(n_cpu, multiprocessing.cpu_count(), len(response))
    pool = multiprocessing.get_context("fork").Pool(n_cpu)
    result = pool.map(fit, triangles)
    #print(result)
    fit_parameters = np.array([r[1] for r in result])
    result = np.array([r[0] for r in result])
    end = time.time()
    print(f"Done.\n\tFitting all triangles took: {end - start:2.2f} s\n\tFor {np.mean(np.isnan(result))*100:2.2f} % of the {n_triangles} Triangles, no fit was found")

    # TODO: refit nan-elements, using the average fit parameters of not-nan fits around them
    max_iter = 0
    i_iter = 0
    while i_iter < max_iter and np.any(np.isnan(result)):
        start = time.time()
        refit_mask = np.isnan(result)
        triangles_refit = []
        # For each triangle, get the average fit result of its neighbors (if there is one)
        for i_triangle in np.arange(n_triangles)[refit_mask]:
            corners = connectivity[i_triangle,:]
            neighbor_parameters = fit_parameters[np.any(np.isin(connectivity, corners), axis=1), :]
            valid_neighbors = ~np.any(np.isnan(neighbor_parameters), axis=1)
            neighbor_parameters = neighbor_parameters[valid_neighbors,:]
            if np.sum(valid_neighbors) == 0:
                p0 = None
            else:
                p0 = np.mean(neighbor_parameters, axis=0)
            triangles_refit += [Triangle(np.squeeze(E[:,i_triangle,:]), response, p0=p0)]

        refit_result = pool.map(fit, triangles_refit)
        fit_parameters[refit_mask,:] = np.array([r[1] for r in refit_result])
        result[refit_mask] = np.array([r[0] for r in refit_result])
        end = time.time()
        print(f"Done.\n\tRefitting nan triangles iteration {i_iter} took: {end - start:2.2f} s\n\tFor {np.mean(np.isnan(result))*100:2.2f} % of the {n_triangles} Triangles, no fit was found")
        i_iter += 1
    # TODO: add refit discontinuities etc.
    return result












# Adapted from 07_calc_r2

parser = argparse.ArgumentParser(description='Creates a new subject in path')
parser.add_argument('-f', '--fn_subject', help='path of new subject', required=True, type=str)
parser.add_argument('-e', '--exp_id', help='id of the experiment', required=True, type=str)
parser.add_argument('-m', '--mesh_id', help='Mesh ID', required=True, type=str)
parser.add_argument('-r', '--roi_id', help='ROI ID', required=True, type=str)
parser.add_argument('-n', '--n_cpu', help='How many cpus to use', type=int, default=21)
parser.add_argument('-s', '--splits', help='Whether to do a 1st half/2nd half split in addition to the default whole-dataset regression (assess intra-subject stability)', action="store_true")
args = parser.parse_args()

fn_subject = os.path.abspath(args.fn_subject)
subject_id = os.path.split(fn_subject)[1]
exp_id = args.exp_id
mesh_id = args.mesh_id
roi_id = args.roi_id
n_cpu = args.n_cpu
do_splits = args.splits

#fn_subject = 'TMS_localization/TMS_loc_results/sub-001/sub-001.hdf5'
#exp_id = 'main'
#mesh_id = 'mesh0'
#roi_id = "midlayer_larger"
#fn_raw_data_table = "sub-001_raw.csv"
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

subject = pynibs.load_subject(fn_subject)
subject_dir = subject.subject_folder
fn_raw_data_table = subject.exp[exp_id]['response_csv']

raw_data = pd.read_csv(os.path.join(subject_dir, fn_raw_data_table))

n_trials = raw_data.shape[0]
trial_indices = {"all": np.arange(n_trials)}
if do_splits:
    print("Doing splits!")
    midpoint = n_trials // 2
    trial_indices["first_half"]  = np.arange(midpoint)
    trial_indices["second_half"] = np.arange(midpoint, n_trials)
else:
    print("Not doing splits!")






functions = [getattr(pynibs, 'sigmoid4')]

# quantities of electric field the fit is performed for
e_qoi_list = ["E_mag"] # "E", 

# map data to whole brain surface
map2surf = False # TODO: check if this brain surface is individual or average! -- with present code, it's the individual brain --- and its unclear how to use this code to get to inflated fs average (although this is mentioned -- but not explained-- in the pynibs documentation)

# TODO: there is some interesting cortical-layer-specific stuff!
# https://direct.mit.edu/imag/article/doi/10.1162/imag_a_00036/118104/Directional-sensitivity-of-cortical-neurons
layerid = None
neuronmodel = None
waveform = 'biphasic'

n_refit = 20
score_type = "R2"

convergence = False # is this needed? -- i think i should compare this between excitability and inhibition.
verbose = True


## Setting up files:
# parent folder containing the output data
fn_results = os.path.join(subject_dir, "results", f"exp_{exp_id}", "r2",
                          f"mesh_{mesh_id}", f"roi_{roi_id}")
fn_result_suffix = "_neuron" if layerid is not None else ""

# file containing the electric field in the ROI [n_trials (?) x n_components (?)]
fn_e_roi = os.path.join(subject_dir, "results", f"exp_{exp_id}", "electric_field",
                       f"mesh_{mesh_id}", f"roi_{roi_id}", f"e_neuron_scaled.hdf5")

if not os.path.exists(fn_e_roi):
    fn_e_roi = os.path.join(subject_dir, "results", f"exp_{exp_id}", "electric_field",
                            f"mesh_{mesh_id}", f"roi_{roi_id}", f"e.hdf5")
    if not os.path.exists(fn_e_roi):
        raise FileNotFoundError(f"No electric field file found!")

fn_result_suffix = f"_neuron_{layerid}_{neuronmodel}" if layerid is not None else ""

subject.exp[exp_id]["fn_exp_hdf5"] = [os.path.join(subject.subject_folder, "exp", exp_id, mesh_id, "experiment.hdf5")]



print(f"Processing subject: {subject.id}, experiment: {exp_id}")
print("=====================================================")

# load ROI
roi = pynibs.load_roi_surface_obj_from_hdf5(subject.mesh[mesh_id]['fn_mesh_hdf5'])[roi_id]
con = roi.node_number_list
print(f' Roi number list = {con}')
mesh_folder = subject.mesh[mesh_id]["mesh_folder"]

convergence_results = {'metric': [],
                       'qoi': [],
                       'fun': [],
                       'nrmsd': [],
                       'geodesic dist': []}



for i_m, m in enumerate(subject.exp[exp_id]['response_columns']):
    for split_name, split_indices in trial_indices.items():
        print(f"\n> Split: {split_name}")
        print(f"> Measure: {m}")

        for i_fun, fun in enumerate(functions):
            if neuronmodel and not (fun == pynibs.sigmoid4 or fun == pynibs.sigmoid4_log):
                raise NotImplementedError("Neuron regression can only be "
                                        "performed with sigmoid4 or sigmoid4_log functions.")

            print(f"> fun: {fun.__name__}")

            gof_dict = dict()

            for i_e_qoi, e_qoi in enumerate(e_qoi_list):
                print(f"> E-QOI: {e_qoi}")

                if e_qoi == "E_norm":
                    select_signed_data = True
                else:
                    select_signed_data = False

                # load electric field
                with h5py.File(fn_e_roi, "r") as f:
                    # Here: Neuron-model/Layer_id stuff was removed, bc. not needed for now.
                    e_matrix = f[e_qoi][:]
                    print(f'Shape of e_matrix for {e_qoi}: {e_matrix.shape}') # (n_trials, n_triangles, 3) for E
                    nodes = roi.node_coord_mid
                    con = roi.node_number_list


                if m.startswith("SIHIscore") or m.startswith("CsE"):
                    muscle = m.split("_")[1]
                else:
                    muscle = m.split("_")[2]

                assert(muscle in ["APB", "FDI", "ADM"], f"Unknown muscle: {muscle}")
                PrI_column_name = f"preinnervation_{muscle}_in_uV"
                PrI_threshold = 50 # uV
                CR_column_name = f'inhibited_mep_{muscle}_in_uV'
                CR_threshold = 40 # uV

                # load MEPs
                mep = np.array(raw_data[m])
                preinnervation = np.array(raw_data[PrI_column_name])

                # Restrict to current split
                mep = mep[split_indices]
                e_matrix = e_matrix[split_indices, :]
                preinnervation = preinnervation[split_indices]
                PrI_too_high = preinnervation > PrI_threshold

                # Reject trials with too little conditioned response!
                if m.startswith('SIHIscore'):
                    conditioned_response = np.array(raw_data[CR_column_name])
                    conditioned_response = conditioned_response[split_indices] # Restrict to current split

                    CR_too_low   = conditioned_response < CR_threshold
                    rejected = np.logical_or(CR_too_low, PrI_too_high)
                    print(f'{split_name}: Rejected {100*np.mean(rejected)} % of trials\n - {100*np.mean(CR_too_low)} % for lack of conditioned response (< {CR_threshold} uV)  andor\n - {100*np.mean(PrI_too_high)} % for pre-innervation (> {PrI_threshold} uV)')
                    mep = mep[np.logical_not(rejected)]
                    e_matrix = e_matrix[np.logical_not(rejected), :]
                else:
                    # Reject trials with too much pre-innervation
                    print(f'{split_name}: Rejected {100*np.mean(PrI_too_high)} % of trials for pre-innervation (> {PrI_threshold} uV)')
                    mep = mep[np.logical_not(PrI_too_high)]
                    e_matrix = e_matrix[np.logical_not(PrI_too_high), :]

                if m.startswith("inhibited"):
                    mep = np.max(mep)-mep

                if len(mep) == 0:
                    print(f"{split_name}: NO TRIALS REMAIN AFTER TRIAL-REJECTION!")
                    continue


                # check for zero e-fields and filter them (problems in FEM!)
                zero_mask = (e_matrix == 0).all(axis=1)

                if zero_mask.any():
                    print(f"\nWarning! {np.sum(zero_mask)} zero e-fields detected in element! Check FEM! Ignoring them for now!\n\n")
                    e_matrix = e_matrix[np.logical_not(zero_mask), :]
                    mep = mep[np.logical_not(zero_mask)]

                # import matplotlib
                # import matplotlib.pyplot as plt
                # matplotlib.use('Qt5Agg')
                # hotspot_idx = 12045
                # # hotspot_idx = 8148
                # # hotspot_idx = 8426
                # # hotspot_idx = 12004
                # # hotspot_idx = 18599
                # plt.scatter(e_matrix_orig[:, hotspot_idx], mep)
                # plt.scatter(e_matrix[:, hotspot_idx], mep)
                    



                if layerid is not None and neuronmodel == "IOcurve":
                    gof_dict[e_qoi] = 1/np.var(e_matrix, axis=0)

                elif e_qoi == 'E':
                    gof_dict[e_qoi] = regress_all_triangles_parallel(e_matrix, mep, con, n_cpu)
                else:
                    gof, fit_result = regress_data(e_matrix=e_matrix,
                                                mep=mep,
                                                fun=fun,
                                                n_cpu=n_cpu,
                                                con=con,
                                                n_refit=n_refit,
                                                return_fits=True,
                                                score_type=score_type,
                                                verbose=verbose,
                                                pool=None,
                                                refit_discontinuities=True,
                                                select_signed_data=select_signed_data)
                    gof_dict[e_qoi] = gof
                    if verbose:
                        i_best = np.argmax(gof)
                        print(f'\n- - - - - - - - -\nBest Triangle with {fit_result[i_best]}\n')
                        

                if convergence:
                    print(f"Computing n-1 results to assess convergence.")
                    if e_qoi == 'E':
                        res_prev = regress_all_triangles_parallel(e_matrix[:-1,:,:], mep[:-1], con, n_cpu)
                    else:
                        res_prev = regress_data(e_matrix=e_matrix[:-1, :],
                                                mep=mep[:-1],
                                                fun=fun,
                                                n_cpu=n_cpu,
                                                con=con,
                                                n_refit=n_refit,
                                                return_fits=False,
                                                score_type=score_type,
                                                verbose=verbose,
                                                pool=None,
                                                refit_discontinuities=True,
                                                select_signed_data=select_signed_data)

                    # compute normalised root-mean-square deviation between n and n-1 solution
                    norm_to_perc = 99
                    res_final = gof_dict[e_qoi] / np.percentile(gof_dict[e_qoi][gof_dict[e_qoi] > 0],
                                                                norm_to_perc)
                    res_prev = res_prev / np.percentile(res_prev[res_prev > 0], norm_to_perc)
                    nrmsd = pynibs.nrmsd(res_final, res_prev)

                    # compute geodesic distance between best element n and n-1 solution
                    best_elm = np.argmax(res_final)
                    best_elm_prev = np.argmax(res_prev)
                    nodes_dist, tris_dist = pynibs.geodesic_dist(nodes,
                                                                con,
                                                                best_elm,
                                                                source_is_node=False)
                    geod_dist = tris_dist[best_elm_prev]

                    convergence_results['metric'].append(m)
                    convergence_results['qoi'].append(e_qoi)
                    convergence_results['fun'].append(fun.__name__)
                    convergence_results['nrmsd'].append(nrmsd)
                    convergence_results['geodesic dist'].append(geod_dist)

            # save results
            if split_name == "all":
                # In this case store in the standard location
                fn_gof_hdf5 = os.path.join(fn_results, m, fun.__name__, f"r2{fn_result_suffix}_roi_data.hdf5")
            else:
                # Only for the new splits: Create new sub-folders
                fn_gof_hdf5 = os.path.join(fn_results, m, fun.__name__, split_name, f"r2{fn_result_suffix}_roi_data.hdf5")

            # create folder
            if not os.path.exists(os.path.split(fn_gof_hdf5)[0]):
                os.makedirs(os.path.split(fn_gof_hdf5)[0])

            # write hdf5 _geo file
            pynibs.write_geo_hdf5_surf(out_fn=os.path.join(os.path.split(fn_gof_hdf5)[0], f"r2{fn_result_suffix}_roi_geo.hdf5"),
                                    points=nodes,
                                    con=con,
                                    replace=True,
                                    hdf5_path='/mesh')

            # write _data file
            data = [gof_dict[e_qoi] for e_qoi in e_qoi_list]
            data_names = [f'c_{e_qoi}' for e_qoi in e_qoi_list]

            print(f"Writing results to {fn_gof_hdf5}.")
            pynibs.write_data_hdf5_surf(data=data,
                                        data_names=data_names,
                                        data_hdf_fn_out=fn_gof_hdf5,
                                        geo_hdf_fn=os.path.join(os.path.split(fn_gof_hdf5)[0], f"r2{fn_result_suffix}_roi_geo.hdf5"),
                                        replace=True)

            if not map2surf or layerid is not None: # skip this step also when mapping on the layers is requested
                continue
            print("> Mapping data to brain surface")
            c_mapped = pynibs.map_data_to_surface(datasets=data,
                                                points_datasets=[roi.node_coord_mid] * len(data),
                                                con_datasets=[roi.node_number_list] * len(data),
                                                fname_fsl_gm=[None, None],
                                                fname_fsl_wm=[None, None],
                                                fname_midlayer=[
                                                    os.path.join(mesh_folder, subject.mesh[mesh_id]['fn_lh_midlayer']),
                                                    os.path.join(mesh_folder, subject.mesh[mesh_id]['fn_rh_midlayer'])
                                                ],
                                                delta=subject.roi[mesh_id][roi_id]['delta'],
                                                input_data_in_center=True,
                                                return_data_in_center=True,
                                                data_substitute=-1)

            # recreate complete midlayer surface to write in .hdf5 geo file
            points_midlayer, con_midlayer = pynibs.make_GM_WM_surface(gm_surf_fname=[subject.mesh[mesh_id]['fn_lh_gm'],
                                                                                    subject.mesh[mesh_id]['fn_rh_gm']],
                                                                    wm_surf_fname=[subject.mesh[mesh_id]['fn_lh_wm'],
                                                                                    subject.mesh[mesh_id]['fn_rh_wm']],
                                                                    mesh_folder=mesh_folder,
                                                                    midlayer_surf_fname=[
                                                                        subject.mesh[mesh_id]['fn_lh_midlayer'],
                                                                        subject.mesh[mesh_id]['fn_rh_midlayer']
                                                                    ],
                                                                    delta=subject.roi[mesh_id][roi_id]['delta'],
                                                                    x_roi=None,
                                                                    y_roi=None,
                                                                    z_roi=None,
                                                                    layer=1,
                                                                    fn_mask=None)

            # save hdf5 _geo file (mapped)
            print(f" > Writing mapped brain and roi to {os.path.join(os.path.split(fn_gof_hdf5)[0], f'r2{fn_result_suffix}_data.hdf5')}")
            pynibs.write_geo_hdf5_surf(out_fn=os.path.join(os.path.split(fn_gof_hdf5)[0], f"r2{fn_result_suffix}_geo.hdf5"),
                                    points=points_midlayer,
                                    con=con_midlayer,
                                    replace=True,
                                    hdf5_path='/mesh')

            # save hdf5 _data file (mapped)

            pynibs.write_data_hdf5_surf(data=c_mapped,
                                        data_names=data_names,
                                        data_hdf_fn_out=os.path.join(os.path.split(fn_gof_hdf5)[0], f"r2{fn_result_suffix}_data.hdf5"),
                                        geo_hdf_fn=os.path.join(os.path.split(fn_gof_hdf5)[0], f"r2{fn_result_suffix}_geo.hdf5"),
                                        replace=True)

print("DONE")

if convergence:
    print("Convergence for n vs n-1 stimulations:")
    convergence_results = pd.DataFrame().from_dict(convergence_results)
    print(convergence_results)
    convergence_fn = os.path.join(fn_results, f"convergence{fn_result_suffix}.csv")
    print(f"Saved to {convergence_fn}")
    convergence_results.to_csv(convergence_fn)
print("====")
