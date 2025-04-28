import argparse, os
from datetime import datetime
import numpy as np
import pandas as pd

import h5py
import pynibs
from pynibs.util.simnibs import calc_e_in_midlayer_roi
from pynibs.regression import regress_data
import pygpc
from pygpc.AbstractModel import AbstractModel

from collections import OrderedDict
from functools import partial

from simnibs import opt_struct, mesh_io
from simnibs.optimization import opt_struct
from simnibs.simulation import fem, sim_struct, save_matlab_sim_struct


import logging


class PipelineModel(AbstractModel):
    # Parameters:
    # "c_wm": Conductivity of white matter
    # "c_gm": Conductivity of grey matter
    # "c_csf": Conductivity of cerebrospinal fluid
    # "dx_<i>": Deviation of real x-coordinate of the coil from the recorded position in trial i; in mm
    # "dy_<i>": Deviation of real y-coordinate of the coil from the recorded position in trial i; in mm
    # "dz_<i>": Deviation of real z-coordinate of the coil from the recorded position in trial i; in mm
    def __init__(self, subject, mesh, roi, raw_data, coil_matrices, config):
        super().__init__()
        self.subject  = subject
        self.raw_data = raw_data
        self.mesh = mesh
        self.roi  = roi
        self.coil_matrices = coil_matrices # a list (!) of 2D-arrays
        self.config = config
        self.how_often_evaluated = 0

    def validate(self):
        pass

    def simulate(self, process_id=None, matlab_engine=None):
        # (1) Run-many E-field simulations, with the given conductivities and coil locations
        # (2) Do regressions:
        #     Reject trials by preinnervation in the respective target muscle!
        #  (2.a) Regress for hotspot  --> R² map --> max location: hx, hy, hz
        #  (2.b) Regress for coldspot --> R² map --> max location: cx, cy, cz
        # (3) Compute distance [in posterolateral direction] --- that's the return

        # Iterate over the parameter-arrays:
        n_simulations = len(self.p["c_csf"])

        for key, value in self.p.items():
            self.p[key] = np.array(value).flatten()


        # Do all the simulations       

        distances_hotcold = np.zeros((n_simulations,)) * np.nan
        for i_simulation in range(n_simulations):

            ### (1): E - F I E L D   S I M U L A T I O N  
            
            t = datetime.now()
            logger.info(f"Starting E-field simulation {i_simulation} at {t.strftime('%Y-%m-%d %H:%M:%S')}")
            fn_out = os.path.join(self.subject.subject_folder, "UQ", f"run-{self.config['timestamp']}", f"exp_{self.config['exp_id']}", "electric_field", f"mesh_{self.config['mesh_id']}", f"roi_{self.config['roi_id']}", f"iteration_{i_simulation}")
            fn_regression_results = os.path.join(self.subject.subject_folder, "UQ", f"run-{self.config['timestamp']}", f"exp_{self.config['exp_id']}", "r2", f"mesh_{self.config['mesh_id']}", f"roi_{self.config['roi_id']}", f"iteration_{i_simulation}")
            if not os.path.exists(fn_out):
                os.makedirs(fn_out)
            fn_e = os.path.join(fn_out, f"e.hdf5")
            fn_hdf5 = fn_e #os.path.join(fn_out, fn_e)
            logger.info(f"fn_hdf5 = {fn_hdf5}")
            if os.path.exists(fn_hdf5):
                os.rename(fn_hdf5, f'{os.path.splitext(fn_hdf5)[0]}.rn{datetime.now().strftime("%Y%m%d_%H%M%S")}.hdf5')  

            tms_opt = opt_struct.TMSoptimize()
            tms_opt.fnamehead = subject.mesh[mesh_id]["fn_mesh_msh"]
            tms_opt.pathfem   = fn_out
            tms_opt.fnamecoil = '/mnt/c/Users/bnplab-admin/SimNIBS-4.0/simnibs_env/Lib/site-packages/simnibs/resources/coil_models/Drakaki_BrainStim_2022/MagVenture_Cool-B35.ccd'
            tms_opt.target = None
            tms_opt.open_in_gmsh = False
            tms_opt.anisotropy_type = "scalar"

            # [0] Scalp
            # [1] Skull
            tms_opt.cond[2].value = self.p["c_csf"][i_simulation]  # CSF
            tms_opt.cond[3].value = self.p["c_gm"][i_simulation]   # Gray matter
            tms_opt.cond[4].value = self.p["c_wm"][i_simulation]   # White matter


            
            matsimnibs_list = [M + np.array([[0,0,0,self.p[f"dx_{i}"]][i_simulation],
                                            [0,0,0,self.p[f"dy_{i}"]][i_simulation],
                                            [0,0,0,self.p[f"dz_{i}"]][i_simulation],
                                            [0, 0, 0, 0]]) 
                                            for i, M in enumerate(self.coil_matrices)]
            

            postpro = partial(
                    calc_e_in_midlayer_roi,
                    roi=self.roi,
                    mesh=self.mesh,
                    qoi=["E", "mag"],
                    layer_gm_wm_info=None,
            )

            didt_list = [0.96*i for i in self.raw_data['Intensity_percentMSO']]
            conductivities = sim_struct.SimuList.cond2elmdata(tms_opt, mesh=self.mesh)

            dataset = "/tmp"
            

            fem.tms_many_simulations(
                mesh=self.mesh,
                cond=conductivities,
                fn_coil=tms_opt.fnamecoil,
                matsimnibs_list=matsimnibs_list,
                didt_list=didt_list,
                fn_hdf5=fn_hdf5,
                dataset=dataset,
                post_pro=postpro,
                solver_options=tms_opt.solver_options,
                n_workers=self.config["n_cpu"],
                field="E")
            logger.info(f"Done. E-field simulation took {datetime.now() - t}")

            with h5py.File(fn_hdf5, 'r') as h5:
                # undo zero padding by storing just as many values as elements on layer
                num_roi_elmts = roi.node_number_list.shape[0]
                E = h5[dataset][:, :num_roi_elmts, 0:3]
                E_mag = h5[dataset][:, :num_roi_elmts, 3]

                qoi_idx = 6

            with h5py.File(fn_hdf5.replace('.hdf5','reshaped.hdf5'), 'w') as h5:
                h5.create_dataset('/E', data=E)
                h5.create_dataset('/E_mag', data=E_mag)

            with h5py.File(fn_e, "a") as f:
                logger.info(f"Reshaping fields.")
                E = f['tmp'][:, :, 0:3]
                E_mag = f['tmp'][:, :, 3]
                f.create_dataset('/E', data=E)
                f.create_dataset('/E_mag', data=E_mag)



            ### (2) R E G R E S S I O N S
            fun = getattr(pynibs, 'sigmoid4')
            e_qoi = "E_mag"
            n_refit = 20
            score_type = "R2"
            verbose = False
            select_signed_data = False

            ### (2.a) CsE_FDI_in_uV
            m = "CsE_FDI_in_uV"

            # load electric field
            with h5py.File(fn_e, "r") as f:
                # Here: Neuron-model/Layer_id stuff was removed, bc. not needed for now.
                e_matrix = f[e_qoi][:]
                logger.info(f'Shape of e_matrix for {e_qoi}: {e_matrix.shape}') # (n_trials, n_triangles, 3) for E
                nodes = roi.node_coord_mid
                con = roi.node_number_list
            
            muscle = m.split("_")[1]
            PrI_column_name = f"preinnervation_{muscle}_in_uV"
            PrI_threshold = 50 # uV

            # load MEPs
            mep = np.array(raw_data[m])
            preinnervation = np.array(raw_data[PrI_column_name])

            PrI_too_high = preinnervation > PrI_threshold
            mep = mep[np.logical_not(PrI_too_high)]
            e_matrix = e_matrix[np.logical_not(PrI_too_high), :]

            # check for zero e-fields and filter them (problems in FEM!)
            zero_mask = (e_matrix == 0).all(axis=1)

            if zero_mask.any():
                logger.warning(f"\nWarning! {np.sum(zero_mask)} zero e-fields detected in element! Check FEM! Ignoring them for now!\n\n")
                e_matrix = e_matrix[np.logical_not(zero_mask), :]
                mep = mep[np.logical_not(zero_mask)]

            gof, fit_result = regress_data(e_matrix=e_matrix,
                                                mep=mep,
                                                fun=fun,
                                                n_cpu=self.config["n_cpu"],
                                                con=con,
                                                n_refit=n_refit,
                                                return_fits=True,
                                                score_type=score_type,
                                                verbose=verbose,
                                                pool=None,
                                                refit_discontinuities=True,
                                                select_signed_data=select_signed_data)
            i_best = np.argmax(gof)
            nodes_of_best_element = nodes[i_best]
            location_hotspot = np.mean(con[nodes_of_best_element,:], axis=0)
            logger.info(f"iteration {i_simulation}: Hotspot found at {location_hotspot} with fit: {fit_result}")

            ### (2.a) SIHIscore_FDI
            m = "SIHIscore_FDI"

            # load electric field
            with h5py.File(fn_e, "r") as f:
                # Here: Neuron-model/Layer_id stuff was removed, bc. not needed for now.
                e_matrix = f[e_qoi][:]
                logger.info(f'Shape of e_matrix for {e_qoi}: {e_matrix.shape}') # (n_trials, n_triangles, 3) for E
                nodes = roi.node_coord_mid
                con = roi.node_number_list
            
            muscle = m.split("_")[1]
            PrI_column_name = f"preinnervation_{muscle}_in_uV"
            PrI_threshold   = 50 # uV
            CR_column_name  = f'inhibited_mep_{muscle}_in_uV'
            CR_threshold    = 40 # uV

            # load MEPs
            mep = np.array(raw_data[m])
            preinnervation = np.array(raw_data[PrI_column_name])
            conditioned_response = np.array(raw_data[CR_column_name])

            PrI_too_high = preinnervation > PrI_threshold
            CR_too_low   = conditioned_response < CR_threshold
            rejected     = np.logical_or(CR_too_low, PrI_too_high)
            mep          = mep[np.logical_not(rejected)]
            e_matrix     = e_matrix[np.logical_not(rejected), :]

            # check for zero e-fields and filter them (problems in FEM!)
            zero_mask = (e_matrix == 0).all(axis=1)

            if zero_mask.any():
                logger.warning(f"\nWarning! {np.sum(zero_mask)} zero e-fields detected in element! Check FEM! Ignoring them for now!\n\n")
                e_matrix = e_matrix[np.logical_not(zero_mask), :]
                mep = mep[np.logical_not(zero_mask)]

            gof, fit_result = regress_data(e_matrix=e_matrix,
                                                mep=mep,
                                                fun=fun,
                                                n_cpu=self.config["n_cpu"],
                                                con=con,
                                                n_refit=n_refit,
                                                return_fits=True,
                                                score_type=score_type,
                                                verbose=verbose,
                                                pool=None,
                                                refit_discontinuities=True,
                                                select_signed_data=select_signed_data)
            i_best = np.argmax(gof)
            nodes_of_best_element = nodes[i_best]
            location_coldspot = np.mean(con[nodes_of_best_element,:], axis=0)
            logger.info(f"iteration {i_simulation}: Coldspot found at {location_coldspot} with fit: {fit_result}")

            ### (3) C O M P U T E   D I S T A N C E
            if hemisphere == "l":
                compare_along = np.array([-1, -1, -1]) / np.sqrt(3)
            else:
                compare_along = np.array([ 1, -1, -1]) / np.sqrt(3)
            connector = location_coldspot - location_hotspot
            
            # project onto posterolateral and get signed distance in posterolateral direction:
            distances_hotcold[i_simulation] = connector @ compare_along
            logger.info(f"Distance in {i_simulation}  =  {distances_hotcold[i_simulation]} mm  in posterolateral direction")
            logger.info(f"Iteration {i_simulation} took {datetime.now() - t}")
            self.how_often_evaluated += 1

        return distances_hotcold

        





















parser = argparse.ArgumentParser(description='Runs generalized polynomial chaos expansion for Uncertainty quantification')
parser.add_argument('-f', '--fn_subject', help='path of subject', required=True, type=str)
parser.add_argument('-e', '--exp_id', help='id of the experiment', required=True, type=str)
parser.add_argument('-n', '--n_cpu', help='How many cpus to use', type=int, default=21)
args = parser.parse_args()



logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
ROOT = "/home/bnplab-admin/TMS_localization/HighDef"
file_handler = logging.FileHandler(f"{ROOT}/gpc_run-{datetime.now().strftime('%Y-%m-%dT%H%M%S')}.log", encoding='utf-8')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(name)s: %(message)s')
file_handler.setFormatter(formatter)
file_handler.setLevel(logging.DEBUG)
logger.addHandler(file_handler)



fn_subject = os.path.abspath(args.fn_subject)
subject_id = os.path.split(fn_subject)[1]
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






model = PipelineModel(subject, mesh_simn, roi, raw_data, pos_matrices, config)

parameters = OrderedDict()
parameters["c_wm"]  = pygpc.Beta(pdf_shape=[3, 3], pdf_limits=[0.1,0.4])
parameters["c_gm"]  = pygpc.Beta(pdf_shape=[3, 3], pdf_limits=[0.1,0.6]) 
parameters["c_csf"] = pygpc.Beta(pdf_shape=[3, 3], pdf_limits=[1.2,1.8])

for i in range(n_trials):
    parameters[f"dx_{i}"] = pygpc.Norm(pdf_shape=[0, 0.350]) # in mm
    parameters[f"dy_{i}"] = pygpc.Norm(pdf_shape=[0, 0.350])
    parameters[f"dz_{i}"] = pygpc.Norm(pdf_shape=[0, 0.350])

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
options["fn_results"] = "tmp/mygpc" 
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



### Turns out: This is too heavy for the PC. The model cannot even be initialized
# CMD out:
# Initializing gPC object...
# Killed

# This is because of the countless (3*~900 =~ 2700) dx,dy,dz parameters for the uncertainty of merely the location of the coil (not even orientation at this point)


