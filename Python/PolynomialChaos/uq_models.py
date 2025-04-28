
import pynibs
from pynibs.regression import regress_data
from pynibs.util.simnibs import calc_e_in_midlayer_roi

from simnibs.simulation import gpc, gPC_regression

from pygpc.AbstractModel import AbstractModel

import numpy as np
from datetime import datetime
import logging

import os, glob

logger = logging.getLogger(__name__)




class SubModel:
    def __init__(self, regression: gPC_regression, roi, cropped):
        self.regression = regression
        self.coefficients = self.regression.expand_quantity(build_get_midlayer_E_mag(roi, cropped))

    def evaluate(self, xi):
        return self.regression.evaluate(self.coefficients, xi)
    

def build_get_midlayer_E_mag(roi, cropped):
    def get_midlayer_E_mag(Es):
        E_mag = np.zeros((Es.shape[0], roi.n_tris))
        for i, E in enumerate(Es):
            res = calc_e_in_midlayer_roi(E, roi=roi, mesh=cropped, qoi=["mag"], layer_gm_wm_info=None)
            assert(len(res.shape) < 2 or res.shape[1] < 2)
            E_mag[i,:] = res.flatten()
        return E_mag
    return get_midlayer_E_mag



class WrapperModel(AbstractModel):
    # Parameters:
    # "c_wm":  Conductivity of white matter           --  follows beta distribution over [-1, 1]
    # "c_gm":  Conductivity of grey matter            --  follows beta distribution over [-1, 1]
    # "c_csf": Conductivity of cerebrospinal fluid    --  follows beta distribution over [-1, 1]
    # "c_cb":  Conductivity of compact bone           --  follows beta distribution over [-1, 1]
    # Note: the actual range of the parameters is defined when creating the wrapped models! I.e. basically "compile" the sub-models with desired ranges, and here only employ them.
    def __init__(self, subject, mesh, roi, raw_data, config):
        super().__init__()
        self.subject  = subject
        self.raw_data = raw_data
        self.mesh = mesh
        self.roi  = roi
        self.config = config
        self.how_often_evaluated = 0

        self.n_trials = raw_data.shape[0]
        self.sub_models = [] # a list of SubModel instances!
        self.load_submodels()

    def load_submodels(self):
        self.sub_models = []
        offset = 0
        current_endpoint = self.n_trials
        for i in range(self.n_trials):
            logger.info(f"gathering gPC for trial {i+1}")
            
            if i == current_endpoint:
                offset = current_endpoint
                current_endpoint = self.n_trials

            partial_path = f'/mnt/d/HighDef-operate/TMS_UQ/raw/{self.subject.id}/exp_{self.config["exp_id"]}'
            if offset == 0:
                startname = "first"
            else:
                startname = str(offset)
            fs = glob.glob(f"{partial_path}/trials_{startname}-*")
            assert(len(fs) == 1, f'Expected exactly one sub-folder with name trials_{startname}-*, but found {len(fs)}: {fs}')
            endpoint = os.path.split(fs[0])[-1].split("-")[-1]
            if endpoint != "last":
                current_endpoint = int(endpoint)
            else:
                current_endpoint = self.n_trials

            fn_hdf5 = f"{fs[0]}/e_{(i+1-offset):04.0f}_gpc.hdf5"
            print(f"gathering gPC for trial {i+1} from {fn_hdf5}")
            # Read the regression object from the HDF5 file
            self.sub_models += [SubModel(gpc.gPC_regression.read_hdf5(fn_hdf5), self.roi, self.mesh)]




            



            



    def validate(self):
        pass

    def simulate(self, process_id=None, matlab_engine=None):
        # (1) Evaluate the wrapped models at the given conductivity values
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
        distances_hotcold = np.zeros((n_simulations,1)) * np.nan
        for i_simulation in range(n_simulations):
            t = datetime.now()
            #### (1) Evaluate the submodels
            conductivities = np.array([self.p["c_wm"], self.p["c_gm"], self.p["c_csf"], self.p["c_cb"]]).reshape((1, -1)) # -> (1, n_params)
            E_mag = np.zeros((self.n_trials, self.roi.n_tris)) * np.nan
            # Go through all per-trial gPC-models, and sample E_mag in ROI (n_roi_elements,) from each
            # for regression, need (n_trials, n_roi_elements)
            for i_trial in range(self.n_trials):
                E_mag[i_trial,:] = self.sub_models[i_trial].evaluate(conductivities)

            ### (2) R E G R E S S I O N S
            location_hotspot, fit_result = self.best_spot("CsE_FDI_in_uV", E_mag.copy())
            logger.info(f"iteration {i_simulation}: Hotspot found at {location_hotspot} with fit: {fit_result}")

            location_coldspot, fit_result = self.best_spot("SIHIscore_FDI", E_mag.copy())
            logger.info(f"iteration {i_simulation}: Coldspot found at {location_coldspot} with fit: {fit_result}")

            ### (3) C O M P U T E   D I S T A N C E
            if self.config["hemisphere"] == "l":
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


    def best_spot(self, response, e_matrix):
        fun = getattr(pynibs, 'sigmoid4')
        n_refit = 20
        score_type = "R2"
        verbose = False
        select_signed_data = False
    
        nodes = self.roi.node_coord_mid
        con = self.roi.node_number_list
        
        muscle = response.split("_")[1]
        PrI_column_name = f"preinnervation_{muscle}_in_uV"
        PrI_threshold = 50 # uV
        CR_column_name = f'inhibited_mep_{muscle}_in_uV'
        CR_threshold = 40 # uV

        # load MEPs
        mep = np.array(self.raw_data[response])
        preinnervation = np.array(self.raw_data[PrI_column_name])

        PrI_too_high = preinnervation > PrI_threshold

        # Reject trials with too little conditioned response!
        if response.startswith('SIHIscore'):
            conditioned_response = np.array(self.raw_data[CR_column_name])

            CR_too_low   = conditioned_response < CR_threshold
            rejected = np.logical_or(CR_too_low, PrI_too_high)
            mep = mep[np.logical_not(rejected)]
            e_matrix = e_matrix[np.logical_not(rejected), :]
        else:
            mep = mep[np.logical_not(PrI_too_high)]
            e_matrix = e_matrix[np.logical_not(PrI_too_high), :]

        if response.startswith("inhibited"):
            mep = np.max(mep)-mep



        # check for zero e-fields and filter them (problems in FEM!)
        zero_mask = (e_matrix == 0).all(axis=1)

        if zero_mask.any():
            logger.warning(f"\nWarning! {np.sum(zero_mask)} zero e-fields detected in element! Check FEM! Ignoring them for now!\n\n")
            e_matrix = e_matrix[np.logical_not(zero_mask), :]
            mep = mep[np.logical_not(zero_mask)]

        print(f"e_matrix: {e_matrix.shape}, mep: {mep.shape}")

        gof, fit_result = regress_data( e_matrix=e_matrix,
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
        nodes_of_best_element = con[i_best,:]
        location = np.mean(nodes[nodes_of_best_element,:], axis=0)
        return location, fit_result