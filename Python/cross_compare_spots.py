# In this script:
# (1) Select the best-fit triangles from the R²-result files (+ADM, +APB, +FDI, -FDI)
# (2) Get optimal coil position of ADM MEP, APB MEP, FDI MEP, FDI SIHI (call these coil positions 8@+ADM, 8@+APB, 8@+FDI, 8@-FDI)  // Figure-of-eight at ...
# (3) Compute E-field at +ADM, +APB, +FDI, -FDI when coil is placed at each of those 4 locations (4x4)
# (4) Collect nearest M samples (closest in |E|) as the samples for rmANOVA
# Do this for all subjects (4,4,n_subjects)

# In the end, i want: a "p-value" for each comparison (+ADM vs +APB, +ADM vs +FDI, +APB vs +FDI, +FDI vs -FDI)
# to this end,

import h5py
import numpy as np
from glob import glob
import pandas as pd
from functools import partial, reduce
import operator
from simnibs import mesh_io
from simnibs.optimization import opt_struct
from simnibs.simulation import fem, sim_struct, save_matlab_sim_struct
import pynibs
from pynibs.util.simnibs import calc_e_in_midlayer_roi
import os
import datetime

import matplotlib.pyplot as plt

fn_config = '/home/bnplab-admin/TMS_localization/config.json'
config = pd.read_json(fn_config)

n_samples = 50

def encode(response_name: str) -> str:
    if response_name.startswith("CsE"):
        return "+" + response_name[4:7]
    else:
        return "-" + response_name[-3:]

subjects  = ["sub-001", "sub-002", "sub-003", "sub-004", "sub-006", "sub-007", "sub-008", "sub-009", "sub-011", "sub-012", "sub-013", "sub-014", "sub-015", "sub-016", "sub-017", "sub-019", "sub-022", "sub-023"]
responses = ["CsE_ADM_in_uV", "CsE_APB_in_uV", "CsE_FDI_in_uV", "SIHIscore_FDI"]
contrasts = reduce(operator.add, [[[A, B] for B in responses if A != B and ((A.startswith("CsE") and B.startswith("CsE")) or ("FDI" in A and "FDI" in B))] for A in responses], [])
experiments = ["map-L2R", "map-R2L"]
print("Contrasting: ")
print("\n".join([f'{p}' for p in contrasts]))

root = '/mnt/d/HighDef-operate/HighDef'

selected_trials    = np.zeros((len(contrasts), 2, n_samples, len(experiments), len(subjects)), dtype=np.int64)
selected_responses = np.zeros((len(contrasts), 2, n_samples, len(experiments), len(subjects))) * np.nan


for i_subject, ID in enumerate(subjects):
    fn_subject = f'{root}/{ID}/{ID}.hdf5'
    subject = pynibs.load_subject(fn_subject)

    print(f'{"_".join(30*["_"])}\n ID: {ID}')

    for i_exp, exp_id in enumerate(experiments):
        roi_id = f'midlayer_{exp_id[-3].lower()}'
        # (0) Read geometry
        # Needed for checking that the correct opt-file exists and is used!
        with h5py.File(f'{root}/{ID}/mesh/roi/{roi_id}/geo.hdf5') as f:
            vertex_coordinates = f['/mesh/nodes/node_coord'][:]
            triangles = f['/mesh/elm/triangle_number_list'][:]

        response_data = pd.read_csv(os.path.join(subject.subject_folder, subject.exp[exp_id]['response_csv']))

        # (1)+(2)
        best_triangle_indices = np.zeros((len(responses),), dtype=np.int64)
        matsimnibs = np.zeros((4, 4, len(responses)))
        for i_response, response in enumerate(responses):
            print(f'{ID} {exp_id} {response}')
            # (1) get the spots:
            with h5py.File(f'{root}/{ID}/results/exp_{exp_id}/r2/mesh_mesh0/roi_{roi_id}/{response}/sigmoid4/r2_roi_data.hdf5') as f:
                qof = f[f'data/tris/c_E_mag'][:]
                if np.max(qof) < config.R2_rejection.only_apply_if_best_R2_is_lower_than:
                    lqof = np.log(f[f'data/tris/c_E_mag'][:])
                    Q3 = np.quantile(lqof, 0.75)
                    iqr = Q3 - np.quantile(lqof, 0.25)

                    banned_elements = lqof > (config.R2_rejection.Cutoff_IQR_factor*iqr + Q3)
                    
                    qof[banned_elements] = 0

                best_triangle_indices[i_response] = np.argmax(qof)
            
            best_triangle_center = np.mean(vertex_coordinates[triangles[best_triangle_indices[i_response],:],:], axis=0)
            
            
            # (2) read in optimal coil locations
            optimized_file = f'{root}/{ID}/opt/{exp_id}/{response}_[{" ".join([f"{c:+.2f}" for c in best_triangle_center])}]/mesh0/E_mag/opt_coil_pos.hdf5'
            print(f'Trying to open {optimized_file}', end='')
            with h5py.File(optimized_file) as f:
                matsimnibs[:, :, i_response] = np.squeeze(f['/matsimnibs'][:])
            print(' --- Found')

        # Now have: 
        # best_triangle_indices  i.e.   +ADM,   +APB,   +FDI,   -FDI
        # matsimnibs             i.e. 8@+ADM, 8@+APB, 8@+FDI, 8@-FDI

        with h5py.File(f'{root}/{ID}/results/exp_{exp_id}/electric_field/mesh_mesh0/roi_{roi_id}/ereshaped.hdf5') as f:
            simulated_E_mag = f['/E_mag'][:]
            # ^ Shape: (n_trials, n_triangles)
            simulated_E_mag = simulated_E_mag[:, best_triangle_indices]
            # ^ Shape: (n_trials, 4) -- one column for each spot.



        mesh_simn = mesh_io.read_msh(subject.mesh['mesh0']['fn_mesh_msh'])

        n_cpu = 28
        anisotropy_type = "scalar"

        tms_opt = opt_struct.TMSoptimize()
        tms_opt.fnamehead = subject.mesh['mesh0']["fn_mesh_msh"]
        tms_opt.pathfem   = '/mnt/d/HighDef-operate/Figures'
        tms_opt.fnamecoil = '/mnt/c/Users/bnplab-admin/SimNIBS-4.0/simnibs_env/Lib/site-packages/simnibs/resources/coil_models/Drakaki_BrainStim_2022/MagVenture_Cool-B35.ccd'
        tms_opt.target = None
        tms_opt.open_in_gmsh = False
        tms_opt.anisotropy_type = anisotropy_type

        print(f"Preparing conductivity tensors")
        conductivities = sim_struct.SimuList.cond2elmdata(tms_opt, mesh=mesh_simn)

        roi = pynibs.load_roi_surface_obj_from_hdf5(subject.mesh['mesh0']['fn_mesh_hdf5'])[roi_id]
        layer_gm_wm_info = None
        postpro = partial(
                calc_e_in_midlayer_roi,
                roi=roi,
                mesh=mesh_simn,
                qoi=['E', 'mag', 'norm', 'tan'],
                layer_gm_wm_info=layer_gm_wm_info,
        )

        pos_matrices = [matsimnibs[:, :, i] for i in range(matsimnibs.shape[-1])]        
        
        # (A) Run sim once with fixed, stereotypical intensity
        I0 = 50
        didt_list = matsimnibs.shape[-1] * [I0]
        dataset = '/tmp'
        fn_hdf5 = f'/mnt/d/HighDef-operate/Figures/E_{ID}_exp_{exp_id}.hdf5'
        if os.path.exists(fn_hdf5):
            os.rename(fn_hdf5, f'{os.path.splitext(fn_hdf5)[0]}_outdated_at_{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}.hdf5')

        # Next: Compute E-field when coil is placed at these 8@s. 
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
            field='E')
        
        with h5py.File(fn_hdf5, 'r') as h5:
            # undo zero padding by storing just as many values as elements on layer
            num_roi_elmts = roi.node_number_list.shape[0]
            E_mag = np.squeeze(h5[dataset][:, :num_roi_elmts, 3])

        slopes = E_mag[:, best_triangle_indices] / I0
        # slopes: (from 8@X, to Y)

        

        # In a contrast (contrastant_1, contrastant_2), the coil is placed on 8@contrastant_1, and on 8@contrastant_2
        # and the response obtained from contrastant_1 is compared between these two locations
        for i_contrast, (referent, other) in enumerate(contrasts):
            print(f'Contrasting {referent} and {other}')
            i_referent = responses.index(referent)
            i_other = responses.index(other)

            # pre-innervation rejection w.r.t. contrastant_1
            target_muscle = referent.split("_")[1]
            print(f"w.r.t. muscle {target_muscle}")
            PrI_column_name = f"preinnervation_{target_muscle}_in_uV"
            PrI_threshold = 50 # uV
            preinnervation = np.array(response_data[PrI_column_name])
            PrI_too_high = preinnervation > PrI_threshold
            accepted = np.logical_not(PrI_too_high)


            # E-fields at best triangle for response 1
            simulated_E_at_ref = simulated_E_mag[accepted, i_referent]

            Intensities = np.linspace(0, 100, 100)
            difference_of_mean_response = np.zeros((len(Intensities),))
            for i_I, I in enumerate(Intensities):
                E_from_ref_at_ref = slopes[i_referent, i_referent] * I
                E_from_other_at_ref = slopes[i_other, i_referent] * I

                # Select those trials which have the closest simulated E field at the reference spot (triangle)
                closest_trials_ref   = np.argsort(np.abs(simulated_E_at_ref - E_from_ref_at_ref))[:n_samples]
                closest_trials_other = np.argsort(np.abs(simulated_E_at_ref - E_from_other_at_ref))[:n_samples]

                difference_of_mean_response[i_I] = np.mean((np.array(response_data[referent])[accepted])[closest_trials_ref]) - np.mean((np.array(response_data[referent])[accepted])[closest_trials_other])

            most_favorable_intensity = Intensities[np.argmax(difference_of_mean_response)]
            print(f' > Highest difference for {most_favorable_intensity} % MSO')
            if most_favorable_intensity < 10.0:
                print(f'   [!] Rejected due to being too low (implausible), setting to 50 % MSO')
                most_favorable_intensity = 50.0
            
            E_from_ref_at_ref   = slopes[i_referent, i_referent] * most_favorable_intensity
            E_from_other_at_ref = slopes[i_other, i_referent]    * most_favorable_intensity

            trial_labels = np.where(accepted)[0]
            selected_trials[i_contrast, 0, :, i_exp, i_subject] = trial_labels[np.argsort(np.abs(simulated_E_at_ref - E_from_ref_at_ref))[:n_samples]]
            selected_trials[i_contrast, 1, :, i_exp, i_subject] = trial_labels[np.argsort(np.abs(simulated_E_at_ref - E_from_other_at_ref))[:n_samples]]
            for i in [0, 1]:
                selected_responses[i_contrast, i, :, i_exp, i_subject] = np.array(response_data[referent])[selected_trials[i_contrast, i, :, i_exp, i_subject]]

        
# Format selected data into long form table:
# Coil_spot, Receiver_spot (this implies the response type!), Response, Subject, Coil_hemisphere, Trial
# // this is CsE:
# +FDI, +FDI, 1000, sub-001, L, 4
# +FDI, +FDI, 1050, sub-001, L, 42
# ...
# -FDI, +FDI,  900, sub-001, L, 100
# -FDI, +FDI,  800, sub-001, L, 12
# ...
# // this is SIHI:
# +FDI, -FDI,  1.5, sub-001, L, 5
# -FDI, -FDI,  2.5, sub-001, L, 50

data = {'Subject': [], 'Hemisphere': [], 'Trial': [], 'Coil': [], 'Receiver': [], 'Response': [], 'Contrast': []}


for i_subject, s in enumerate(subjects):
    for i_exp, e in enumerate(experiments):
        hemi = e.split('-')[-1][0]
        for i_contrast, contrast in enumerate(contrasts):
            ref = contrast[0]
            for i, c in enumerate(contrast):
                data['Coil']       += n_samples * [encode(c)]
                data['Receiver']   += n_samples * [encode(ref)]
                data['Subject']    += n_samples * [s]
                data['Hemisphere'] += n_samples * [hemi]
                data['Trial']      += selected_trials[i_contrast, i, :, i_exp, i_subject].tolist()
                data['Response']   += selected_responses[i_contrast, i, :, i_exp, i_subject].tolist()
                data['Contrast']   += n_samples * ['_vs_'.join([encode(v) for v in contrast])]

df = pd.DataFrame(data)
df.to_csv(f'/mnt/d/HighDef-operate/Figures/contrasts.csv')

print(df)


