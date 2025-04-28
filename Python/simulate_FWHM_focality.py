# In this script:
# (1) Get the best fitting triangle from the SIHI-map
#     Compute the fit (E -> Response) on this triangle with the real data
#     Further compute some model of the noise (check if heteroskedastic, if so, maybe some histogram based thing)
#   + Compute the FWHM of the raw data
# (2) From the E-field in each trial in this best triangle, compute the "clean response" (directly from fit) and add error (following distribution above) as new "response"
# (3) Run the Weise regression script on this simulated dataset
#     Compute FWHM of result and collect this!
# (4) Plot the distribution of the simulated FWHM values and the FWHM of the raw data



import h5py
import numpy as np
import pandas as pd
from time import sleep
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.io import savemat
import pynibs
from pynibs.regression import regress_data
import argparse
import multiprocessing


parser = argparse.ArgumentParser(description='Creates a new subject in path')
parser.add_argument('-s', '--subject', help='id of subject', required=True, type=str)
parser.add_argument('-e', '--exp_id', help='id of the experiment', required=True, type=str)
parser.add_argument('-m', '--mesh_id', help='Mesh ID', required=False, type=str, default="mesh0")
parser.add_argument('-r', '--roi_id', help='ROI ID', required=False, type=str, default=None)
parser.add_argument('-o', '--response', help='Response (either CsE_FDI_in_uV or SIHIscore_FDI)', required=False, type=str, default='SIHIscore_FDI')
args = parser.parse_args()


exp_id = args.exp_id
mesh_id = args.mesh_id
roi_id = args.roi_id

n_refit = 20
score_type = "R2"


sub_id = args.subject
exp_id = args.exp_id
if args.roi_id is None:
    roi_id = f"midlayer_{exp_id.split('-')[-1][0].lower()}"
else:
    roi_id = args.roi_id
mesh_id = args.mesh_id
response = args.response #f"SIHIscore_{muscle}" # f"CsE_{muscle}_in_uV" # 
muscle = response.split("_")[1]

PrI_column_name = f"preinnervation_{muscle}_in_uV"
PrI_threshold = 50 # uV
CR_column_name = f'inhibited_mep_{muscle}_in_uV'
CR_threshold = 40 # uV

if response.startswith("SIHI"):
    # For SIHI_score
    relative_kernel_width = 0.5
    relative_noise_width  = 0.05
else:
    # For CsE:
    relative_kernel_width = 0.05
    relative_noise_width  = 0.05


print(f'Data: {sub_id}\n      {exp_id}\n      {mesh_id}\n      {roi_id}\n  y = {response}\n')

# Structure inside that file:
# '/E'      : 3 x nTriangles x nTrials
# '/E_mag'  :     nTriangles x nTrials
# '/E_norm' :     nTriangles x nTrials
# '/E_tan'  :     nTriangles x nTrials
ROOT = "/mnt/d/HighDef-operate/HighDef"
e_field_file = f"{ROOT}/{sub_id}/results/exp_{exp_id}/electric_field/mesh_{mesh_id}/roi_{roi_id}/e.hdf5"
print(f'Opening: {e_field_file}')
with h5py.File(e_field_file, 'r') as e_field_data:
    predictors = e_field_data['/E_mag'][:]
    # this is: (nTrials, nTriangles)

# Structure:
# '/data/tris/c_E_mag' : nTriangles
r_squared_file = f"{ROOT}/{sub_id}/results/exp_{exp_id}/r2/mesh_{mesh_id}/roi_{roi_id}/{response}/sigmoid4/r2_roi_data.hdf5"
with h5py.File(r_squared_file, 'r') as r_squared_data:
    qof = r_squared_data['/data/tris/c_E_mag'][:]
    # this is: (nTriangles,)

response_file = f'{ROOT}/{sub_id}/{sub_id}_{exp_id.split("-")[-1]}_raw.csv'
response_data = pd.read_csv(response_file)

# "response_data" and "predictors" are used below, so rejecting here is sufficient!


best_R2 = np.max(qof)
print(f'Weise best fit:\t\tR² = {best_R2}')
best_triangle = np.argmax(qof)


preinnervation = np.array(response_data[PrI_column_name])
PrI_too_high = preinnervation > PrI_threshold

if response.startswith('SIHIscore'):
    CR_too_low = response_data[CR_column_name] < CR_threshold # uV
    rejected = np.logical_or(CR_too_low, PrI_too_high)
else:
    rejected = PrI_too_high

accepted = np.logical_not(rejected)

x = np.array(predictors[accepted, best_triangle])
y = np.array(response_data[response])[accepted]



FWHM_raw = np.mean(qof > 0.5*np.max(qof))
print(f'FWHM on actual data:\tFWHM = {FWHM_raw}')





# [1] Fit the model

def sigmoid(x, L ,x0, k, b):
    y = L / (1 + np.exp(-k*(x-x0))) + b
    return (y)

p0 = [max(y), np.median(x), (max(y) - min(y)) / (max(x) - min(x)), min(y)] # this is an mandatory initial guess
popt, pcov = curve_fit(sigmoid, x, y, p0, method='lm')

residuals = y - sigmoid(x,*popt)

print(f'Std of residuals: \to² = {np.std(residuals)}')
R2_sigmoid_direct = 1 - (np.mean(residuals**2) / np.var(y))
print(f'Curve fit:\t\tR² = {R2_sigmoid_direct}')

"""
plt.subplot(2,2,1)
q = np.linspace(min(x), max(x), 100)
plt.plot(x, y, 'k+')
plt.plot(q, sigmoid(q, *popt), 'r-')
plt.xlabel("x")
plt.ylabel("y")

plt.subplot(2,2,2)
plt.plot(x, residuals, 'k+')
plt.xlabel("x")
plt.ylabel("Residual")
"""



def sample_heteroskedastic_noise_at(xq, x, residuals, relative_kernel_width=0.2, relative_noise_width=0.5, rng=None):
    """
    xq: At which value (single value!) to sample from the heteroskedastic noise specified by x and residuals
    x: x-value at which the respective residual was observed
    residuals: residual at the respective x-value
    relative_kernel_width: width (std) of the kernel along the x-axis, that specifies the weights of the residuals
    relative_noise_width: width (std) of the normal noise added to the sampled residual, expressed as a ratio of the range of the residuals within +- kernel_width
    """
    if rng is None:
        rng = np.random.default_rng()

    range_x = max(x) - min(x)
    kernel_width = range_x * relative_kernel_width
    p = np.exp(-(0.5/kernel_width**2) * (x - xq)**2)
    p = p / np.sum(p)

    mask = np.logical_and(x > xq - kernel_width, x < xq + kernel_width)
    noise_width = relative_noise_width * (max(residuals[mask]) - min(residuals[mask]))

    return rng.choice(residuals, size=(1,), p=p) + rng.normal(loc=0, scale=noise_width, size=(1,))



rng = np.random.default_rng()
"""
# Variant (b): Heteroskedastic noise (risk of overfitting)
plt.subplot(2,2,3)
yq = [sigmoid(xq, *popt) + sample_heteroskedastic_noise_at(xq, x, residuals, relative_noise_width=relative_noise_width, relative_kernel_width=relative_kernel_width, rng=rng) for xq in x]
plt.plot(x, yq, 'g+')
yq = [sigmoid(xq, *popt) + sample_heteroskedastic_noise_at(xq, x, residuals, relative_noise_width=relative_noise_width, relative_kernel_width=relative_kernel_width, rng=rng) for xq in x]
plt.plot(x, yq, 'c+')
yq = [sigmoid(xq, *popt) + sample_heteroskedastic_noise_at(xq, x, residuals, relative_noise_width=relative_noise_width, relative_kernel_width=relative_kernel_width, rng=rng) for xq in x]
plt.plot(x, yq, 'y+')
plt.plot(q, sigmoid(q, *popt), 'r-')
plt.xlabel("x")
plt.ylabel("y")


# Variant (a): Homoskedastic noise (simplification)
plt.subplot(2,2,4)
plt.plot(x, sigmoid(x, *popt) + rng.normal(loc=0, scale=np.std(residuals), size=x.shape), 'k+')
plt.plot(q, sigmoid(q, *popt), 'r-')
plt.xlabel("x")
plt.ylabel("y")


plt.show()
"""

#noise_mode = "homoskedastic"
noise_mode = "heteroskedastic"


subject = pynibs.load_subject(f'{ROOT}/{sub_id}/{sub_id}.hdf5')
roi     = pynibs.load_roi_surface_obj_from_hdf5(subject.mesh[mesh_id]['fn_mesh_hdf5'])[roi_id]
nodes   = roi.node_coord_mid
con     = roi.node_number_list

#plt.ion()
#fig, ax = plt.subplots()
#ax.plot(1, FWHM_raw, 'r+')
#plt.pause(0.1)
#plt.draw()

MC_iterations = round(1e2)

gofs = np.zeros((MC_iterations,))
fwhms = np.zeros((MC_iterations,))

n_cpu = 30

n_cpu_available = multiprocessing.cpu_count()
n_cpu = min(n_cpu, n_cpu_available)
pool = multiprocessing.get_context("fork").Pool(n_cpu)


for it in range(MC_iterations):
    PrI_too_high = np.array(response_data[PrI_column_name]) > PrI_threshold
    if response.startswith("SIHI"):
        CR_too_low = response_data[CR_column_name] < CR_threshold # uV
        rejected = np.logical_or(CR_too_low, PrI_too_high)
    else:
        rejected = PrI_too_high
    accepted = np.logical_not(rejected)

    # Simulate the responses under the assumption that they only depend on the E-field at the best-matching triangle
    E_mag_at_best_triangle = np.array(predictors[accepted, best_triangle])

    if noise_mode == "homoskedastic":
        noise = rng.normal(loc=0, scale=np.std(residuals), size=E_mag_at_best_triangle.shape)
    else:
        noise = np.reshape(np.array([sample_heteroskedastic_noise_at(xq, E_mag_at_best_triangle, residuals, relative_noise_width=relative_noise_width, relative_kernel_width=relative_kernel_width, rng=rng) for xq in E_mag_at_best_triangle]), E_mag_at_best_triangle.shape)

    response_from_focal_source = sigmoid(E_mag_at_best_triangle, *popt) + noise
    # This is the e_field_matrix
    e_matrix = predictors[accepted, :]
    mep = response_from_focal_source # E_mag_at_best_triangle is already rejected


    # check for zero e-fields and filter them (problems in FEM!)
    zero_mask = (e_matrix == 0).all(axis=1)

    if zero_mask.any():
        print(f"\nWarning! {np.sum(zero_mask)} zero e-fields detected in element! Check FEM! Ignoring them for now!\n\n")
        e_matrix = e_matrix[np.logical_not(zero_mask), :]
        mep = mep[np.logical_not(zero_mask)]
    
    gof, fit_result = regress_data(e_matrix=e_matrix,
                                            mep=mep,
                                            fun=getattr(pynibs, 'sigmoid4'),
                                            n_cpu=n_cpu,
                                            con=con,
                                            n_refit=n_refit,
                                            return_fits=True,
                                            score_type=score_type,
                                            verbose=True,
                                            pool=pool,
                                            refit_discontinuities=True,
                                            select_signed_data=False)
    gofs[it] = np.max(gof)
    fwhms[it] = np.mean(gof > 0.5 * np.max(gof))
    
    i_best = np.argmax(gof)
    print(f'\n- - - - - - - - -\n [{it}/{MC_iterations}] Best Triangle with {fit_result[i_best]}\n FWHM = {fwhms[it]} (vs. true {FWHM_raw})\n')

    #ax.plot(0, fwhms[it], 'k.')
    #plt.pause(0.1)
    #plt.draw()


#plt.close()
savemat(f'{ROOT}/{sub_id}/{sub_id}_{exp_id}_{response}_simulated_FWHMs_{noise_mode}.mat', {'simulated_fwhms': fwhms, 'FWHM_raw': FWHM_raw})
#np.savetxt(f'/home/bnplab-admin/TMS_localization/HighDef/{sub_id}/{sub_id}_{exp_id}_{response}_simulated_FWHMs_{noise_mode}.txt', fwhms)


counts, bins = np.histogram(fwhms)

plt.stairs(counts / np.sum(counts), bins)
plt.vlines(FWHM_raw, 0, 1, 'r')
plt.ylim(0,1.1*np.max(counts/np.sum(counts)))
plt.xlabel('FWHM')
plt.ylabel('Probability under assumption of full focality')

plt.savefig(f'/mnt/d/HighDef-operate/Figures/{sub_id}_{exp_id}_{response}_simulated_FWHM_{noise_mode}.png')










