import h5py
import numpy as np
import pandas as pd
from time import sleep
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import argparse

parser = argparse.ArgumentParser(description='Evaluates all experiments registered for given subject (also reruns create_subject)')
parser.add_argument('-s', '--subject', help='subject ID', default='sub-001', type=str)
parser.add_argument('-e', '--exp_id', help='id of the experiment', default='map-L2R', type=str)
args = parser.parse_args()

sub_id = args.subject
exp_id = args.exp_id
roi_id = f"midlayer_{exp_id.split('-')[-1][0].lower()}"
mesh_id = "mesh0"
cse_name = "CsE_FDI_in_uV"
sihi_name = "SIHIscore_FDI"

print(f'Data: {sub_id}\n      {exp_id}\n      {mesh_id}\n      {roi_id}\n  y = {cse_name} ; {sihi_name}\n')
# Structure inside that file:
# '/E'      : 3 x nTriangles x nTrials
# '/E_mag'  :     nTriangles x nTrials
# '/E_norm' :     nTriangles x nTrials
# '/E_tan'  :     nTriangles x nTrials
e_field_file = f"/home/bnplab-admin/TMS_localization/HighDef/{sub_id}/results/exp_{exp_id}/electric_field/mesh_{mesh_id}/roi_{roi_id}/e.hdf5"
with h5py.File(e_field_file, 'r') as e_field_data:
    predictors = e_field_data['/E_mag'][:]
    # this is: (nTrials, nTriangles)

# Structure:
# '/data/tris/c_E_mag' : nTriangles
r_squared_file = f"/home/bnplab-admin/TMS_localization/HighDef/{sub_id}/results/exp_{exp_id}/r2/mesh_{mesh_id}/roi_{roi_id}/{cse_name}/sigmoid4/r2_roi_data.hdf5"
with h5py.File(r_squared_file, 'r') as r_squared_data:
    qof = r_squared_data['/data/tris/c_E_mag'][:]
    # this is: (nTriangles,)
best_R2_cse = np.max(qof)
print(f'Weise best fit {cse_name}:\t\t R² = {best_R2_cse})')
best_triangle_cse = np.argmax(qof)

# Structure:
# '/data/tris/c_E_mag' : nTriangles
r_squared_file = f"/home/bnplab-admin/TMS_localization/HighDef/{sub_id}/results/exp_{exp_id}/r2/mesh_{mesh_id}/roi_{roi_id}/{sihi_name}/sigmoid4/r2_roi_data.hdf5"
with h5py.File(r_squared_file, 'r') as r_squared_data:
    qof = r_squared_data['/data/tris/c_E_mag'][:]
    # this is: (nTriangles,)
best_R2_sihi = np.max(qof)
print(f'Weise best fit {cse_name}:\t\t R² = {best_R2_sihi})')
best_triangle_sihi = np.argmax(qof)

response_file = f'/home/bnplab-admin/TMS_localization/HighDef/{sub_id}/{sub_id}_{exp_id.split("-")[-1]}_raw.csv'
response_data = pd.read_csv(response_file)



x = predictors[:, best_triangle_cse]
cse = response_data[cse_name]
sihi = response_data[sihi_name]

# [A.1] Sigmoid fit
def sigmoid(x, L ,x0, k, b):
    y = L / (1 + np.exp(-k*(x-x0))) + b
    return (y)
p0 = [max(cse), np.median(x), (max(cse) - min(cse)) / (max(x) - min(x)), min(cse)] # this is an mandatory initial guess
popt, pcov = curve_fit(sigmoid, x, cse, p0, method='lm')
R2_cse_raw = 1 - (np.mean((cse - sigmoid(x,*popt))**2) / np.var(cse))
print(f'Curve fit {cse_name}:\t\t R² = {R2_cse_raw}')

# Formula from: https://statproofbook.github.io/D/snr.html
SNR_cse_raw = np.var(sigmoid(x,*popt)) / np.var(cse - sigmoid(x,*popt))
print(f'Curve fit {cse_name}:\t\tSNR = {SNR_cse_raw}')

p0 = [max(sihi), np.median(x), (max(sihi) - min(sihi)) / (max(x) - min(x)), min(sihi)] # this is an mandatory initial guess
popt, pcov = curve_fit(sigmoid, x, sihi, p0, method='lm')
R2_sihi_raw = 1 - (np.mean((sihi - sigmoid(x,*popt))**2) / np.var(sihi))
print(f'Curve fit {sihi_name}:\t\t R² = {R2_sihi_raw}')

# Formula from: https://statproofbook.github.io/D/snr.html
SNR_sihi_raw = np.var(sigmoid(x,*popt)) / np.var(sihi - sigmoid(x,*popt))
print(f'Curve fit {sihi_name}:\t\tSNR = {SNR_sihi_raw}')





rng = np.random.default_rng()
n_query_points = 100
R2  = np.zeros((n_query_points,)) * np.nan
SNR = np.zeros((n_query_points,)) * np.nan

for i, sigma in enumerate(np.linspace(0, max(cse)-min(cse), n_query_points)):
    print(f"it {i} / {n_query_points}   (scale={sigma})")
    y_ = cse + rng.normal(loc=0, scale=sigma, size=cse.shape)
    p0 = [max(y_), np.median(x), (max(y_) - min(y_)) / (max(x) - min(x)), min(y_)] # this is an mandatory initial guess
    try:
        popt, pcov = curve_fit(sigmoid, x, y_, p0, method='lm')
        y_pred = sigmoid(x,*popt)
        R2[i]  = 1 - (np.mean((y_ - y_pred)**2) / np.var(y_))
        SNR[i] = np.var(y_pred) / np.var(y_ - y_pred)
    except RuntimeError as e:
        print(f"FIT FAILED in  {i} / {n_query_points}   (scale={sigma})")
        pass
    

#xq = np.linspace(min(SNR), max(SNR), 100)
#plt.plot(xq, 1-1/xq, color="#f0f0f0")
# ^ the line is in fact not that!
plt.plot(SNR, R2, 'k.', markersize=1)
plt.plot(SNR_cse_raw, R2_cse_raw, 'r+')
plt.plot(SNR_sihi_raw, R2_sihi_raw, 'b+')
plt.xlabel("SNR")
plt.ylabel("R²")
plt.show()
