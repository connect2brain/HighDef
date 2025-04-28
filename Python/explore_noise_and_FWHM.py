import h5py
import numpy as np
import pandas as pd
from time import sleep
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import argparse
from multiprocessing import Pool

def sigmoid(x, L ,x0, k, b):
    y = L / (1 + np.exp(-k*(x-x0))) + b
    return (y)


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

response_file = f'/home/bnplab-admin/TMS_localization/HighDef/{sub_id}/{sub_id}_{exp_id.split("-")[-1]}_raw.csv'
response_data = pd.read_csv(response_file)


n_trials, n_triangles = predictors.shape

cse = response_data[cse_name]
sihi = response_data[sihi_name]

rng = np.random.default_rng()
n_query_points = 30

R2  = np.zeros((n_query_points,n_triangles)) * np.nan
SNR = np.zeros((n_query_points,n_triangles)) * np.nan

for i, sigma in enumerate(np.linspace(0, 3*(max(cse)-min(cse)), n_query_points)):
    print(f'...')
    added_noise = rng.normal(loc=0, scale=sigma, size=cse.shape)
    for j in range(n_triangles):
        x = predictors[:,j]
        y_ = cse + added_noise
        p0 = [max(y_), np.median(x), (max(y_) - min(y_)) / (max(x) - min(x)), min(y_)] # this is an mandatory initial guess
        try:
            popt, pcov = curve_fit(sigmoid, x, y_, p0, method='lm')
            y_pred = sigmoid(x,*popt)
            R2[i,j]  = 1 - (np.mean((y_ - y_pred)**2) / np.var(y_))
            SNR[i,j] = np.var(y_pred) / np.var(y_ - y_pred)
        except RuntimeError as e:
            pass
    print(f'\nFinished {i+1} / {n_query_points}')


np.savetxt('.backup_R2', R2)
np.savetxt('.backup_SNR', SNR)







# Structure:
# '/data/tris/c_E_mag' : nTriangles
r_squared_file = f"/home/bnplab-admin/TMS_localization/HighDef/{sub_id}/results/exp_{exp_id}/r2/mesh_{mesh_id}/roi_{roi_id}/{cse_name}/sigmoid4/r2_roi_data.hdf5"
with h5py.File(r_squared_file, 'r') as r_squared_data:
    qof_cse = r_squared_data['/data/tris/c_E_mag'][:]
    # this is: (nTriangles,)
best_triangle_cse = np.argmax(qof_cse)

x = predictors[:, best_triangle_cse]
p0 = [max(cse), np.median(x), (max(cse) - min(cse)) / (max(x) - min(x)), min(cse)] # this is an mandatory initial guess
popt, pcov = curve_fit(sigmoid, x, cse, p0, method='lm')
R2_cse_raw = 1 - (np.mean((cse - sigmoid(x,*popt))**2) / np.var(cse))
SNR_cse_raw = np.var(sigmoid(x,*popt)) / np.var(cse - sigmoid(x,*popt))
N_cse = np.sum(qof_cse > 0.5*R2_cse_raw)

# Structure:
# '/data/tris/c_E_mag' : nTriangles
r_squared_file = f"/home/bnplab-admin/TMS_localization/HighDef/{sub_id}/results/exp_{exp_id}/r2/mesh_{mesh_id}/roi_{roi_id}/{sihi_name}/sigmoid4/r2_roi_data.hdf5"
with h5py.File(r_squared_file, 'r') as r_squared_data:
    qof_sihi = r_squared_data['/data/tris/c_E_mag'][:]
    # this is: (nTriangles,)
best_triangle_sihi = np.argmax(qof_sihi)

x = predictors[:, best_triangle_sihi]
p0 = [max(sihi), np.median(x), (max(sihi) - min(sihi)) / (max(x) - min(x)), min(sihi)] # this is an mandatory initial guess
popt, pcov = curve_fit(sigmoid, x, sihi, p0, method='lm')
R2_sihi_raw = 1 - (np.mean((sihi - sigmoid(x,*popt))**2) / np.var(sihi))
SNR_sihi_raw = np.var(sigmoid(x,*popt)) / np.var(sihi - sigmoid(x,*popt))
N_sihi = np.sum(qof_sihi > 0.5*R2_sihi_raw)


SNR_of_best_triangle = np.diag(SNR[:,np.nanargmax(R2, axis=1)])
n_triangles_with_any_R2 = np.sum(~np.isnan(R2), axis=1)
n_triangles_with_good_R2 = np.sum(R2 > np.tile(0.5 * np.nanmax(R2, axis=1), (n_triangles, 1)).T, axis=1)


plt.plot(SNR_of_best_triangle, n_triangles_with_good_R2 / n_triangles_with_any_R2, 'k.', markersize=1)
plt.plot(SNR_cse_raw,  N_cse  / np.sum(~np.isnan(qof_cse)),  'r+')
plt.plot(SNR_sihi_raw, N_sihi / np.sum(~np.isnan(qof_sihi)), 'b+')

plt.xlabel("SNR")
plt.ylabel("FWHM")
plt.show()
