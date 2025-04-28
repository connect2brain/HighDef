# use conda env embeddings

import h5py
import numpy as np
import pandas as pd
from time import sleep
import matplotlib.pyplot as plt
import pacmap
from sklearn.decomposition import PCA
from scipy.optimize import curve_fit
from sklearn.ensemble import RandomForestRegressor, HistGradientBoostingRegressor
from sklearn.model_selection import KFold
from sklearn.inspection import permutation_importance
from scipy.io import savemat


sub_id = "sub-003"
exp_id = "map-L2R"
roi_id = f"midlayer_{exp_id.split('-')[-1][0].lower()}"
mesh_id = "mesh0"
response = "CsE_FDI_in_uV"

print(f'Data: {sub_id}\n      {exp_id}\n      {mesh_id}\n      {roi_id}\n  y = {response}\n')

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
r_squared_file = f"/home/bnplab-admin/TMS_localization/HighDef/{sub_id}/results/exp_{exp_id}/r2/mesh_{mesh_id}/roi_{roi_id}/{response}/sigmoid4/r2_roi_data.hdf5"
with h5py.File(r_squared_file, 'r') as r_squared_data:
    qof = r_squared_data['/data/tris/c_E_mag'][:]
    # this is: (nTriangles,)

response_file = f'/home/bnplab-admin/TMS_localization/HighDef/{sub_id}/{sub_id}_{exp_id.split("-")[-1]}_raw.csv'
response_data = pd.read_csv(response_file)








best_R2 = np.max(qof)
print(f'Weise best fit:\t\tR² = {best_R2})')
best_triangle = np.argmax(qof)
x = predictors[:, best_triangle]
y = response_data[response]


# --- #
# [A] #   Do a first fit to get an R²
# --- #

# [A.1] Sigmoid fit
def sigmoid(x, L ,x0, k, b):
    y = L / (1 + np.exp(-k*(x-x0))) + b
    return (y)

p0 = [max(y), np.median(x), (max(y) - min(y)) / (max(x) - min(x)), min(y)] # this is an mandatory initial guess
popt, pcov = curve_fit(sigmoid, x, y, p0, method='lm')
R2_sigmoid_direct = 1 - (np.mean((y - sigmoid(x,*popt))**2) / np.var(y))
print(f'Curve fit:\t\tR² = {R2_sigmoid_direct}')

# [A.1.b] Trying out k-fold to check R²
n_folds = 10
kf = KFold(n_splits=n_folds, shuffle=True)
R2 = np.zeros((n_folds,))
for j, (train_index, test_index) in enumerate(kf.split(x)):
    x_train = x[train_index]
    y_train = y[train_index]
    x_test  = x[test_index]
    y_test  = y[test_index]
    p0 = [max(y_train), np.median(x_train), (max(y_train) - min(y_train)) / (max(x_train) - min(x_train)), min(y_train)] # this is an mandatory initial guess
    popt, pcov = curve_fit(sigmoid, x_train, y_train, p0, method='lm')
    R2[j] = 1 - (np.mean((y_test - sigmoid(x_test,*popt))**2) / np.var(y_test))
R2_kfold_sigmoid = np.mean(R2)
print(f'Curve fit:\t\t{n_folds}-fold CV: R² = {R2_kfold_sigmoid}')


max_depth = 3 # allow trees only up to depth
min_percent_of_total_samples_for_leaf = 0.02
max_bins = 16
n_jobs = 26
n_trees = 100

# [A.2] Random Forest
random_forest = RandomForestRegressor(max_depth=max_depth, min_samples_leaf=min_percent_of_total_samples_for_leaf, n_jobs=n_jobs, n_estimators=n_trees)

x_ = np.reshape(x, (len(x),1))
random_forest.fit(x_,y)
R2 = 1 - (np.mean((y - random_forest.predict(x_))**2) / np.var(y))
print(f'Random Forest fit:\tR² = {R2}')

# [A.3] Histogram Gradient boosted Random Forest

hgbt_regressor = HistGradientBoostingRegressor(max_depth=max_depth, min_samples_leaf=round(min_percent_of_total_samples_for_leaf * len(x_)), max_bins=max_bins)
hgbt_regressor.fit(x_, y)
R2 = 1 - (np.mean((y - hgbt_regressor.predict(x_))**2) / np.var(y))
print(f'HGBT fit:\t\tR² = {R2}')



fig = plt.figure()
#plt.ion()

ax = fig.add_subplot(1,1,1)
ax.scatter(x,y, s=2)
q = np.linspace(np.min(x), np.max(x), 200)
q_ = np.reshape(q, (len(q), 1))
ax.plot(q, sigmoid(q, *popt), c='k', label='Sigmoid least squares')
ax.plot(q, random_forest.predict(q_), c='r', label='Random forest regression')
ax.plot(q, hgbt_regressor.predict(q_), c='b', label='Histogram gradient boosted regression tree')
ax.legend()

plt.show()


















#mask = qof > np.quantile(qof, 97.5e-2)
#mask = qof > np.quantile(qof, 90e-2)
#mask = qof > 0.5*np.max(qof)
mask = qof > 0.1*np.max(qof)

# Data has to be (N, D) with N = number of samples, D = dimensionality of each sample
# i.e.: N = nTrials, D = nTrianglesSelected
X = predictors[:,mask]

print(f'Data has {X.shape[0]} samples and {X.shape[1]} features')





# initializing the pacmap instance
# Setting n_neighbors to "None" leads to a default choice shown below in "parameter" section
embedding = pacmap.PaCMAP(n_components=3, n_neighbors=None, MN_ratio=0.5, FP_ratio=2.0) 

# fit the data (The index of transformed data corresponds to the index of the original data)
X_transformed = embedding.fit_transform(X, init="pca")

s = 5
# visualize the embedding
fig = plt.figure()
#plt.ion()

ax = fig.add_subplot(2,2,1, projection='3d')
ax.scatter(X_transformed[:, 0], X_transformed[:, 1], X_transformed[:, 2], s=s, c=response_data[response])
ax.set_title("PAC-Map")

ax = fig.add_subplot(2,2,3, projection='3d')
embedding = pacmap.PaCMAP(n_components=2, n_neighbors=None, MN_ratio=0.5, FP_ratio=2.0) 
X_transformed = embedding.fit_transform(X, init="pca")
ax.scatter(X_transformed[:, 0], X_transformed[:, 1], response_data[response], s=s, c=response_data[response])
ax.set_title("PAC-Map (2 comps)")
ax.set_zlabel(response)


pca= PCA(n_components=3)
pca.fit(X)
X_transformed = pca.transform(X)
ax = fig.add_subplot(2,2,2, projection='3d')
ax.scatter(X_transformed[:, 0], X_transformed[:, 1], X_transformed[:, 2], s=s, c=response_data[response])
ax.set_title("PCA")

ax = fig.add_subplot(2,2,4, projection='3d')
pca= PCA(n_components=2)
pca.fit(X)
X_transformed = pca.transform(X)
ax.scatter(X_transformed[:, 0], X_transformed[:, 1], response_data[response], s=s, c=response_data[response])
ax.set_title("PCA (2 comps)")
ax.set_zlabel(response)







plt.show()
#plt.draw()
#while plt.fignum_exists(fig.number):
#    plt.pause(1.)
#    print("Extending figure life")
#    sleep(0.8)
    
#print("Display ended")

y = response_data[response]
#max_depth = 2 # allow trees only up to depth
#min_percent_of_total_samples_for_leaf = 0.3
#max_bins = 16

n_component_values = np.append(np.arange(1,7,1), np.append(np.arange(8, 31, 2), np.arange(32, 81, 4)))

kf = KFold(n_splits=n_folds, shuffle=True)
qof_random_forest = np.zeros((len(n_component_values), n_folds))
qof_hgbt = np.zeros((len(n_component_values), n_folds))


# Use PCA as a preprocessor for model:
for i, n_components in enumerate(n_component_values):
    
    feature_importance_means = np.zeros((n_folds, n_components)) * np.nan
    feature_importance_stds  = np.zeros((n_folds, n_components)) * np.nan
    n_selected_triangles = X.shape[-1]
    pca_weights = np.zeros((n_folds, n_components, n_selected_triangles))
    selected_triangles = np.where(mask)

    for j, (train_index, test_index) in enumerate(kf.split(X)):
        X_train = X[train_index,:]
        y_train = y[train_index]
        X_test = X[test_index,:]
        y_test = y[test_index]

        pca = PCA(n_components=n_components)
        pca.fit(X_train)
        C_train = pca.transform(X_train)
        # Use the PCA fitted on X_train to project X_test!
        C_test = pca.transform(X_test)

        # (a)
        random_forest = RandomForestRegressor(max_depth=max_depth, min_samples_leaf=min_percent_of_total_samples_for_leaf, n_jobs=n_jobs, n_estimators=n_trees)
        random_forest.fit(C_train, y_train)
        R2 = 1 - (np.mean((y_test - random_forest.predict(C_test))**2) / np.var(y_test))
        qof_random_forest[i,j] = R2

        result = permutation_importance(random_forest, C_test, y_test, n_repeats=10, n_jobs=n_jobs)
        
        # TODO: For each number of components, and each k-fold split, collect:
        # > the feature importances (n_features,)
        # > the std of feature importances (n_features,)
        # > the PCA-component weights to the selected triangles (n_components, n_selected_triangles)
        # (n_component_values, n_folds, ...)
        # So that i can then in matlab color the ROI/selected triangles by the PCA-weight (and the feature importance)
        # and so that i can in matlab plot the feature importances and stds, if needed
        feature_importance_means[j,:] = result.importances_mean[:]
        feature_importance_stds[j,:]  = result.importances_std[:]
        pca_weights[j,:,:] = pca.components_
        
        # (b)
        #hgbt_regressor = HistGradientBoostingRegressor(max_depth=max_depth, min_samples_leaf=round(min_percent_of_total_samples_for_leaf * C_train.shape[0]), max_bins=max_bins)
        #hgbt_regressor.fit(C_train, y_train)
        #R2 = 1 - (np.mean((y_test - hgbt_regressor.predict(C_test))**2) / np.var(y_test))
        #qof_hgbt[i,j] = R2
        
    print(f'n_features={n_components}  ->  Random Forest fit:\tmean R² = {np.mean(qof_random_forest[i,:])}')
    print(f'                  HGBT fit:\t\tmean R² = {np.mean(qof_hgbt[i,:])}')

    matlab_dict = {'feature_importance_means': feature_importance_means, 'feature_importance_stds': feature_importance_stds, 'pca_weights': pca_weights, 'selected_triangles': selected_triangles}
    savemat(f'/home/bnplab-admin/TMS_localization/Figures/regression_forest___{sub_id}_{exp_id}_{mesh_id}_{roi_id}_{response}_n_components-{n_components}.mat', matlab_dict)

fig = plt.figure()
#plt.ion()

np.savetxt(f'/home/bnplab-admin/TMS_localization/Figures/regression_forest_all_qof___{sub_id}_{exp_id}_{mesh_id}_{roi_id}_{response}.mat', qof_random_forest)

ax = fig.add_subplot(1,1,1)
ax.errorbar(n_component_values-0.1, np.mean(qof_random_forest, axis=1), np.std(qof_random_forest, axis=1), c='r', label='Random forest')
ax.errorbar(n_component_values+0.1, np.mean(qof_hgbt, axis=1), np.std(qof_hgbt, axis=1), c='b', label='Histogram gradient boosted')
ax.hlines(R2_kfold_sigmoid, 0, np.max(n_component_values), colors='k', linestyles='dashed', label=f"Sigmoid {n_folds}-fold CV")
ax.hlines(R2_sigmoid_direct, 0, np.max(n_component_values), colors='k', linestyles='dotted', label=f"Sigmoid")
ax.set_xlabel('number of PCA components')
ax.set_ylabel('R²')
ax.legend()
plt.show()