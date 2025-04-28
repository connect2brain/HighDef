import os
import pandas as pd
import numpy as np
import argparse
import h5py
import pynibs
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import dendrogram
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from fastcluster import linkage
from sklearn.cluster import AgglomerativeClustering

# From: https://scikit-learn.org/stable/auto_examples/cluster/plot_agglomerative_dendrogram.html
def plot_dendrogram(model, **kwargs):
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack(
        [model.children_, model.distances_, counts]
    ).astype(float)

    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, **kwargs)


# From: https://gmarti.gitlab.io/ml/2017/09/07/how-to-sort-distance-matrix.html

def seriation(Z,N,cur_index):
    '''
        input:
            - Z is a hierarchical tree (dendrogram)
            - N is the number of points given to the clustering process
            - cur_index is the position in the tree for the recursive traversal
        output:
            - order implied by the hierarchical tree Z
            
        seriation computes the order implied by a hierarchical tree (dendrogram)
    '''
    if cur_index < N:
        return [cur_index]
    else:
        left = int(Z[cur_index-N,0])
        right = int(Z[cur_index-N,1])
        return (seriation(Z,N,left) + seriation(Z,N,right))
    
def compute_serial_matrix(dist_mat,method="ward"):
    '''
        input:
            - dist_mat is a distance matrix
            - method = ["ward","single","average","complete"]
        output:
            - seriated_dist is the input dist_mat,
              but with re-ordered rows and columns
              according to the seriation, i.e. the
              order implied by the hierarchical tree
            - res_order is the order implied by
              the hierarhical tree
            - res_linkage is the hierarhical tree (dendrogram)
        
        compute_serial_matrix transforms a distance matrix into 
        a sorted distance matrix according to the order implied 
        by the hierarchical tree (dendrogram)
    '''
    N = len(dist_mat)
    flat_dist_mat = squareform(dist_mat)
    res_linkage = linkage(flat_dist_mat, method=method,preserve_input=True)
    res_order = seriation(res_linkage, N, N + N-2)
    seriated_dist = np.zeros((N,N))
    a,b = np.triu_indices(N,k=1)
    seriated_dist[a,b] = dist_mat[ [res_order[i] for i in a], [res_order[j] for j in b]]
    seriated_dist[b,a] = seriated_dist[a,b]
    
    return seriated_dist, res_order, res_linkage





parser = argparse.ArgumentParser(description='Creates a new subject in path')
parser.add_argument('-f', '--fn_subject', help='path of new subject', required=True, type=str)
parser.add_argument('-e', '--exp_id', help='id of the experiment', required=True, type=str)
parser.add_argument('-m', '--mesh_id', help='Mesh ID', required=True, type=str)
parser.add_argument('-r', '--roi_id', help='ROI ID', required=True, type=str)
args = parser.parse_args()

fn_subject = os.path.abspath(args.fn_subject)
subject_id = os.path.split(fn_subject)[1]
exp_id = args.exp_id
mesh_id = args.mesh_id
roi_id = args.roi_id


subject = pynibs.load_subject(fn_subject)
subject_dir = subject.subject_folder
fn_raw_data_table = subject.exp[exp_id]['response_csv']

raw_data = pd.read_csv(os.path.join(subject_dir, fn_raw_data_table))
fn_e_roi = os.path.join(subject_dir, "results", f"exp_{exp_id}", "electric_field",
                            f"mesh_{mesh_id}", f"roi_{roi_id}", f"e.hdf5")
with h5py.File(fn_e_roi, "r") as f:
    E_matrix = f['E_mag'][:]

print('Shape of E_matrix: ', E_matrix.shape)

high_SI_trials = raw_data['Intensity_percentMSO'] == np.max(raw_data['Intensity_percentMSO'])
high_SI_data = raw_data.loc[high_SI_trials,'SIHIscore_FDI']
print(high_SI_data)

cutoff = np.percentile(high_SI_data, 75)
selected_trials = high_SI_data > cutoff

E_selected = E_matrix[high_SI_trials, :]
E_selected = E_selected[selected_trials, :]
print('Shape of E_selected: ', E_selected.shape)

N = np.sum(selected_trials)
dist_mat = squareform(pdist(E_selected, 'correlation'))
print('Shape of dist_mat: ', dist_mat.shape)



model = AgglomerativeClustering(n_clusters=None, affinity='precomputed', distance_threshold=0, linkage='complete')

model = model.fit(dist_mat)
plt.title("Hierarchical Clustering Dendrogram")
# plot the top three levels of the dendrogram
plot_dendrogram(model, truncate_mode="level", p=3)
plt.xlabel("Number of points in node (or index of point if no parenthesis).")
plt.show()









# divnorm=colors.TwoSlopeNorm(vmin=0., vcenter=1., vmax=2.)  # np.max(dist_mat)

# methods = ["ward","single","average","complete"]
# for method in methods:
#     print("Method:\t",method)
    
#     ordered_dist_mat, res_order, res_linkage = compute_serial_matrix(dist_mat,method)
    
#     plt.pcolormesh(ordered_dist_mat, norm=divnorm, cmap='RdGy_r') # RdGy_r or Spectral_r
#     plt.xlim([0,N])
#     plt.ylim([0,N])
#     plt.title(method)
#     plt.colorbar()
#     plt.show()