import subprocess
import numpy as np
import pandas as pd
import nibabel as nib
from scipy.spatial import KDTree

# TODO: use correct subject

hemisphere = "dominant"
sub_id = "sub-014"

fs_surface = "inflated"

if sub_id != "sub-006":
    if hemisphere == "dominant":
        pial_name = f"lh.{fs_surface}"
        exp_id = "L2R"
    else:
        pial_name = f"rh.{fs_surface}"
        exp_id = "R2L"
else:
    if hemisphere == "dominant":
        pial_name = f"rh.{fs_surface}"
        exp_id = "R2L"
    else:
        pial_name = f"lh.{fs_surface}"
        exp_id = "L2R"

result = subprocess.run(f"mri_info --vox2ras TMS_localization/HighDef/{sub_id}/mesh/T1.nii.gz", capture_output=True, shell=True)
vox2ras_matrix = np.array([[float(v) for v in s.split()] for s in result.stdout.decode("utf-8").split("\n") if s != ''])
print("vox2ras matrix:")
print(vox2ras_matrix)

df = pd.read_csv(f"TMS_localization/HighDef/{sub_id}/{sub_id}_{exp_id}_raw.csv", delimiter=",")

locations = np.vstack([df["p1"].to_numpy(), df["p2"].to_numpy(), df["p3"].to_numpy(), np.ones((df.shape[0],))])
transformed_locations = locations #vox2ras_matrix @ locations

print(locations.shape, "->", transformed_locations.shape)

fs_subject_dir = '/usr/local/freesurfer/7.4.1/subjects'

pial_vertices, pial_faces = nib.freesurfer.read_geometry(f"{fs_subject_dir}/fsaverage/surf/{pial_name}")
pial_tree = KDTree(pial_vertices)

distances, indices = pial_tree.query(transformed_locations[0:-1,:].T)
closest_pial_points = pial_vertices[indices]

print(closest_pial_points.shape)


import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(projection="3d")
ax.plot_trisurf(pial_vertices[:,0], pial_vertices[:,1], pial_vertices[:,2], triangles=pial_faces)
ax.scatter(transformed_locations[0,:], transformed_locations[1,:], transformed_locations[2,:], marker="o")
ax.scatter(closest_pial_points[:,0], closest_pial_points[:,1], closest_pial_points[:,2], marker="^")
plt.show()
