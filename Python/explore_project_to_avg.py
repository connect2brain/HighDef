
import numpy as np
import nibabel as nib
import h5py
import matplotlib.pyplot as plt
import h5py
import subprocess
import argparse
from time import sleep


parser = argparse.ArgumentParser(description='Projects data from given subject, experiment etc back to freesurfer inflated average brain')
parser.add_argument('-s', '--subject', help='ID (!) of subject --- not its path!', required=True, type=str)
parser.add_argument('-e', '--exp_id', help='id of the experiment (map-L2R or map-R2L)', required=True, type=str)
parser.add_argument('-m', '--mesh_id', help='Mesh ID', required=True, type=str)
parser.add_argument('-r', '--roi_id', help='ROI ID (midlayer_l/midlayer_r or small_l/small_r)', required=True, type=str)
parser.add_argument('-o', '--response', help='Response (either CsE_FDI_in_uV or SIHIscore_FDI)', required=False, type=str, default='SIHIscore_FDI')
parser.add_argument('-p', '--part', help="Which split/part to perform the back-projection on. Default: Full experiment", required=False, default=None, type=str)
args = parser.parse_args()

sub_id = args.subject
exp_id = args.exp_id
mesh_id = args.mesh_id
roi_id = args.roi_id
response = args.response
split = args.part

suffix = ""
if split is not None:
    suffix = f"_{split}"

#sub_id = "sub-001"
base_path = f"/mnt/d/HighDef-operate/HighDef/{sub_id}"
#exp_id = "map-L2R"
#roi_id = "small_l"
#response = "CsE_FDI_in_uV"
#response = "SIHIscore_FDI"

print(f"{'_'.join(20*['_'])}\n\n Projecting {sub_id} back to average brain\n Experiment: {exp_id}\n Mesh:       {mesh_id}\n ROI:        {roi_id}\n Response:   {response}\n")

individual_surface = nib.load(f'{base_path}/mesh/m2m_{sub_id}/surfaces/{"l" if exp_id == "map-L2R" else "r"}h.central.gii')
individual_vertices = individual_surface.darrays[0].data
individual_faces = individual_surface.darrays[1].data

print(f"Individual freesurfer surface has {individual_vertices.shape[0]} vertices and {individual_faces.shape[0]} faces")

# Important question: Does this align directly with the data from the hdf5?
# i.e 
# (1) same vertices? -- total number will be different, bc. hdf5 is only ROI!
# (2) same triangles? -- this should be the case

# Basically check: Are all vertices of the ROI-mesh (hdf5) also present
# Load hdf5
with h5py.File(f'{base_path}/mesh/roi/{roi_id}/geo.hdf5') as file:
    h5_vertices  = file['/mesh/nodes/node_coord'][:]
    h5_triangles = file['/mesh/elm/triangle_number_list'][:]

print(f"ROI from hdf5 has {h5_vertices.shape[0]} vertices and {h5_triangles.shape[0]} faces")

# As far as i can see, they in fact coincide!
# (h5_triangles == individual_faces[:,None]).all(-1).any() 
# ^ Apparently (https://stackoverflow.com/questions/64219670/check-if-any-row-in-a-numpy-array-is-part-of-another-array), this checks it.

# So now:
# Load the data for the ROI:
if split is not None:
    with h5py.File(f'{base_path}/results/exp_{exp_id}/r2/mesh_mesh0/roi_{roi_id}/{response}/sigmoid4/{split}/r2_roi_data.hdf5') as file:
        h5_data = file['/data/tris/c_E_mag'][:]
else:
    with h5py.File(f'{base_path}/results/exp_{exp_id}/r2/mesh_mesh0/roi_{roi_id}/{response}/sigmoid4/r2_roi_data.hdf5') as file:
        h5_data = file['/data/tris/c_E_mag'][:]

print(f"There are {h5_data.shape[0]} datapoints")

#distance = 1e-3
#for i_vertex, vertex in enumerate(h5_vertices):
#    if np.min(np.sqrt(np.sum((individual_vertices - vertex)**2, axis=-1))) > distance:
#        print(f"Vertex {i_vertex} (@{vertex[0]}, {vertex[1]}, {vertex[2]}) has no match within {distance}")
# ^ This did not find any mismatches!


# Each data point is associated with the corresponding triangle:
# h5_triangles[i,:] <- h5_data[i]
# And: h5_triangles[i,:] contains also directly the indices of the correct individual_vertices to which to assign these data
# From the botched result, it seems the indices don't actually correspond perfectly.
# Use nearest neighbor.

vertex_maxdata = np.zeros((individual_vertices.shape[0],))
vertex_sumdata = np.zeros((individual_vertices.shape[0],))
vertex_counts = np.zeros((individual_vertices.shape[0],))

n_tris = h5_triangles.shape[0]
for i_triangle in range(n_tris):
    val = h5_data[i_triangle]
    for vertex in h5_triangles[i_triangle,:]:
        vertex_coordinates = h5_vertices[vertex]
        which_is_that_in_the_gii = np.argmin(np.sqrt(np.sum((individual_vertices - vertex_coordinates)**2, axis=-1)))
        vertex_maxdata[which_is_that_in_the_gii]  = max(vertex_maxdata[which_is_that_in_the_gii], val)
        vertex_sumdata[which_is_that_in_the_gii] += val
        vertex_counts[which_is_that_in_the_gii]  += 1
    if (100*(i_triangle / n_tris)) % 10 == 9:
        print(f"  {i_triangle} / {n_tris}")

vertex_counts[vertex_counts == 0] = 1
vertex_meandata = vertex_sumdata / vertex_counts

# Save to mgh file.
image_maxdata = vertex_maxdata[:, np.newaxis, np.newaxis].astype(np.float32)
mgh_img = nib.MGHImage(image_maxdata, affine=None)
output_path = f'{base_path}/{sub_id}_{exp_id}_{response}{suffix}_vertmax.mgh'
nib.save(mgh_img, output_path)

image_meandata = vertex_meandata[:, np.newaxis, np.newaxis].astype(np.float32)
mgh_img = nib.MGHImage(image_meandata, affine=None)
output_path = f'{base_path}/{sub_id}_{exp_id}_{response}{suffix}_vertavg.mgh'
nib.save(mgh_img, output_path)


# Then run mri_surf2surf:

bash_script_path = f"{base_path}/project_back_{response}{suffix}.sh"
with open(bash_script_path, "w") as f:
      target_path = f"{base_path}/mesh"
      f.write("export LD_LIBRARY_PATH=/home/bnplab-admin/miniconda3/envs/tms_loco/lib:$LD_LIBRARY_PATH\n")
      f.write(f"export SUBJECTS_DIR='{target_path}'\n")
      f.write(f"cd {target_path}\n")
      f.write(f"mri_surf2surf --srcsubject m2m_{sub_id} "
              f"--srcsurfval '{base_path}/{sub_id}_{exp_id}_{response}{suffix}_vertmax.mgh' "
              f"--trgsurfval '{base_path}/projected_{sub_id}_{exp_id}_{response}{suffix}_vertmax.mgh' "
              f"--hemi {'l' if exp_id == 'map-L2R' else 'r'}h "
              f"--trgsubject fsaverage --trgsurfreg sphere.reg --srcsurfreg sphere.reg.gii\n")

sleep(0.3)
print(f"\nWrote commands to {bash_script_path}; running it now")
subprocess.run(['bash', bash_script_path])


