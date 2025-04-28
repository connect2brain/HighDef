# This script just serves to produce a visualization of the coil placed according to a single trial

from simnibs import sim_struct, run_simnibs, mni2subject_coords
from datetime import datetime
from scipy.spatial.transform import Rotation
import numpy as np
import pandas as pd

ID = 'sub-003'

IN_BASE = '/home/bnplab-admin/TMS_localization/HighDef'
IN = f'{IN_BASE}/{ID}'
OUT = f'/home/bnplab-admin/TMS_localization/TMS_loc_results/sim-{datetime.now().strftime("%Y-%m-%dT%H%M%S")}'


# Initalize a session
s = sim_struct.SESSION()
# Name of head mesh
s.subpath = f'{IN}/mesh/m2m_{ID}'
# Output folder
s.pathfem = OUT


# Initialize a list of TMS simulations
tmslist = s.add_tmslist()
# Select coil
tmslist.fnamecoil = '/mnt/c/Users/bnplab-admin/SimNIBS-4.0/simnibs_env/Lib/site-packages/simnibs/resources/coil_models/Drakaki_BrainStim_2022/MagVenture_Cool-B35.ccd'
#tmslist.fnamecoil = 'C:/Users/David Admin/SimNIBS-4.0/simnibs_env/Lib/site-packages/simnibs/resources/coil_models/Drakaki_BrainStim_2022/MagVenture_Cool-B35.ccd'


## Initialize a coil position
#pos = tmslist.add_position()
## Select coil centre
#pos.centre = 'C3'
## Select coil direction
#pos.pos_ydir = 'FC1'



selected_trial = 47 - 1 # 47 in Matlab


transforms_file = f'{IN}/{ID}_L2R_raw.csv'
raw_data = pd.read_csv(transforms_file)

pos = tmslist.add_position()
    
M = np.eye(4)
for j, colname in enumerate(["x", "y", "z", "p"]):
    for i in range(3):
        M[i,j] = raw_data[f"{colname}{i+1}"][selected_trial]
        
# Localite to SimNIBS:
# x <-> z; -y
M_simnibs = M
M_simnibs[:,0] =  M[:,2]
M_simnibs[:,2] =  M[:,0]
M_simnibs[:,1] = -M[:,1]
pos.matsimnibs = M_simnibs

print(f"Selected trial {selected_trial}:\n\tADM = {raw_data['CsE_ADM_in_uV'][selected_trial]}\n\tFDI = {raw_data['CsE_FDI_in_uV'][selected_trial]}\n\tAPB = {raw_data['CsE_ADM_in_uV'][selected_trial]}\n")


s.open_in_gmsh = False

# Big PC:
#run_simnibs(s, cpus=16)
run_simnibs(s)
