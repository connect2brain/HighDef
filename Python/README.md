# HighDef Python code

This code needs to be run on a Linux PC or windows subsystem

TODO: Add instructions for how to set up WSL subsystem


## NEW commands:

```python Installation/tmsloc_proto-0.2023.8/scripts/01_create_subject_structure.py -f /mnt/d/HighDef-operate/HighDef/sub-999```

Then modify the `create_sub-XXX.py` file, putting in the correct MRI names etc.

```python /mnt/d/HighDef-operate/HighDef/sub-001/create_sub-001.py```

```python Installation/tmsloc_proto-0.2023.8/scripts/02_make_msh_from_mri_simnibs4.py -s /mnt/d/HighDef-operate/HighDef/sub-001/sub-001.hdf5 -m mesh0 --charm_ini /mnt/c/Users/bnplab-admin/SimNIBS-4.0/simnibs_env/Lib/site-packages/simnibs/charm.ini --charm_cmd_params forceqform```

```python TMS_localization/evaluate.py -f /mnt/d/HighDef-operate/HighDef/sub-001/sub-001.hdf5 -e map-L2R```


## OLD Commands:
These commands assume that the directory Installation contains the Weise-pipeline code as tmsloc_proto-0.2023.8

```python Installation/tmsloc_proto-0.2023.8/scripts/01_create_subject_structure.py -f TMS_localization/HighDef/sub-001```

```python TMS_localization/HighDef/sub-001/create_sub-001.py```

```python Installation/tmsloc_proto-0.2023.8/scripts/02_make_msh_from_mri_simnibs4.py -s TMS_localization/HighDef/sub-001/sub-001.hdf5 -m mesh0 --charm_ini /mnt/c/Users/bnplab-admin/SimNIBS-4.0/simnibs_env/Lib/site-packages/simnibs/charm.ini --charm_cmd_params forceqform```

```python TMS_localization/0506_custom.py -f TMS_localization/HighDef/sub-001/sub-001.hdf5 -e map-L2R -m mesh0 -r midlayer_l```

```python TMS_localization/07_custom.py -f TMS_localization/HighDef/sub-001/sub-001.hdf5 -e map-R2L -m mesh0 -r midlayer_r```

```python TMS_localization/08_calc_opt_coil_pos.py -c '/mnt/c/Users/bnplab-admin/SimNIBS-4.0/simnibs_env/Lib/site-packages/simnibs/resources/coil_models/Drakaki_BrainStim_2022/MagVenture_Cool-B35.ccd' -s TMS_localization/HighDef/sub-002/sub-002.hdf5 -m mesh0 -e map-R -n 25 -q "E_mag" -t TMS_localization/HighDef/sub-002/results/exp_map-R/r2/mesh_mesh0/roi_midlayer_r/CsE_FDI_in_uV/sigmoid4/r2_roi_data.hdf5 -a "scalar"```

## Installations:
### To install the MCR (Matlab):
sudo FREESURFER_HOME=$FREESURFER_HOME /usr/local/freesurfer/7.4.1/bin/fs_install_mcr R2019b