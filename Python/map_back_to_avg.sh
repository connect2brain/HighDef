export SUBJECTS_DIR='/home/bnplab-admin/TMS_localization/HighDef/sub-001/mesh'
cd /home/bnplab-admin/TMS_localization/HighDef/sub-001/mesh
mri_surf2surf --srcsubject m2m_sub-001 --srcsurfval '/home/bnplab-admin/TMS_localization/HighDef/sub-001/mesh/roi/small_l/mask_leftmotor.mgh' --trgsurfval '/home/bnplab-admin/TMS_localization/HighDef/sub-001/trying_out.mgh' --hemi lh --trgsubject fsaverage --trgsurfreg sphere.reg --srcsurfreg sphere.reg.gii
