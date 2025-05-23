# Which script does what?

### Read raw recordings and format for analysis

 - `conventional_initiation_extract_raw.m`: For single-pulse mapping session (initiation session, yields sub-XXX_R_raw.csv and sub-XXX_L_raw.csv)

 - `conventional_paired_extract_raw.m`: For paired-pulse mapping sessions (ses-2, ses-3, yields sub-XXX_R2L_raw.csv and sub-XXX_L2R_raw.csv)

## Figures

 - `visualize_overview_figure.m`: -> **Figure 1**

 - `plot_3d_result_mapping_joint.m`: -> **Figure 2**

 - `conventional_project_to_surface.m`: Does the Borghetti/Julkunen/Pitkänen projection

 - `compare_SSAM_and_projection.m`: -> **p-value for Projection approach hot vs cold (avg)** and **p-value distances Projection vs SSAM** 
    Compares the distances on the avg. brain between Projection approach (Borghetti, Julkunen) and the source-space association mapping (Weise and Numssen, et al.)

 - `visualize_distances_avg_and_individual.m`: -> **Figure 3** and **p-value for SSAM hot vs cold (ind)**, requires the result of `compare_SSAM_and_projection.m`

 - `visualize_fwhm_focality.m`: -> Table of focality assumption


# For post-hoc effect size in SIHI terms
 - On linux: `cross_compare_spots.py` -- this does the E-field simulations from hotspot on coldspot etc (saved as contrasts.csv)
 - Then run (in R): `cross-effectiveness.R` (this writes `2023-01 HighDef/Results/Evaluation/group_level_spot_cross_relevance_dominant.csv`) -> **Supplementary Figure 1**
 - And for visualization, run: `cross_efficacy_all_spots.m`  -> **Figure 4**


