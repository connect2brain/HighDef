# TODO: In this script, run all the evaluation and analysis steps on all complete subjects
# (1) Create headmodel [optional: only if not yet done] and Regions of interest [obligatory!]
# (2) For all experiments:
#   (2.1) E field simulation
#   (2.2) Regression
#   (2.3) Optimal coil position
# (3) Project FDI hotspot and coldspot onto average brain, for map-R2L and map-L2R
# (4) Run FWHM analysis for FDI hotspot and coldspot (for map-R2L and map-L2R)

import os
from os import system as run_command
import pynibs
from datetime import datetime, timedelta
import time
import logging

GLOBAL_START = time.time()

ROOT = "/home/bnplab-admin/TMS_localization/HighDef"

logger = logging.getLogger(__name__)
logging.basicConfig(filename=f"{ROOT}/run_focality-{datetime.now().strftime('%Y-%m-%dT%H%M%S')}.log", encoding='utf-8', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

subjects = sorted([s for s in os.listdir(ROOT) if s.startswith("sub") and not os.path.isfile(s)])
logger.info(f'Subjects: {subjects}')

for subject_id in subjects:
    fn_subject = f"{ROOT}/{subject_id}/{subject_id}.hdf5"
    subject = pynibs.load_subject(fn_subject)
    logger.info(f"Started {subject_id}")

    # [0:ii]  Check that the subject is complete (i.e. all map-R, map-L, map-R2L, map-L2R are present in the subject-hdf5) -- otherwise skip    
    if not all([e in subject.exp.keys() for e in ["map-L", "map-R", "map-R2L", "map-L2R"]]):
        print(f">>] {subject_id} incomplete: Only found {subject.exp.keys()}, but need all 'map-L', 'map-R', 'map-L2R', 'map-R2L' --- thus skipping this subject")
        logger.warning(f'{subject_id} is incomplete: SKIPPING')
        continue

    for i_exp, exp_id in enumerate(["map-L2R", "map-R2L"]):
        roi_id = f'midlayer_{exp_id.split("-")[-1].split("2")[0].lower()}'
        logger.info(f'Post-evaluation for {subject_id} exp={exp_id} roi={roi_id}')
        
        # [4] Run FWHM
        start = time.time()
        print(f'\nCollecting FWHM for {subject_id} {exp_id}')
        logger.info(f'Collecting SIHIscore_FDI FWHM for {subject_id} exp={exp_id} roi={roi_id}')
        run_command(f"python TMS_localization/simulate_FWHM_focality.py -s {subject_id} -e {exp_id} -o SIHIscore_FDI")
        logger.info(f'Completed collecting SIHIscore_FDI FWHM for {subject_id} exp={exp_id} roi={roi_id}')

        logger.info(f'Collecting CsE_FDI_in_uV FWHM for {subject_id} exp={exp_id} roi={roi_id}')
        run_command(f"python TMS_localization/simulate_FWHM_focality.py -s {subject_id} -e {exp_id} -o CsE_FDI_in_uV")
        logger.info(f'Completed collecting CsE_FDI_in_uV FWHM for {subject_id} exp={exp_id} roi={roi_id}')
        print(f'\nDone collecting FWHM for {subject_id} {exp_id}: Took {timedelta(seconds=time.time() - start)}\n\n')

    logger.info(f'Completed {subject_id}')

logger.info(f"\n\nDone. Took {timedelta(seconds=time.time() - GLOBAL_START)}.")
print(f"\n\nDone. Took {timedelta(seconds=time.time() - GLOBAL_START)}.")




