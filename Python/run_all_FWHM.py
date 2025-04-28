
import os
from os import system as run_command
import pynibs
from datetime import datetime, timedelta
import time
import logging
import argparse

GLOBAL_START = time.time()

ROOT = "/mnt/d/HighDef-operate/HighDef"

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

file_handler = logging.FileHandler(f"{ROOT}/run_FWHM_all-{datetime.now().strftime('%Y-%m-%dT%H%M%S')}.log", encoding='utf-8')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(name)s: %(message)s')
file_handler.setFormatter(formatter)
file_handler.setLevel(logging.DEBUG)
logger.addHandler(file_handler)


parser = argparse.ArgumentParser(description='Does full FWHM analysis for all specified subjects')
parser.add_argument('-s', '--subject', help='id of subject, "all", or comma-separated list (e.g. -s "sub-001, sub-002, sub-005")', required=False, type=str, default="all")
args = parser.parse_args()

if args.subject == "all":
    subjects = sorted([s for s in os.listdir(ROOT) if s.startswith("sub") and not os.path.isfile(s)])
else:
    subjects = [s.strip() for s in args.subject.split(",")]
logger.info(f'Subjects: {subjects}')

for subject_id in subjects:
    fn_subject = f"{ROOT}/{subject_id}/{subject_id}.hdf5"
    subject = pynibs.load_subject(fn_subject)

    # [0:ii]  Check that the subject is complete (i.e. all map-R, map-L, map-R2L, map-L2R are present in the subject-hdf5) -- otherwise skip    
    if not all([e in subject.exp.keys() for e in ["map-L", "map-R", "map-R2L", "map-L2R"]]):
        print(f">>] {subject_id} incomplete: Only found {subject.exp.keys()}, but need all 'map-L', 'map-R', 'map-L2R', 'map-R2L' --- thus skipping this subject")
        logger.warning(f'{subject_id} is incomplete: SKIPPING')
        continue

    for i_exp, exp_id in enumerate(["map-L2R", "map-R2L"]):
        roi_id = f'midlayer_{exp_id.split("-")[-1].split("2")[0].lower()}'
        logger.info(f'Post-evaluation for {subject_id} exp={exp_id} roi={roi_id}')
        
        start = time.time()
        print(f'\nCollecting FWHM for {subject_id} {exp_id}')
        logger.info(f'Collecting SIHIscore_FDI FWHM for {subject_id} exp={exp_id} roi={roi_id}')
        run_command(f"python /home/bnplab-admin/TMS_localization/simulate_FWHM_focality.py -s {subject_id} -e {exp_id} -o SIHIscore_FDI")
        logger.info(f'Completed collecting SIHIscore_FDI FWHM for {subject_id} exp={exp_id} roi={roi_id}')

        logger.info(f'Collecting CsE_FDI_in_uV FWHM for {subject_id} exp={exp_id} roi={roi_id}')
        run_command(f"python /home/bnplab-admin/TMS_localization/simulate_FWHM_focality.py -s {subject_id} -e {exp_id} -o CsE_FDI_in_uV")
        logger.info(f'Completed collecting CsE_FDI_in_uV FWHM for {subject_id} exp={exp_id} roi={roi_id}')
        print(f'\nDone collecting FWHM for {subject_id} {exp_id}: Took {timedelta(seconds=time.time() - start)}\n\n')

    logger.info(f'Completed {subject_id}')

print(f"\n\nDone subjects {subjects}. Took {timedelta(seconds=time.time() - GLOBAL_START)}.")



