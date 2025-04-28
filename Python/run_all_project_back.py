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
import argparse

GLOBAL_START = time.time()

#ROOT = "/home/bnplab-admin/TMS_localization/HighDef"
ROOT = "/mnt/d/HighDef-operate/HighDef"

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

file_handler = logging.FileHandler(f"{ROOT}/batch_backproject-{datetime.now().strftime('%Y-%m-%dT%H%M%S')}.log", encoding='utf-8')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(name)s: %(message)s')
file_handler.setFormatter(formatter)
file_handler.setLevel(logging.DEBUG)
logger.addHandler(file_handler)

parser = argparse.ArgumentParser(description='Projects results for all specified subjects to fsavg')
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

    for i_exp, exp_id in enumerate([v for v in subject.exp.keys() if len(v) > 5]):
        roi_id = f'midlayer_{exp_id.split("-")[-1].split("2")[0].lower()}'
        logger.info(f'Post-evaluation for {subject_id} exp={exp_id} roi={roi_id}')
        
        # [3] Project to average headmodel (only for map-R2L and map-L2R)
        logger.info(f'Projecting back to avg. for {subject_id} exp={exp_id} roi={roi_id}')
        start = time.time()
        print(f'\nProjecting {subject_id} {exp_id} back to average brain')
        run_command(f"python TMS_localization/explore_project_to_avg.py -s {subject_id} -e {exp_id} -r {roi_id} -o CsE_FDI_in_uV -m mesh0")
        run_command(f"python TMS_localization/explore_project_to_avg.py -s {subject_id} -e {exp_id} -r {roi_id} -o SIHIscore_FDI -m mesh0")
        print(f'\nDone projecting {subject_id} {exp_id} back to average: Took {timedelta(seconds=time.time() - start)}\n\n')
        logger.info(f'Completed projecting back to avg. for {subject_id} exp={exp_id} roi={roi_id}, took {timedelta(seconds=time.time() - start)}')

    logger.info(f'Completed {subject_id}')

print(f"\n\nDone. Took {timedelta(seconds=time.time() - GLOBAL_START)}.")




