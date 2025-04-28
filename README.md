# HighDef
Code base to accompany publication "Short-interval interhemispheric inhibition does not originate from the motor hotspot"

## Overview
The code is sorted by programming language. Two systems were used to run the analysis: 
 - The python code was run on a linux subsystem (installation see below)
 - the Matlab and R code meanwhile under windows

The Matlab code is split into:
 - `Experiment`: The code used to run the experiment in a setup with neurone, localite and bossdevice
 - `Analysis`: The code for analysis before and after the analysis done in python on the linux subsystem. Includes plotting scripts


Further details are given in the respective README.md files in the sub-folders



# How to install the ubuntu subsystem

## Basis:

Run as admin with windows command line:
```
wsl --install -d Ubuntu-22.04
```


**R E B O O T !**



After reboot, check that the distribution runs on WSL 2: 
```
wsl -l -v
```
If not, run:
```
wsl --set-version Ubuntu-22.04 2
```

Start wsl, set username and password (if not already set during installation)

Then run:
```
sudo apt-get update
sudo apt-get upgrade
```
	
Install anaconda (or miniconda; see https://www.anaconda.com/docs/getting-started/miniconda/install#linux):
```
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh

source ~/miniconda3/bin/activate

conda init --all
```

Follow the instructions these commands print (such as closing and opening the cmd window again).
	
	
	
	
Install Weise-Numssen-pipeline as:
```
conda env update -f https://github.com/simnibs/simnibs/releases/download/v4.1.0/environment_linux.yml --name tms_loco
conda activate tms_loco
pip install https://github.com/simnibs/simnibs/releases/download/v4.1.0/simnibs-4.1.0-cp39-cp39-linux_x86_64.whl
pip install pynibs
```

Download Weise-Numssen release 0.2023.8 from https://gitlab.gwdg.de/tms-localization/papers/tmsloc_proto/-/releases
and put extract it into an "Installation" folder, under (~/Installation); this gives access to the scripts


## Freesurfer

Download freesurfer 7.4.1, and follow https://surfer.nmr.mgh.harvard.edu/fswiki/LinuxInstall (see also: https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall5.3)
That is, first run:
```
wget https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/7.4.1/freesurfer-linux-ubuntu22_amd64-7.4.1.tar.gz
sudo tar -C /usr/local -xzvf freesurfer-linux-ubuntu22_amd64-7.4.1.tar.gz
```
then add the following to your .bashrc file (located in home/USERNAME/.bashrc):
```
export FREESURFER_HOME=/usr/local/freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh
```
Copy your freesurfer license key into the FREESURFER_HOME directory

Close the WSL window and open it again. This should now display a bunch of information about your freesurfer installation (such as FREESURFER_HOME)

You can now delete the freesurfer tar file


## VS Code (optional; for developing)
Install VS Code on your PC (outside linux subsystem), if you don't have it already
Add the following extensions:
 - Pylance
 - Python
 - Python Debugger
 - WSL

The WSL extension allows you to operate the scripts on the WSL subsystem.


