
% This script initializes the HIGHDEF Experiment % 
% Requirements: MATLAB 2017, Activated (old) BOSS device, Activated and
% connected NeurOne (needed for BOSS Device to work), Hardware setup 

% (1) Reset workspace
% (2) Set Correct Paths 
% (3) Setup BOSS Device Settings
% (4) Setup and initialize MAGIC Toolbox to communicate with stimulators


%% (1) Reset workspace
close all
clear
clc

rng(now) % initialize the random number generator (needed for randomizing stimulaiton intervals)


%% (2) Set Correct Paths 

% Path to MAGIC toolbox:
addpath(genpath('Z:/Projects/2023-01 HighDef/libraries/MAGIC-toolbox')) 
% Path to bnp for SIHI-curve script
addpath(genpath('Z:/Projects/2023-01 HighDef/libraries/bnp'))
% change to current directory
currentDoc = matlab.desktop.editor.getActive; cd(fileparts(currentDoc.Filename));


%% (3) Setup BOSS Device Settings
bd = bossdevice;
bd.calibration_mode = 'no';
bd.armed = 'no';

% Disable eye blink detection
setparam(bd.tg, 'QLY', 'eye_artifact_threshold', 1e6)

% Disable IF stability thresholds
% IF = instantaneous frequency, i.e. the derivative of instantaneous phase
%this undoes an experimental parameter Christoph used for theta detection.
%this turns it off for our purposes.
max_instability = 1e6;
setparam(bd.tg, 'QLY', 'inst_freq_max_instability', [1e6 1e6 1e6 1e6 1e6 1e6])

% Define channels to be used
bd.aux_channels = 6; % 2 emg channels each side %aj: should this be 8?
bd.eeg_channels = 10; % C3 and C4 Hjorth channels. This needs to be correctly configured in the NeurOne

% Define range of interest for EMG data (Default: EMG data from -100 to
% +100 ms)
bd.scope_emg.NumSamples = 1000;
bd.scope_emg.NumPrePostSamples = -500;

% Set EEG artifact detection threshold
eeg_artifact_threshold = 1e6;
setparam(bd.tg, 'QLY', 'eeg_artifact_threshold', [eeg_artifact_threshold eeg_artifact_threshold])



%% (4) Setup and initialize MAGIC Toolbox to communicate with stimulators
try
    stimulators = [];
    stimulators.left = magventure('COM2');
    stimulators.right = magventure('COM1');
    stimulators.left.connect();
    stimulators.right.connect();
    stimulators.left.setTrig_ChargeDelay(0,0,1000);
catch
    fprintf(2,['Problem with the MAGIC toolbox, please restart Matlab 2017.\n'...
        '(can be run only once for each Matlab session.)\n']);
end
