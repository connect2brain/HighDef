% This scripts runs the HIGHDEF Experiment with conventional coils %

% (1) Initialize HIGHDEF Experiment by calling initialize_experiment
% (2) Define Experimental Design settings
% (4) Define subject details and output location
% (5.1) Hotspot-search LEFT coil
% (5.2) Hotspot-search RIGHT coil
% (6) Enter estimated  RMTs here:
% (7.1) IO Curve LEFT coil
% (7.2) IO Curve RIGHT coil
% (8) Define intensity for mapping from IO Curves LEFT and RIGHT
% (9) SIHI IO Curve using intensities from step (8)
% (10) Mapping: Blocks using intensities from step (8)

% Functions for running I/O curves and Blocks

%% Experimental design thoughts (David):
% With 17 Blocks x 30 SP:60 PP x 3s ITI, we'd get >= 1000 PPs --
% however, this would take at least 1:20

% If we instead use SP:PP of 1:5 (e.g. 15:75), we could get >= 1000 PPs in
% about 1:05 min (plus placing coils etc)

% At SP:PP of 1:8 (e.g. 10:80), we could get to >= 1000 PPs in 1h -- which
% i would prefer.
% I don't think we should go below 10 single pulses per block. At 10:80
% SP:PP, we'd measure 10 SP in 4:30min, which seems reasonable (approximately 1SP every 45s).

%% (1) Initialize HIGHDEF Experiment by calling initialize_experiment %%
    % (1) Reset workspace
    % (2) Set Correct Paths 
    % (3) Setup BOSS Device Settings
    % (4) Setup and initialize MAGIC Toolbox to communicate with stimulators
rng('shuffle')
initialize_experiment


%% (2) Define Experimental Design settings

% How many Blocks?
nBlocksReducedIntensity = 4;
nBlocksDisplacement = 10;
nBlocks = nBlocksReducedIntensity + nBlocksDisplacement;
% aj: evaluate Blocknr. and if Blocks should use different intensities and include rotation

% What is the Ratio between single (test-pusle) and randomly interleaved paired pulses (test-pulse + conditioning pulse)? 
% aj: compare with Literature (i.e. Ferbert et al.), maybe nSinglepulses = 20 would be better
nSinglePulses = 20;
nPairedPulses = 80;

% What is the inter-trial interval (ITI) and ITI-range (jitter for randomized ITI)
ITI = 4; %s
JITTER = 1; %s, plusminus, i.e. Unif([3,5])

% What is the time delay between test-pulse and conditioning-pulse (default for SIHI = 10 ms - approx. 20 ms,  Ferbert et al. 1992) 
PAIRED_DELAY = 10e-3; % 10ms


% (3) Calculate & Print Experimental Design settings as defined in (2)

% Calculate total duration of experiment
expectedDuration = seconds((nSinglePulses + nPairedPulses)*ITI); expectedDuration.Format = 'mm:ss';
stdDuration = seconds(sqrt(((nPairedPulses+nSinglePulses)*((2*JITTER)^2)/12))); stdDuration.Format = 's';
% Variance of uniform distr over [a,b]: ((b-a)^2)/12; here: b-a=2*JITTER;
% Variance of sum of iid RV = sum of variances of iid RV; sqrt(VAR) = sd;

% Print Experimental design settings
fprintf(' Total # paired pulses:    %d\n\n', nBlocks*nPairedPulses)
fprintf(' Estimated time per block: %s min +/- %s\n', expectedDuration, stdDuration)
% Print total duration: expectedDuration.Format = 'hh:mm'
fprintf('         Total duration >= %s min\n\n', nBlocks*expectedDuration)


blockCondition = [repmat({'Reduced Intensity'}, 1, nBlocksReducedIntensity) repmat({'Displacement'}, 1, nBlocksDisplacement)];
blockCondition = blockCondition(randperm(length(blockCondition)));

%% (4) Define subject details and output location

% What is the subject ID?
SUBJECT = arrayfun(@num2str, zeros(1,3));
name = flip(input('Enter subject ID > ', 's')); % for pilots: prefix with p
SUBJECT(1:length(name)) = name;
SUBJECT = flip(SUBJECT);

% Output location
OUT_ = 'Z:/Experimental Data/2023-01 HighDef/Pilots/Conventional';
OUT = sprintf('%s/sub-%s', OUT_, SUBJECT);
mkdir(OUT)

iBlock = 1;

% Display warning before starting Hotspot search
warning('Check  C O I L and T R A C K E R orientation! Orange button should be on the left!')
warning('Check Stimulator settings! W A V E F O R M  needs to be biphasic!')


%% (5.1) Hotspot-search LEFT coil

% initialize left coil (starting amplitude: 65% MSO) 
stimulators.left.arm;
stimulators.left.setAmplitude(65); 

% observe EMG and abort using Ctrl+C
fprintf('\nHotspot search. Observe EMG, change  Press Ctrl+C to abort\n')

% Initializing Hotspot search with randomly triggered pulses (1.7 - 2.7 seconds) to avoid modulation effect
while(true), bd.sendPulse(1), fprintf('.'), pause(1.7+rand()), end

% Suggested method of Hotspot search
% % start: coil held tangentially on the scalp at an angle of 45? to the midline in PA (posterior-anterior) position. % Bashir et al 2013
% % Hotspot: shift coil by 0,5 cm and find highest MEP in 5x 5 cm area
% % Orientation: choose point of highest MEP and rotate 45 degrees left and 45 degrees right (increments of 5-10 degreees)
% % RMT estimated RMT for I/O Curve by consecutively reducing MSO by 1-3% such that approximately 50% elicit a MEP.


%% (5.2) Hotspot-search RIGHT coil

% initialize right coil (starting amplitude: 65% MSO)
stimulators.right.arm;
stimulators.right.setAmplitude(63);

% observe EMG and abort using Ctrl+C
fprintf('\nHotspot search. Press Ctrl+C to abort\n')

% Initializing Hotspot search with randomly triggered pulses (1.7 - 2.7 seconds) to avoid modulation effect
while(true), bd.sendPulse(2), fprintf('.'), pause(1.7+rand()), end


%% (6) Enter estimated  RMTs here:
% Count as MEP iff the amplitude is >= 50 uV
% RMT: 5/10 (i.e. ~half) of pulses yield an MEP over 50 uV

rmtLeft  = input('Enter RMT for  L E F T  coil > ');
rmtRight = input('Enter RMT for  R I G H T  coil > ');

fprintf('\n + + + + + +\n + R  M  T +\n + L = %d%% +\n + R = %d%% +\n + + + + + +\n\n', rmtLeft, rmtRight)


%% (7.1) IO Curve LEFT coil
% IO Curve uses 5 different intensities with range RMT and 140% RMT
maxNumPulsesIO = 30; % nr. of pulses per intensity (default: 5)
[intensities, meps, fig, fitresult_left_coil] = runIOCurve(bd, stimulators.left, 1, rmtLeft, maxNumPulsesIO, 1, 2) 

save([OUT filesep sprintf('sub-%s_mep_io-left_hand.mat', SUBJECT)], 'intensities', 'meps', 'fitresult_left_coil')% save to mat file
print(fig, sprintf('%s/sub-%s_mep_io-left_hand.pdf', OUT, SUBJECT), '-dpdf', '-r0') % print to pdf
fitresult_left_coil


%% (7.2) IO Curve RIGHT coil
% IO Curve uses 5 different intensities with range RMT and 140% RMT
[intensities, meps, fig, fitresult_right_coil] = runIOCurve(bd, stimulators.right, 2, rmtRight, maxNumPulsesIO, 3, 4) 

save([OUT filesep sprintf('sub-%s_mep_io-right_hand.mat', SUBJECT)], 'intensities', 'meps', 'fitresult_right_coil') % save to mat file
print(fig, sprintf('%s/sub-%s_mep_io-right_hand.pdf', OUT, SUBJECT), '-dpdf', '-r0') % print to pdf
fitresult_right_coil


%% (8) Define intensity for mapping from IO Curves LEFT and RIGHT
% %we should use a sigmoidal fit here instead (need like 25 pulses)

% Intensity is defined as inflection point of IO curve (chosen arbitrarily) 
intensityLeft  = round(fitresult_left_coil.s) % intensity conditioning pulse
intensityRight = round(fitresult_right_coil.s) % intensity test pulse

% print stimulation intensities
fprintf('\n Stimulation intensities:\n > L = %d%% \n > R = %d%% \n\n', intensityLeft, intensityRight)


%% (9) SIHI IO Curve using intensities
% % SIHI curve using rmtLeft and intensityRight 
FDIl = 5;
maxNumPairedPulses = 40;
sc = bnp_sihi_curve(bd, rmtLeft,intensityRight, stimulators.left, stimulators.right, FDIl, maxNumPairedPulses, 1, 2);

% save SIHI Curve
save(sprintf('%s/sub-%s_sihi_curve-%s.mat', OUT, SUBJECT, datestr(now,'YYYY-mm-dd_HHMMSS')), 'SUBJECT', 'sc')



%% Addendum (2024-01-16) based on Weise et al.:
% They stimulate at 150% RMT for hotspot-identification
% Thus, we use one high intensity (130 % RMT), and one low intensity (100 %
% RMT -- or less?) to sample the intensity space of both MEP-I/O-curve and SIHI-curve
% sufficiently well
% We sample a little asymmetrically, in that we do more high-intensity
% trials
% Important! Check that these values make sense in the context of the
% MEP-I/O-curve (high intensity should be towards the upper plateau)
% and of the SIHI-curve (low intensity should be towards left upper plateau
% of SIHI-I/O-curve)

highIntensityLeftCoil = round(1.3 * rmtLeft);
lowIntensityLeftCoil = round(0.75 * rmtLeft);

%intensityRight = round(1.2 * rmtRight);
%intensityLeft = round(1.2 * rmtLeft);
%relativeIntensities = [0.9 1.0 1.1 1.2 1.3 1.4];
%intensitiesLeft = rmtLeft .* relativeIntensities;




%% Show next block condition:
condition = blockCondition{iBlock};
fprintf('\n\n Condition = %s\n\n', condition)

%% (10) Mapping: Blocks using intensities from step (8)

if iBlock > length(blockCondition)
    condition = input('\n Adding block. Specify condition (Displacement/Reduced Intensity) > ', 's');
else
    condition = blockCondition{iBlock};
end

if strcmpi(condition, 'Displacement')
    intensityLeft = highIntensityLeftCoil;
elseif strcmpi(condition, 'Reduced Intensity')
    intensityLeft = lowIntensityLeftCoil;
else
    error('Unknown condition: %s', condition)
end


% setup Block
[times, channel, markers] = setupBlock(nSinglePulses, nPairedPulses, ITI, JITTER, PAIRED_DELAY);

filename = sprintf('%s/sub-%s_main_pulse_sequence_block-%d_%s', OUT, SUBJECT, iBlock, datetime('now','Format','yyyy-mm-dd-HHMMSS'));
save([filename '.mat'], ...
    "times", "channel", "markers", "condition") % save to mat
planning = table(times, channel, markers, repmat({condition}, length(times), 1), repmat(intensityLeft, length(times), 1), repmat(intensityRight, length(times), 1), 'VariableNames', {'Time', 'Stimulator', 'Marker', 'Condition', 'Intensity_LCoil_MSO', 'Intensity_RCoil_MSO'});
writetable(planning, [filename '.csv'])

fprintf('Running %s Block %d of %d   (~ %d s)\n', condition, iBlock, nBlocks, times(end))

pause(10) % get ready to map with conditioning pulse!

% Run block


runBlock(bd, times, channel, markers, stimulators.left, stimulators.right, intensityLeft, intensityRight)
% Block counter
iBlock = iBlock + 1;





%% STOP SEQUENCE: Run this to abort the sequence 
bd.stop


%% Add additional Rotation block (will be next block):
addWhere = iBlock; 

rotationOrder = [blockCondition(1:iBlock-1) {'Rotation'} blockCondition(iBlock:end)]
nBlocks = nBlocks+1;



%% FUNCTIONS 






%% setupBlock
function [times, channel, markers] = setupBlock(nSinglePulses, nPairedPulses, ITI, JITTER, PAIRED_DELAY)
%arguments
%    nSinglePulses (1,1) int {mustBePositive};
%    nPairedPulses (1,1) int {mustBePositive};
%end

% Create experiment plan (time_point_marker)
singlePulses = repmat({'single'}, nSinglePulses, 1);
pairedPulses = repmat({'paired'}, nPairedPulses, 1);

conditions = [singlePulses; pairedPulses];
nTrials = size(conditions, 1);
randomTrialOrder = randperm(nTrials);
conditions = conditions(randomTrialOrder);
codeSingleRight = 4; % instead of 0 maybe 4? % In p001, it still was 0
codePairedRight = 1;
codePairedLeft  = 8;

markersRight = [repmat(codeSingleRight, nSinglePulses, 1); repmat(codePairedRight, nPairedPulses, 1)];
markersRight = markersRight(randomTrialOrder);
markersLeft = repmat(codePairedLeft, nPairedPulses, 1);
markers = [markersRight; markersLeft]; % i.e. block [RIGHT; LEFT] -- then use the times generated below to bring these into order

% WARNING: DO NOT USE THIS WAY OF PICKING ITIs!: timesRight = ((1:nTrials) .* ITI) + (rand(1,nTrials) .* JITTER);
% -> This has been used widely in the lab, but yields a TRIANGULAR ITI distribution!
% Instead, uniformly pick the ITIs:
% An ITI is: Base-ITI (e.g. 3s) plusminus JITTER (uniform plusminus), i.e. ITI ~ Unif([2 4])
minITI = ITI - JITTER;
ITIrange = 2*JITTER;
ITIs = minITI + (rand(1,nTrials) .* ITIrange);
timesRight = cumsum(ITIs);
timesConditioning = timesRight(strcmpi(conditions, 'paired')) - PAIRED_DELAY;

channelRight = 2.*ones(nTrials, 1);
channelConditioning = ones(nPairedPulses, 1);
channel = [channelRight; channelConditioning];


% Bring arrays into shared order:
[times, trialOrder] = sort([timesRight'; timesConditioning']);
channel = channel(trialOrder);
markers = markers(trialOrder);
end