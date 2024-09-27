rng('shuffle')
initialize_experiment
sessionName = struct(); sessionName.R = 'R2L'; sessionName.L = 'L2R';

%% (2) Define Experimental Design settings

% How many Blocks?
nBlocksReducedIntensity = 2;
nBlocksDisplacement = 10;
nBlocks = nBlocksReducedIntensity + nBlocksDisplacement;

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
name = flip(input('Enter subject ID > sub-', 's')); % for pilots: prefix with p
SUBJECT(1:length(name)) = name;
SUBJECT = flip(SUBJECT);
ID = sprintf('sub-%s', SUBJECT);


T = readtable('//storage.neurologie.uni-tuebingen.de/bbnp_lab/Experimental Data/2023-01 HighDef/Conventional/Planning.xlsx');
subjectRow = T(strcmpi(T.Subject, ID), :);

whichSession = input('Which mapping session [1,2] > ');
assert(ismember(whichSession, [1, 2]))
if whichSession == 1
    map_hemisphere = subjectRow.x1st_bilateral{:};
else
    map_hemisphere = subjectRow.x2nd_bilateral{:};
end


% Output location
OUT_ = 'Z:/Experimental Data/2023-01 HighDef/Conventional';
OUT = sprintf('%s/sub-%s/map-%s', OUT_, SUBJECT, sessionName.(map_hemisphere));
mkdir(OUT)

iBlock = 1;



assert(size(subjectRow, 1) == 1)


% Display warning before starting Hotspot search
warning('Check  C O I L and T R A C K E R orientation! Orange button should be on the left!')
warning('Check Stimulator settings! W A V E F O R M  needs to be biphasic!')


%% Search/Check Left hotspot
searchHotspot(stimulators.left, bd, 1)

%% Search/Check Right hotspot
searchHotspot(stimulators.right, bd, 2)

%% Enter estimated  RMTs here:
% Count as MEP iff the amplitude is >= 50 uV
% RMT: 5/10 (i.e. ~half) of pulses yield an MEP over 50 uV
rmtLeft  = input('Enter RMT for  L E F T  coil > ');
rmtRight = input('Enter RMT for  R I G H T  coil > ');
fprintf('\n[RMT]\n      L = %d %% MSO\n      R = %d %% MSO\n\n', rmtLeft, rmtRight)



%% (7.1) IO Curve LEFT coil
maxNumPulsesIO = 30;
[intensities, meps, fig, fitresult_left_coil] = runIOCurve(bd, stimulators.left, 1, rmtLeft, maxNumPulsesIO, 1, 2);

save([OUT filesep sprintf('%s_mep_io-left_hand.mat', ID)], 'intensities', 'meps', 'fitresult_left_coil')% save to mat file
print(fig, sprintf('%s/%s_mep_io-left_hand.pdf', OUT, ID), '-dpdf', '-r0') % print to pdf
fitresult_left_coil


%% (7.2) IO Curve RIGHT coil
[intensities, meps, fig, fitresult_right_coil] = runIOCurve(bd, stimulators.right, 2, rmtRight, maxNumPulsesIO, 3, 4);

save([OUT filesep sprintf('%s_mep_io-right_hand.mat', ID)], 'intensities', 'meps', 'fitresult_right_coil') % save to mat file
print(fig, sprintf('%s/%s_mep_io-right_hand.pdf', OUT, ID), '-dpdf', '-r0') % print to pdf
fitresult_right_coil

%% (8) Define intensity for mapping from IO Curves LEFT and RIGHT
shift_up_TS = 5;
if strcmpi(map_hemisphere, 'L')
    rmtCS = rmtLeft;
    intensityTS = round(fitresult_right_coil.s) + shift_up_TS;
    CSstimulator = stimulators.left;
    TSstimulator = stimulators.right;
    inhibitedFDI = 4;
    portCS = 1;
    portTS = 2;
else
    rmtCS = rmtRight;
    intensityTS = round(fitresult_left_coil.s) + shift_up_TS;
    CSstimulator = stimulators.right;
    TSstimulator = stimulators.left;
    inhibitedFDI = 2;
    portCS = 2;
    portTS = 1;
end

highIntensityCS = min([round(1.5 * rmtCS) 100]);
lowIntensityCS = round(0.75 * rmtCS);

% print stimulation intensities
fprintf('\n Mapping %s hemisphere with:\n > CS = %d %% / %d %% \n > TS = %d%% \n\n', map_hemisphere, lowIntensityCS, highIntensityCS, intensityTS)


%% (9) SIHI IO Curve using these intensities

maxNumPairedPulses = 40;
sc = bnp_sihi_curve(bd, rmtCS, intensityTS, CSstimulator, TSstimulator, inhibitedFDI, maxNumPairedPulses, portCS, portTS);

% save SIHI Curve
save(sprintf('%s/sub-%s_sihi_curve-%s.mat', OUT, SUBJECT, datestr(now,'YYYY-mm-dd_HHMMSS')), 'SUBJECT', 'sc')



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
    intensityCS = highIntensityCS;
elseif strcmpi(condition, 'Reduced Intensity')
    intensityCS = lowIntensityCS;
else
    error('Unknown condition: %s', condition)
end



% setup Block
[times, channel, markers] = setupBlock(nSinglePulses, nPairedPulses, ITI, JITTER, PAIRED_DELAY, portCS, portTS);

filename = sprintf('%s/sub-%s_main_pulse_sequence_block-%d_%s', OUT, SUBJECT, iBlock, datetime('now','Format','yyyy-mm-dd-HHMMSS'));
save([filename '.mat'], ...
    "times", "channel", "markers", "condition") % save to mat
planning = table(times, channel, markers, repmat({condition}, length(times), 1), repmat(intensityCS, length(times), 1), repmat(intensityTS, length(times), 1), 'VariableNames', {'Time', 'Stimulator', 'Marker', 'Condition', 'Intensity_CS_MSO', 'Intensity_TS_MSO'});
writetable(planning, [filename '.csv'])

fprintf('Running %s Block %d of %d   (~ %d s)\n', condition, iBlock, nBlocks, times(end))

pause(10) % get ready to map with conditioning pulse!

% Run block


runBlock(bd, times, channel, markers, CSstimulator, TSstimulator, intensityCS, intensityTS)
% Block counter
iBlock = iBlock + 1;





%% STOP SEQUENCE: Run this to abort the sequence 
bd.stop




%% FUNCTIONS 






%% setupBlock
function [times, channel, markers] = setupBlock(nSinglePulses, nPairedPulses, ITI, JITTER, PAIRED_DELAY, portCS, portTS)
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

channelTest = portTS.*ones(nTrials, 1);
channelConditioning = portCS.*ones(nPairedPulses, 1);
channel = [channelTest; channelConditioning];


% Bring arrays into shared order:
[times, trialOrder] = sort([timesRight'; timesConditioning']);
channel = channel(trialOrder);
markers = markers(trialOrder);
end