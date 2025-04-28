rng('shuffle')
initialize_experiment
%%
nBlocks = 4;
nTrialsPerBlock = 100;

ITI = 4; %s
JITTER = 1; %s, plusminus, i.e. Unif([3,5])

expectedDuration = seconds(nTrialsPerBlock*ITI); expectedDuration.Format = 'mm:ss';
stdDuration = seconds(sqrt(((nTrialsPerBlock)*((2*JITTER)^2)/12))); stdDuration.Format = 's';
% Variance of uniform distr over [a,b]: ((b-a)^2)/12; here: b-a=2*JITTER;
% Variance of sum of iid RV = sum of variances of iid RV; sqrt(VAR) = sd;

% Print Experimental design settings
fprintf(' Total # pulses:    %d\n\n', nBlocks*nTrialsPerBlock)
fprintf(' Estimated time per block: %s min +/- %s\n', expectedDuration, stdDuration)
% Print total duration: expectedDuration.Format = 'hh:mm'
fprintf('         Total duration >= %s min\n\n', nBlocks*expectedDuration)


SUBJECT = arrayfun(@num2str, zeros(1,3));
name = flip(input('Enter subject ID > sub-', 's')); % for pilots: prefix with p
SUBJECT(1:length(name)) = name;
SUBJECT = flip(SUBJECT);
ID = sprintf('sub-%s', SUBJECT);

% Output location
OUT_ = 'Z:/Experimental Data/2023-01 HighDef/Conventional';
OUT = sprintf('%s/sub-%s/initiation', OUT_, SUBJECT);
mkdir(OUT)


% Display warning before starting Hotspot search
warning('Check  C O I L and T R A C K E R orientation! Orange button should be on the left!')
warning('Check Stimulator settings! W A V E F O R M  needs to be biphasic!')


%% Search Left hotspot
searchHotspot(stimulators.left, bd, 1)

%% Search Right hotspot
searchHotspot(stimulators.right, bd, 2)

%% Enter estimated  RMTs here:
% Count as MEP iff the amplitude is >= 50 uV
% RMT: 5/10 (i.e. ~half) of pulses yield an MEP over 50 uV
rmtLeft  = input('Enter RMT for  L E F T  coil > ');
rmtRight = input('Enter RMT for  R I G H T  coil > ');
fprintf('\n[RMT]\n      L = %d %% MSO\n      R = %d %% MSO\n\n', rmtLeft, rmtRight)


intensityLeft  = round(1.5 * rmtLeft);
intensityRight = round(1.5 * rmtRight);

% print stimulation intensities
fprintf('\n[SI]\n      L = %d %% MSO\n      R = %d %% MSO\n\n', intensityLeft, intensityRight)



%% (7.1) IO Curve LEFT coil
% IO Curve uses 5 different intensities with range RMT and 140% RMT
maxNumPulsesIO = 30; % nr. of pulses per intensity (default: 5)
[intensities, meps, fig, fitresult_left_coil] = runIOCurve(bd, stimulators.left, 1, rmtLeft, maxNumPulsesIO, 1, 2);

save([OUT filesep sprintf('%s_mep_io-left_coil_%s.mat', ID, datestr(now,'YYYY-mm-dd_HHMMSS'))], 'intensities', 'meps', 'fitresult_left_coil')% save to mat file
print(fig, sprintf('%s/%s_mep_io-left_coil_%s.pdf', OUT, ID, datestr(now,'YYYY-mm-dd_HHMMSS')), '-dpdf', '-r0') % print to pdf
fitresult_left_coil


%% (7.2) IO Curve RIGHT coil
% IO Curve uses 5 different intensities with range RMT and 140% RMT
[intensities, meps, fig, fitresult_right_coil] = runIOCurve(bd, stimulators.right, 2, rmtRight, maxNumPulsesIO, 3, 4);

save([OUT filesep sprintf('%s_mep_io-right_coil_%s.mat', ID, datestr(now,'YYYY-mm-dd_HHMMSS'))], 'intensities', 'meps', 'fitresult_right_coil') % save to mat file
print(fig, sprintf('%s/%s_mep_io-right_hand_%s.pdf', OUT, ID, datestr(now,'YYYY-mm-dd_HHMMSS')), '-dpdf', '-r0') % print to pdf
fitresult_right_coil


%%
intensitySIHI_TS_right = round(fitresult_right_coil.s)+5;
intensitySIHI_TS_left = round(fitresult_left_coil.s)+5;


%% SIHI-curve Left -> Right
% % SIHI curve using rmtLeft and intensityRight 
FDIl = 4;
maxNumPairedPulses = 40;
sc = bnp_sihi_curve(bd, rmtLeft,intensitySIHI_TS_right, stimulators.left, stimulators.right, FDIl, maxNumPairedPulses, 1, 2);

% save SIHI Curve
save(sprintf('%s/sub-%s_sihi_curve_L2R_SITS-%d-%s.mat', OUT, SUBJECT, intensitySIHI_TS_right, datestr(now,'YYYY-mm-dd_HHMMSS')), 'SUBJECT', 'sc')


%% SIHI-curve Right -> Left
FDIr = 2;
maxNumPairedPulses = 40;
sc = bnp_sihi_curve(bd, rmtRight, intensitySIHI_TS_left, stimulators.right, stimulators.left, FDIr, maxNumPairedPulses, 2, 1);

% save SIHI Curve
save(sprintf('%s/sub-%s_sihi_curve_R2L_SITS-%d-%s.mat', OUT, SUBJECT, intensitySIHI_TS_left, datestr(now,'YYYY-mm-dd_HHMMSS')), 'SUBJECT', 'sc')


%% Blocks
%T = readtable('//storage.neurologie.uni-tuebingen.de/bbnp_lab/Experimental Data/2023-01 HighDef/Conventional/Planning.xlsx');
T = readtable('C:/Users/BNP Lab/Documents/HighDef-temp/Planning.csv');
subjectRow = T(strcmpi(T.Subject, ID), :);

assert(size(subjectRow, 1) == 1)

if strcmpi(subjectRow.x1st_initiation{:}, 'R')
    fprintf('Beginning with mapping of Right hemisphere!\n')
    channel_1st = 2;
    channel_2nd = 1;
    intensity_1st = intensityRight;
    intensity_2nd = intensityLeft;
    stimulator_1st = stimulators.right;
    stimulator_2nd = stimulators.left;
else
    fprintf('Beginning with mapping of Left hemisphere!\n')
    channel_1st = 1;
    channel_2nd = 2;
    intensity_1st = intensityLeft;
    intensity_2nd = intensityRight;
    stimulator_1st = stimulators.left;
    stimulator_2nd = stimulators.right;
end

intensity_1st = min([100 intensity_1st]);
intensity_2nd = min([100 intensity_2nd]);
iBlock1st = 1;
iBlock2nd = 1;
condition_1st = subjectRow.x1st_initiation{:};
condition_2nd = subjectRow.x2nd_initiation{:};


%% Mapping first hemisphere:
while iBlock1st <= nBlocks
    input(sprintf('Start block %d? [press any key]', iBlock1st), 's')
    
    [times, channels, markers] = setupBlock(nTrialsPerBlock, ITI, JITTER, channel_1st);
    filename = sprintf('%s/sub-%s_task-%s_run-%d_%s_planning', OUT, SUBJECT, subjectRow.x1st_initiation{:}, iBlock1st, datetime('now','Format','yyyy-mm-dd-HHMMSS'));
    save([filename '.mat'], "times", "channels", "markers", "iBlock1st", "condition_1st", "intensity_1st")
    planning = table(times, channels, markers, repmat({condition_1st}, length(times), 1), repmat(intensity_1st, length(times), 1), 'VariableNames', {'Time', 'Stimulator', 'Marker', 'Condition', 'Intensity_MSO'});
    writetable(planning, [filename '.csv'])

    fprintf('Running %s Block %d of %d   (~ %d s)\n', condition_1st, iBlock1st, nBlocks, times(end))

    pause(10) % get ready to map with conditioning pulse!

    % Run block
    runBlock(bd, times, channels, markers, stimulator_1st, stimulator_2nd, intensity_1st, intensity_2nd)

    
    iBlock1st = iBlock1st + 1;    
end


%% Mapping second hemisphere:
while iBlock2nd <= nBlocks
    input(sprintf('Start block %d? [press any key]', iBlock2nd), 's')
    
    [times, channels, markers] = setupBlock(nTrialsPerBlock, ITI, JITTER, channel_2nd);
    filename = sprintf('%s/sub-%s_task-%s_run-%d_%s_planning', OUT, SUBJECT, subjectRow.x2nd_initiation{:}, iBlock2nd, datetime('now','Format','yyyy-mm-dd-HHMMSS'));
    save([filename '.mat'], "times", "channels", "markers", "iBlock2nd", "condition_2nd", "intensity_2nd")
    planning = table(times, channels, markers, repmat({condition_2nd}, length(times), 1), repmat(intensity_2nd, length(times), 1), 'VariableNames', {'Time', 'Stimulator', 'Marker', 'Condition', 'Intensity_MSO'});
    writetable(planning, [filename '.csv'])

    fprintf('Running %s Block %d of %d   (~ %d s)\n', condition_2nd, iBlock2nd, nBlocks, times(end))

    pause(10) % get ready to map with conditioning pulse!

    % Run block
    runBlock(bd, times, channels, markers, stimulator_1st, stimulator_2nd, intensity_1st, intensity_2nd)

    
    iBlock2nd = iBlock2nd + 1;    
end








function [times, channels, markers] = setupBlock(nTrials, ITI, JITTER, channel)
minITI = ITI - JITTER;
ITIrange = 2*JITTER;
ITIs = minITI + (rand(1,nTrials) .* ITIrange);
times = cumsum(ITIs)';
channels = channel.*ones(nTrials, 1);
markers = channel.*ones(nTrials, 1);
end
