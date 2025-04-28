%% Evaluate a HighDef conventional session
addpath('B:\Projects\2020-12 VETTERPHD Project\libraries\neurone_tools_for_matlab_1.1.3.11_mod')
addpath(genpath('B:/Projects/2023-01 HighDef/libraries/vetter'))
%addpath(genpath('B:/Projects/2023-01 HighDef/libraries/visualization'))
%addpath(genpath('B:/Projects/2020-12 VETTERPHD Project/004-MoCsEFC/MoCsEFC/src'))

x = get(0,'screensize');
SCREEN_WIDTH = x(3); SCREEN_HEIGHT=x(4);
clear('x');
VIS = false;

USmarker = 4;
CSmarker = 8;
TSmarker = 1;

%%
OUT = 'B:/Projects/2023-01 HighDef/Results';
ROOT = 'B:/Experimental Data/2023-01 HighDef/Conventional';
tableFile = sprintf('%s/Sessions.xlsx', ROOT);
T = readtable(tableFile);
subjects = unique(T.Subject);
fprintf('Listed subjects: %s\n\n', sprintf(' %s', subjects{:}))

ID = input('Choose subject > ', 's');
assert(ismember(ID, subjects), sprintf('Please choose only one of the listed subjects \n(update table at %s if needed)', tableFile))

IN = sprintf('%s/%s', ROOT, ID);
ROWS = T(strcmpi(T.Subject, ID),:);
ROWS = ROWS(strcmpi(ROWS.Condition, 'Map'),:);

%% Read and process all data

% Each events_... struct will have fields .US, .CS, .TS, labelling the
% events as unconditioned stimulus, conditinING stimulus, test stimulus
% respectively, and field .time
events_planned  = initialize_events();
intensities = [];
events_localite = initialize_events();
transforms = []; % From Localite recording
coil_seen = [];
events_neurone  = initialize_events();
meps = []; % From neurone recording
meps_apb = [];



for i = 1:size(ROWS, 1)
    row = ROWS(i,:);
    fprintf('  Block [%d]\n', row.Block);
    % Get planned events and add to structure:
    fprintf('   >> Planning\n');
    planning = read_plan(ROOT, ID, row.Block);
    nTrialsPlanned = size(planning, 1);
    block = repmat(row.Block, 1, nTrialsPlanned);
    events_planned = append_events(events_planned, block, planning.Time', ...
        (planning.Marker == USmarker)', (planning.Marker == CSmarker)', (planning.Marker == TSmarker)');
    intensities = [intensities planning.Intensity_LCoil_MSO'];

    % Get localite events:
    fprintf('   >> Localite Recording\n');
    [localite_times, transforms_] = read_localite(sprintf('%s/%s/%s', ROOT, ID, row.LocaliteFolder{:}), 0, row.LocaliteTimeStamp{:});
    [US, CS, TS] = annotate_localite(localite_times);
    missing = coil_positions_missing(transforms_);
    events_localite = append_events(events_localite, block, localite_times, US, CS, TS);
    coil_seen = [coil_seen ~missing];
    transforms = cat(3, transforms, transforms_);

    % Get neurOne events:
    fprintf('   >> NeurOne Recording\n');
    [neurone_times, US, CS, TS, meps_, mepsAPB_] = read_neurone(ROOT, ID, row.NeurOneFolder{:}, row.NeurOneIndex);
    events_neurone = append_events(events_neurone, block, neurone_times, US, CS, TS);
    meps = [meps meps_];
    meps_apb = [meps_apb mepsAPB_];
end

coil_seen = logical(coil_seen); % for some reason

%% QC: visually inspect coil positions:
if VIS
    x = squeeze(transforms(1,4,coil_seen));
    y = squeeze(transforms(2,4,coil_seen));
    z = squeeze(transforms(3,4,coil_seen));
    %plot3(x, y, z, 'k.')
    
    o = zeros(sum(coil_seen), 1);
    plot3((x + [o squeeze(transforms(1,1,coil_seen))])', (y + [o squeeze(transforms(2,1,coil_seen))])', (z + [o squeeze(transforms(3,1,coil_seen))])', 'r')
    hold on
    plot3((x + [o squeeze(transforms(1,2,coil_seen))])', (y + [o squeeze(transforms(2,2,coil_seen))])', (z + [o squeeze(transforms(3,2,coil_seen))])', 'g')
    plot3((x + [o squeeze(transforms(1,3,coil_seen))])', (y + [o squeeze(transforms(2,3,coil_seen))])', (z + [o squeeze(transforms(3,3,coil_seen))])', 'b')
    xlim([-20 85]); ylim([-20 85]); zlim([-20 85])
    
    axis vis3d
    xlabel('x'); ylabel('y'); zlabel('z');
end


% Approach for aligning:
% Goal: In each row of the table, want:
% - Intensity in % MSO (to be translated into dI/dt in [A/µs])
% - Coil transformation (16-4=12 entries)
% - Responses: MEP-amplitude (XR), SIHI-score, conditioned MEP-amplitude
%
% For aligning the neurOne trials and the localite trials, we can use the
% trial matching best in time, and reject trial for which the closest
% localite time-stamp is too far off the NeurOne marker time.
% If the coil position was e.g. not picked up during the conditioning
% pulse, but during the test-pulse (delivered by the other coil), this
% position should still be fine (10 ms later).
%
% For now, since the data seem to align well, just use direct
% cross-indexing. This may need to be robustified later.

%% Check alignment of events:
fprintf('[!] Checking correspondence between events:\n')
assert(all(events_planned.US == events_localite.US))
assert(all(events_planned.CS == events_localite.CS))    % A = B
assert(all(events_planned.TS == events_localite.TS))

assert(all(events_neurone.US == events_localite.US))
assert(all(events_neurone.CS == events_localite.CS))    % B = C
assert(all(events_neurone.TS == events_localite.TS))
% If A = B and B = C, then A = C (don't need to check that)

% Check timing too:
for block = unique(events_planned.block)
    tPlanned  = events_planned.time( events_planned.block  == block) .* 1e3;
    tLocalite = events_localite.time(events_localite.block == block);
    tNeurOne  = events_neurone.time( events_neurone.block  == block) .* 1e3;
    startTimePlanned  = min(tPlanned);
    startTimeLocalite = min(tLocalite);
    startTimeNeurOne  = min(tNeurOne);
    AE1 = abs((tPlanned - startTimePlanned) - (tLocalite - startTimeLocalite));
    AE2 = abs((tPlanned - startTimePlanned) - (tNeurOne - startTimeNeurOne));
    AE3 = abs((tLocalite - startTimeLocalite) - (tNeurOne - startTimeNeurOne));
    fprintf('     block %2.d: Planned  vs. Localite <= %.1f ms / mean = %4.1f ms\n', block, max(AE1), mean(AE1))
    fprintf('               Planned  vs. NeurOne  <= %.1f ms / mean = %4.1f ms\n', max(AE2), mean(AE2))
    fprintf('               Localite vs. NeurOne  <= %.1f ms / mean = %4.1f ms\n', max(AE3), mean(AE3))
end


fprintf(' ^^ correspondence verified\n\n')


%%
% Extract Unconditioned Responses, ConditionED responses, and eXcitability
% responses as follows:
UR = meps(events_neurone.US);
CR = meps(events_neurone.TS & coil_seen);
XR = meps(events_neurone.CS & coil_seen);


Intensity_percentMSO = intensities(events_planned.CS & coil_seen)';
x1 = squeeze(transforms(1,1,events_localite.CS & coil_seen));
x2 = squeeze(transforms(2,1,events_localite.CS & coil_seen));
x3 = squeeze(transforms(3,1,events_localite.CS & coil_seen));
y1 = squeeze(transforms(1,2,events_localite.CS & coil_seen));
y2 = squeeze(transforms(2,2,events_localite.CS & coil_seen));
y3 = squeeze(transforms(3,2,events_localite.CS & coil_seen));
z1 = squeeze(transforms(1,3,events_localite.CS & coil_seen));
z2 = squeeze(transforms(2,3,events_localite.CS & coil_seen));
z3 = squeeze(transforms(3,3,events_localite.CS & coil_seen));
p1 = squeeze(transforms(1,4,events_localite.CS & coil_seen));
p2 = squeeze(transforms(2,4,events_localite.CS & coil_seen));
p3 = squeeze(transforms(3,4,events_localite.CS & coil_seen));

CsE_in_uV = XR';
SIHIscore = -log(CR ./ mean(UR))';
inhibited_mep_in_uV = CR';

Result = table(Intensity_percentMSO, x1, x2, x3, y1, y2, y3, z1, z2, z3, p1, p2, p3, CsE_in_uV, SIHIscore, inhibited_mep_in_uV);
writetable(Result, sprintf('%s/%s_raw.csv', OUT, ID))















%% %% %% %% %% %% %% %% %% F U N C T I O N S %% %% %% %% %% %% %% %% %%



function planning = read_plan(basepath, subject, block)
files = dir(sprintf('%s/%s/%s_main_pulse_sequence_block-%d_*.csv', basepath, subject, subject, block));
if length(files) > 1
    warning('Multiple planning files found, using last: %s', files(end).name)
end
planning = readtable(sprintf('%s/%s', files(end).folder, files(end).name));

end




function [US, CS, TS] = annotate_localite(time)
LocalitePairedPulseTolerance = 33; % there is a fair bit of delay; used for annotation
timeToNextPulse = [diff(time) inf]; % add infinite time after last pulse
timeSinceLastPulse = [inf diff(time)]; % add infinite time before first pulse

US = (timeToNextPulse > 500) & (timeSinceLastPulse > 500); % Single pulses have substantial distance in both directions!
TS = (timeSinceLastPulse >= 0 & timeSinceLastPulse < LocalitePairedPulseTolerance & timeToNextPulse > 500);
CS  = (timeSinceLastPulse > 500 & timeToNextPulse >= 0 & timeToNextPulse < LocalitePairedPulseTolerance);
end


function [time, US, CS, TS, mepFDI, mepAPB] = read_neurone(basepath, subject, folder, index)
subjData = module_read_neurone(sprintf("%s/%s/%s", basepath, subject, folder), sessionPhaseNumber=index);
markerLabels = subjData.markers.type;

if ismember('144', unique(markerLabels))
    markerLabels = markerLabels(~strcmpi(markerLabels, '144'));
    markerLabels(~isnan(str2double(markerLabels))) = num2cell(cellfun(@(s) num2str(str2double(s) - 144), markerLabels(~isnan(str2double(markerLabels)))));
end

codeSingleRight = '4';
codePairedRight = '1';
codePairedLeft  = '8';

outMarkers = find(strcmpi(markerLabels, 'Out'));
truncatedOutMarkerNext = arrayfun(@(x) min([x length(markerLabels)]), outMarkers+1);

CS = (strcmpi(markerLabels(outMarkers-1), codePairedLeft) | strcmpi(markerLabels(truncatedOutMarkerNext), codePairedLeft))';
TS = (strcmpi(markerLabels(outMarkers-1), codePairedRight) | strcmpi(markerLabels(truncatedOutMarkerNext), codePairedRight))';
US = (strcmpi(markerLabels(outMarkers-1), codeSingleRight) | strcmpi(markerLabels(truncatedOutMarkerNext), codeSingleRight))';
triggerIdcs = subjData.markers.index(strcmpi(subjData.markers.type, 'Out')); % This is NOT the same as outMarkers! markerLabels that leads to outMarkers is edited.
time = subjData.markers.time(strcmpi(subjData.markers.type, 'Out'))';

MEPwindow_in_s = [0.021, 0.036];
MEPwindow = MEPwindow_in_s .* subjData.properties.samplingRate;
MEPwindow = round(MEPwindow(1)):round(MEPwindow(2));
MEPwindow_cs = MEPwindow + subjData.properties.samplingRate*0.010; % Shift by 10 ms

mepFDI = nan(1, length(outMarkers));
mepFDI(US) = range(subjData.signal.FDIl.data(triggerIdcs(US) + MEPwindow)');
mepFDI(CS) = range(subjData.signal.FDIr.data(triggerIdcs(CS) + MEPwindow_cs)');
mepFDI(TS) = range(subjData.signal.FDIl.data(triggerIdcs(TS) + MEPwindow)');

mepAPB = nan(1, length(outMarkers));
mepAPB(US) = range(subjData.signal.APBl.data(triggerIdcs(US) + MEPwindow)');
mepAPB(CS) = range(subjData.signal.APBr.data(triggerIdcs(CS) + MEPwindow_cs)');
mepAPB(TS) = range(subjData.signal.APBl.data(triggerIdcs(TS) + MEPwindow)');
end



function events = initialize_events()
events = struct();
events.block = [];
events.time = [];
events.US = false(0);
events.CS = false(0);
events.TS = false(0);
end

function events = append_events(events, block, time, US, CS, TS)
events = add2field(events, 'block', block);
events = add2field(events, 'time', time);
events = add2field(events, 'US', logical(US));
events = add2field(events, 'CS', logical(CS));
events = add2field(events, 'TS', logical(TS));
end

function target = add2field(target, field, newData, dim)
arguments
    target
    field
    newData
    dim = -1
end
if dim == -1
    dim = length(size(target.(field)));
end
target.(field) = cat(dim, target.(field), newData);
end