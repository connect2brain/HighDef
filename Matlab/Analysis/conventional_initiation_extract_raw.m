%% Evaluate a HighDef conventional initiation session
addpath('B:\Projects\2020-12 VETTERPHD Project\libraries\neurone_tools_for_matlab_1.1.3.11_mod')
addpath(genpath('B:/Projects/2023-01 HighDef/libraries/vetter'))

x = get(0,'screensize');
SCREEN_WIDTH = x(3); SCREEN_HEIGHT=x(4);
clear('x');
VIS = true;

coil_ids = struct();
coil_ids.R = 1;
coil_ids.L = 0;

contralateral_hand = struct();
contralateral_hand.R = 'l';
contralateral_hand.L = 'r';

%%
OUT = 'B:/Projects/2023-01 HighDef/Results';
ROOT = 'B:/Experimental Data/2023-01 HighDef/Conventional';
tableFile = sprintf('%s/Sessions.xlsx', ROOT);
T = readtable(tableFile);
subjects = unique(T.Subject);
fprintf('Listed subjects: %s\n\n', sprintf(' %s', subjects{:}))

ID = input('Choose subject > ', 's');
assert(ismember(ID, subjects), sprintf('Please choose only one of the listed subjects \n(update table at %s if needed)', tableFile))

IN = sprintf('%s/%s/initiation', ROOT, ID);
ALL_ROWS = T(strcmpi(T.Subject, ID),:);
ALL_ROWS = ALL_ROWS(strcmpi(ALL_ROWS.Session, 'initiation'),:);

%%
% I need to generate two csv files, one for left hemisphere mapping, one
% for right hemisphere mapping
% Each row: Intensity_percent_MSO, ...coil transform..., MEP-amplitude

for hemisphere = {'L', 'R'}

    fprintf('\n\nHemisphere: %s\n', hemisphere{:})
    planned_time = [];
    intensities = [];
    localite_time = [];
    coordinate_space = {};
    transforms = [];
    coil_seen = [];
    neurone_time = [];
    meps = [];

    PI = [];

    rows = ALL_ROWS(strcmpi(ALL_ROWS.Condition, sprintf('Map-%s', hemisphere{:})),:);
    for i = 1:size(rows, 1)
        row = rows(i,:);
        fprintf('  Block [%d]\n', row.Block);
        % Get planned events and add to structure:
        fprintf('   >> Planning\n');
        planning = read_plan(ROOT, ID, hemisphere{:}, row.Block);
        planned_time_ = planning.Time .* 1e3; % in ms

        % Get localite events:
        fprintf('   >> Localite Recording\n');

        [localite_time_, transforms_, coordinate_space_] = read_localite(sprintf('%s/%s', IN, row.LocaliteFolder{:}), coil_ids.(hemisphere{:}), row.LocaliteTimeStamp{:});
        localite_time_ = localite_time_'; % natively in ms
        missing = coil_positions_missing(transforms_);

        % Get neurOne events:
        fprintf('   >> NeurOne Recording\n');
        [neurone_time_, meps_struct, preinnervation_struct] = read_neurone(IN, row.NeurOneFolder{:}, row.NeurOneIndex);
        neurone_time_ = neurone_time_ .* 1e3; % in ms

        % Time consistency check
        startTimePlanned  = min(planned_time_);
        startTimeLocalite = min(localite_time_);
        startTimeNeurOne  = min(neurone_time_);
        AE1 = abs((planned_time_ - startTimePlanned) - (localite_time_ - startTimeLocalite));
        AE2 = abs((planned_time_ - startTimePlanned) - (neurone_time_ - startTimeNeurOne));
        AE3 = abs((localite_time_ - startTimeLocalite) - (neurone_time_ - startTimeNeurOne));
        fprintf('     block %2.d: Planned  vs. Localite <= %.1f ms / mean = %4.1f ms\n', i, max(AE1), mean(AE1))
        fprintf('               Planned  vs. NeurOne  <= %.1f ms / mean = %4.1f ms\n', max(AE2), mean(AE2))
        fprintf('               Localite vs. NeurOne  <= %.1f ms / mean = %4.1f ms\n', max(AE3), mean(AE3))


        % Append to collected data:
        planned_time = [planned_time; planned_time_];
        intensities = [intensities; planning.Intensity_MSO];
        localite_time = [localite_time; localite_time_];
        coordinate_space = [coordinate_space; repmat({coordinate_space_}, length(localite_time_), 1)];
        transforms = cat(3, transforms, transforms_);
        coil_seen = [coil_seen; ~missing'];
        neurone_time = [neurone_time; neurone_time_];

        hand = contralateral_hand.(hemisphere{:});
        meps_ = [meps_struct.(sprintf('APB%s', hand)), meps_struct.(sprintf('FDI%s', hand)), meps_struct.(sprintf('ADM%s', hand))];
        meps = [meps; meps_];
        PI_ = [preinnervation_struct.(sprintf("APB%s", hand)), preinnervation_struct.(sprintf("FDI%s", hand)), preinnervation_struct.(sprintf("ADM%s", hand))];
        PI = [PI; PI_];
    end

    coil_seen = logical(coil_seen); % for some reason

    % Write data:
    Intensity_percentMSO = intensities(coil_seen);
    x1 = squeeze(transforms(1,1,coil_seen));
    x2 = squeeze(transforms(2,1,coil_seen));
    x3 = squeeze(transforms(3,1,coil_seen));
    y1 = squeeze(transforms(1,2,coil_seen));
    y2 = squeeze(transforms(2,2,coil_seen));
    y3 = squeeze(transforms(3,2,coil_seen));
    z1 = squeeze(transforms(1,3,coil_seen));
    z2 = squeeze(transforms(2,3,coil_seen));
    z3 = squeeze(transforms(3,3,coil_seen));
    p1 = squeeze(transforms(1,4,coil_seen));
    p2 = squeeze(transforms(2,4,coil_seen));
    p3 = squeeze(transforms(3,4,coil_seen));
    Coordinate_space = coordinate_space(coil_seen);

    if VIS
        figure(Name=hemisphere{:});
        o = zeros(sum(coil_seen), 1);
        plot3((p1 + [o x1])', (p2 + [o x2])', (p3 + [o x3])', 'r')
        hold on
        plot3((p1 + [o y1])', (p2 + [o y2])', (p3 + [o y3])', 'g')
        plot3((p1 + [o z1])', (p2 + [o z2])', (p3 + [o z3])', 'b')
        xlim([-100 100]); ylim([-100 100]); zlim([-100 100])

        axis vis3d
        xlabel('x'); ylabel('y'); zlabel('z');
        title(hemisphere{:})
    end


    CsE_APB_in_uV = meps(coil_seen,1);
    CsE_FDI_in_uV = meps(coil_seen,2);
    CsE_ADM_in_uV = meps(coil_seen,3);

    preinnervation_APB_in_uV = PI(coil_seen, 1);
    preinnervation_FDI_in_uV = PI(coil_seen, 2);
    preinnervation_ADM_in_uV = PI(coil_seen, 3);

    Result = table(Intensity_percentMSO, x1, x2, x3, y1, y2, y3, z1, z2, z3, p1, p2, p3, Coordinate_space, CsE_APB_in_uV, CsE_FDI_in_uV, CsE_ADM_in_uV, preinnervation_APB_in_uV, preinnervation_FDI_in_uV, preinnervation_ADM_in_uV);
    writetable(Result, sprintf('%s/%s_%s_raw.csv', OUT, ID, hemisphere{:}))
    fprintf("Number of trials for %s %s: %d\n", ID, hemisphere{:}, size(Result, 1))
end








function planning = read_plan(basepath, subject, hemisphere, block)
look_for = sprintf('%s/%s/initiation/%s_task-%s_run-%d_*.csv', basepath, subject, subject, hemisphere, block);
files = dir(look_for);
if length(files) > 1
    warning('Multiple planning files found, using last: %s', files(end).name)
end
if length(files) < 1
    warning('No planning file found while looking for: %s', look_for)
end
planning = readtable(sprintf('%s/%s', files(end).folder, files(end).name));
end

function events = initialize_events()
events = struct();
events.block = [];
events.time = [];
end


function [time, mep_struct, preinnervation_struct] = read_neurone(basepath, folder, index)
subjData = module_read_neurone(sprintf("%s/%s", basepath, folder), sessionPhaseNumber=index);

triggerIdcs = subjData.markers.index(strcmpi(subjData.markers.type, 'Out')); % This is NOT the same as outMarkers! markerLabels that leads to outMarkers is edited.
time = subjData.markers.time(strcmpi(subjData.markers.type, 'Out'));

MEPwindow = window2index(0.021, 0.036, subjData.properties.samplingRate);
preinnervation_window = window2index(-0.510, -0.010, subjData.properties.samplingRate);

signalNames = fieldnames(subjData.signal);
signalNames = signalNames(endsWith(signalNames, ["l", "r"]))';

mep_struct = struct();
preinnervation_struct = struct();
for m = signalNames
    muscle = m{:};
    mep_struct.(muscle) = range(subjData.signal.(muscle).data(triggerIdcs + MEPwindow)')';
    preinnervation_struct.(muscle) = range(detrend(subjData.signal.(muscle).data(triggerIdcs + preinnervation_window)', 1))';
end
end

function indices = window2index(from, to, samplingrate)
indices = round(from*samplingrate):round(to*samplingrate);
end

