%% Evaluate a HighDef conventional paired pulse session
addpath('B:\Projects\2020-12 VETTERPHD Project\libraries\neurone_tools_for_matlab_1.1.3.11_mod')
addpath(genpath('B:/Projects/2023-01 HighDef/libraries/vetter'))

x = get(0,'screensize');
SCREEN_WIDTH = x(3); SCREEN_HEIGHT=x(4);
clear('x');
VIS = true;

USmarker = 4;
CSmarker = 8;
TSmarker = 1;

resultName = struct();
resultName.R = 'R2L';
resultName.L = 'L2R';
coil_ids = struct();
coil_ids.R = 1;
coil_ids.L = 0;

%%
OUT = 'B:/Projects/2023-01 HighDef/Results';
ROOT = 'B:/Experimental Data/2023-01 HighDef/Conventional';
tableFile = sprintf('%s/Sessions.xlsx', ROOT);
T = readtable(tableFile);
subjects = unique(string(T.Subject));
fprintf('Listed subjects:\n%s\n\n', sprintf('   > %s\n', subjects))

ID = string(input('Choose subject > ', "s"));

if strcmpi(ID, "all")
    IDs = subjects;
else
    assert(ismember(ID, subjects), sprintf('Please choose only one of the listed subjects \n(update table at %s if needed)', tableFile))
    IDs = [ID];
end

for ID = IDs'

    IN = sprintf('%s/%s', ROOT, ID);
    ROWS = T(strcmpi(T.Subject, ID),:);
    ROWS = ROWS(startsWith(ROWS.Session, 'map') & ismember(ROWS.Condition, {'High', 'Low'}) & ~endsWith(ROWS.Session, 'old'),:);

    sessionTypes = unique(ROWS.Session)';

    %% Read and process all data

    for c_sessionType = sessionTypes
        sessionType = char(c_sessionType{:});
        suffix = "";
        if strlength(extractAfter(sessionType, '2')) > 1
            suffix = extractAfter(sessionType, '2');
            suffix = extractAfter(suffix, 1);
        end

        hemisphere = extractBefore(extractAfter(sessionType, '-'), '2');
        fprintf('\n\n Extracting all data of %s for map of %s hemisphere \n\n', ID, hemisphere)

        % Each events_... struct will have fields .US, .CS, .TS, labelling the
        % events as unconditioned stimulus, conditinING stimulus, test stimulus
        % respectively, and field .time
        events_planned  = initialize_events();
        intensities = [];
        events_localite = initialize_events();
        transforms = []; % From Localite recording
        coil_seen = [];
        coordinate_space = {};
        events_neurone  = initialize_events();
        meps_fdi = []; % From neurone recording
        meps_apb = [];
        meps_adm = [];

        PI_fdi = []; % From neurone recording
        PI_apb = [];
        PI_adm = [];

        ROWS_hemi = ROWS(endsWith(ROWS.Session, sessionType),:);

        for i = 1:size(ROWS_hemi, 1)
            row = ROWS_hemi(i,:);
            fprintf('  Block [%d]\n', row.Block);
            % Get planned events and add to structure:
            fprintf('   >> Planning\n');
            planning = read_plan(ROOT, ID, sessionType, row.Block);
            nTrialsPlanned = size(planning, 1);
            block = repmat(row.Block, 1, nTrialsPlanned);
            events_planned = append_events(events_planned, block, planning.Time', ...
                (planning.Marker == USmarker)', (planning.Marker == CSmarker)', (planning.Marker == TSmarker)');
            intensities = [intensities planning.Intensity_CS_MSO'];

            % Get neurOne events:
            fprintf('   >> NeurOne Recording\n');
            [neurone_times, US, CS, TS, mepsFDI_, mepsAPB_, mepsADM_, preInnervationFDI, preInnervationAPB, preInnervationADM] = ...
                read_neurone(ROOT, ID, sessionType, row.NeurOneFolder{:}, row.NeurOneIndex);
            events_neurone = append_events(events_neurone, block, neurone_times, US, CS, TS);
            meps_fdi = [meps_fdi mepsFDI_];
            meps_apb = [meps_apb mepsAPB_];
            meps_adm = [meps_adm mepsADM_];
            PI_fdi = [PI_fdi preInnervationFDI];
            PI_apb = [PI_apb preInnervationAPB];
            PI_adm = [PI_adm preInnervationADM];

            % Get localite events:
            fprintf('   >> Localite Recording\n');
            [localite_times, transforms_, coordinate_space_] = read_localite(sprintf('%s/%s/%s/%s', ROOT, ID, sessionType, row.LocaliteFolder{:}), coil_ids.(hemisphere), row.LocaliteTimeStamp{:});
            [US, CS, TS] = annotate_localite(localite_times);
            missing = coil_positions_missing(transforms_, 100);

            % Reject Localite events if no matching NeurOne event exists
            % (if orange button on coil was pressed accidentally)
            localite_times_wrt_1st_trial_in_s = (localite_times - min(localite_times)) / 1000;
            neurone_times_wrt_1st_trial_in_s  = (neurone_times - min(neurone_times));
            
            delay_to_nearest_match = min(abs(localite_times_wrt_1st_trial_in_s - neurone_times_wrt_1st_trial_in_s'));
            localite_event_found_in_neurone = delay_to_nearest_match < 0.05; % Only works if the additional events are sufficiently distant from a pulse
            if length(localite_times) ~= length(neurone_times) && any(~localite_event_found_in_neurone)
                error("\nSUPERFLUOUS localite events detected!\n Localite events %s at %s have no match\n\n", sprintf(" %d", find(~localite_event_found_in_neurone)), sprintf(" %d ms", localite_times(~localite_event_found_in_neurone)));
            end

            events_localite = append_events(events_localite, block, localite_times, US, CS, TS);
            coil_seen = [coil_seen ~missing];
            transforms = cat(3, transforms, transforms_);
            coordinate_space = [coordinate_space; repmat({coordinate_space_}, length(localite_times), 1)];
        end


        fprintf("FDI preinnervation:\n TS: %3.0f %% of trials showed pre-innervation > 50 µV\n US: %3.0f %% of trials > 50 µV\n", 100*mean(PI_fdi(events_neurone.TS) > 50), 100*mean(PI_fdi(events_neurone.US) > 50))
        fprintf("FDI preinnervation:\n TS: %3.0f %% of trials > 100 µV\n US: %3.0f %% of trials > 100 µV\n", 100*mean(PI_fdi(events_neurone.TS) > 100), 100*mean(PI_fdi(events_neurone.US) > 100))


        coil_seen = logical(coil_seen); % for some reason
        % Additionally, it can happen that only CS xor TS is seen
        % In that case, reject BOTH (alternative: transfer the transformation
        % matrix to the unseen)
        maskCS = coil_seen & events_neurone.CS;
        maskTS = coil_seen & events_neurone.TS;
        followingTSseen = maskCS & [maskTS(2:end) true];
        precedingCSseen = maskTS & [true maskCS(1:end-1)];
        coil_seen(maskCS) = followingTSseen(maskCS);
        coil_seen(maskTS) = precedingCSseen(maskTS);


        % QC: visually inspect coil positions:
        if VIS
            figure(Name=sessionType);
            x = squeeze(transforms(1,4,coil_seen));
            y = squeeze(transforms(2,4,coil_seen));
            z = squeeze(transforms(3,4,coil_seen));

            x_missing = squeeze(transforms(1,4,~coil_seen));
            y_missing = squeeze(transforms(2,4,~coil_seen));
            z_missing = squeeze(transforms(3,4,~coil_seen));
            %plot3(x, y, z, 'k.')

            o = zeros(sum(coil_seen), 1);
            plot3((x + [o squeeze(transforms(1,1,coil_seen))])', (y + [o squeeze(transforms(2,1,coil_seen))])', (z + [o squeeze(transforms(3,1,coil_seen))])', 'r')
            hold on
            plot3((x + [o squeeze(transforms(1,2,coil_seen))])', (y + [o squeeze(transforms(2,2,coil_seen))])', (z + [o squeeze(transforms(3,2,coil_seen))])', 'g')
            plot3((x + [o squeeze(transforms(1,3,coil_seen))])', (y + [o squeeze(transforms(2,3,coil_seen))])', (z + [o squeeze(transforms(3,3,coil_seen))])', 'b')
            plot3(x_missing,y_missing,z_missing,'kx')
            limit_range = max([range(xlim) range(ylim) range(zlim)]);
            window = 0.5.*[-limit_range limit_range];
            xlim(mean(xlim) + window); ylim(mean(ylim) + window); zlim(mean(zlim) + window)
            axis vis3d
            xlabel('x'); ylabel('y'); zlabel('z');


            fig_prein = figure(Name=sprintf("Preinnervation %s", sessionType));
            ax1 = subplot(1,3,[1 2]);
            plot(PI_fdi, 'r.')
            hold on
            plot(PI_adm, 'k.')
            plot(PI_apb, 'b.')
            ax2 = subplot(1,3,3);
            histogram(PI_fdi, FaceColor="r", Orientation="horizontal")
            hold on
            histogram(PI_adm, FaceColor="k", Orientation="horizontal")
            histogram(PI_apb, FaceColor="b", Orientation="horizontal")
            linkaxes([ax1 ax2], "y")
            exportgraphics(fig_prein, sprintf("%s/QC/qc_preinnervation_%s_%s_%s.pdf", OUT, ID, resultName.(hemisphere), suffix))
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

        % Check alignment of events:
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

        fprintf("\nOrganizing APB (%s, %s)\n", ID, sessionType)
        [CsE_APB_in_uV, SIHIscore_APB, inhibited_mep_APB_in_uV, preinnervation_APB_in_uV] = organizeMEPsForTable(meps_apb, PI_apb, events_neurone, coil_seen, VIS, 200);
        set(gcf, 'Name', sprintf('APB for %s', sessionType));

        fprintf("\nOrganizing FDI (%s, %s)\n", ID, sessionType)
        [CsE_FDI_in_uV, SIHIscore_FDI, inhibited_mep_FDI_in_uV, preinnervation_FDI_in_uV] = organizeMEPsForTable(meps_fdi, PI_fdi, events_neurone, coil_seen, VIS, 50);
        set(gcf, 'Name', sprintf('FDI for %s', sessionType));

        fprintf("\nOrganizing ADM (%s, %s)\n", ID, sessionType)
        [CsE_ADM_in_uV, SIHIscore_ADM, inhibited_mep_ADM_in_uV, preinnervation_ADM_in_uV] = organizeMEPsForTable(meps_adm, PI_adm, events_neurone, coil_seen, VIS, 200);
        set(gcf, 'Name', sprintf('ADM for %s', sessionType));
        drawnow;



        Block = events_planned.block(events_planned.CS & coil_seen)';
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


        Coordinate_space = coordinate_space(events_localite.CS & coil_seen);

        Result = table(Block, Intensity_percentMSO, x1, x2, x3, y1, y2, y3, z1, z2, z3, p1, p2, p3, Coordinate_space, ...
            CsE_APB_in_uV, CsE_FDI_in_uV, CsE_ADM_in_uV, ...
            SIHIscore_APB, SIHIscore_FDI, SIHIscore_ADM, ...
            inhibited_mep_APB_in_uV, inhibited_mep_FDI_in_uV, inhibited_mep_ADM_in_uV, ...
            preinnervation_APB_in_uV, preinnervation_FDI_in_uV, preinnervation_ADM_in_uV);

        writeWhere = sprintf('%s/%s_%s%s_raw.csv', OUT, ID, resultName.(hemisphere), suffix);


        %writeWhere = sprintf('%s/%s_%s_raw.csv', OUT, ID, sessionType);
        writetable(Result, writeWhere)

        fprintf('<< Wrote %d rows to %s\n%s\n', size(Result, 1), writeWhere, repmat('_', 1, 80))

    end

    close all
end










%% %% %% %% %% %% %% %% %% F U N C T I O N S %% %% %% %% %% %% %% %% %%

function [CsE_in_uV, SIHIscore, inhibited_mep_in_uV, preinnervation_paired_pulse] = organizeMEPsForTable(meps, preinnervation, events_neurone, coil_seen, VIS, pi_cutoff)

% Preinnervation is only used to reject trials from the computation of the
% SIHI-reference (mean unconditioned) at this point!

UR = meps(events_neurone.US);
nBlocks = length(unique(events_neurone.block));
UR_intercepts = nan(nBlocks, 1);
UR_slopes = nan(nBlocks, 1);

SIHIscore = [];


for b = unique(events_neurone.block)
    blockMask = events_neurone.block == b;
    blockMaskUS = events_neurone.block(events_neurone.US) == b;
    R = table(UR(blockMaskUS)', events_neurone.time(events_neurone.US & blockMask)', preinnervation(blockMaskUS)', ...
        VariableNames={'UR', 'Time', 'PI'});

    fprintf("[INFO] UR-rejection:\n")
    fprintf("starting from\t %d trials\n", size(R,1))
    R = R(R.UR > 50, :); % Reject all trials without MEP!
    fprintf("   MEP > 50µV\t %d trials\n", size(R,1))
    R = R(R.PI <= pi_cutoff, :); % Reject all trials with too much preinnervation!
    fprintf("    PI < %dµV\t %d trials\n", pi_cutoff, size(R,1))


    mdl = fitlm(R, 'UR~1+Time');
    mdl0 = fitlm(R, 'UR~1');
    LR = 2*(mdl.LogLikelihood - mdl0.LogLikelihood); % has a X2 distribution with a df equals to number of constrained parameters, here: 1
    pval = 1 - chi2cdf(LR, 1);

    if pval < 0.05
        UR_intercepts(b) = table2array(mdl.Coefficients(1,1));
        UR_slopes(b) = table2array(mdl.Coefficients(2,1));
    else
        UR_intercepts(b) = mean(R.UR);
        UR_slopes(b) = 0;
    end

    CR_in_this_block = meps(events_neurone.TS & coil_seen & blockMask);
    time_of_CR_in_this_block = events_neurone.time(events_neurone.TS & coil_seen & blockMask);
    mean_UR_at_that_time = UR_intercepts(b) + UR_slopes(b)*time_of_CR_in_this_block;
    mean_UR_at_that_time(mean_UR_at_that_time < 50) = 50; % minimum MEP response is hard-enforced to avoid numeric errors.
    SIHIscore = [SIHIscore; -log(CR_in_this_block ./ mean_UR_at_that_time)'];
end

CR = meps(events_neurone.TS & coil_seen);
XR = meps(events_neurone.CS & coil_seen);
CsE_in_uV = XR';
inhibited_mep_in_uV = CR';

preinnervation_paired_pulse = preinnervation(events_neurone.TS & coil_seen)';

if VIS
    figure;
    subplot(2,3,[1 2 3])
    nBlocks = length(unique(events_neurone.block));
    colors = [linspace(0, 1, nBlocks); linspace(0, 1, floor(nBlocks/2)) linspace(1, 0, ceil(nBlocks/2)); linspace(1,0,nBlocks)]';
    colors = colors(randperm(nBlocks),:);
    for b = unique(events_neurone.block)
        hold on;
        blockMaskUS = events_neurone.block == b;
        visualOffset = max(events_neurone.time)*b;
        x = events_neurone.time(blockMaskUS & events_neurone.US);
        y = meps(blockMaskUS & events_neurone.US);
        plot(x+visualOffset, y,'x', Color=colors(b,:))
        plot(mean(x)+visualOffset, mean(y), 'o', Color=colors(b,:), MarkerFaceColor=colors(b,:))
        plot([min(x) max(x)]+visualOffset, UR_intercepts(b) + UR_slopes(b).*[min(x) max(x)], '-', LineWidth=2, Color=colors(b,:))
        % TODO: plot block regression
    end
    ylabel('UR amplitude [\muV]'); xlabel('Progress through session');

    bins = linspace(0, max([max(XR) max(CR) max(UR)]), 50);

    subplot(2,3,4)
    histogram(XR, BinEdges=bins, FaceColor=[1 0 0], Normalization="probability")
    legend({'Exctiability responses'})
    subplot(2,3,5)
    histogram(CR, BinEdges=bins, FaceColor=[0 0 1], Normalization="probability")
    hold on
    histogram(UR, BinEdges=bins, FaceColor=[0 0 0], Normalization="probability")
    xline(mean(UR), 'k:')
    legend({'Conditioned', 'Unconditioned', 'mean unconditioned'})
    subplot(2,3,6)
    histogram(SIHIscore, FaceColor=[1 0 0], Normalization="probability")
    legend({'SIHIscore'})
end
end




function planning = read_plan(basepath, subject, session, block)
files = dir(sprintf('%s/%s/%s/%s_main_pulse_sequence_block-%d_*.csv', basepath, subject, session, subject, block));
if length(files) > 1
    warning('Multiple planning files found, using last: %s', files(end).name)
end
planning = readtable(sprintf('%s/%s', files(end).folder, files(end).name));

end




function [US, CS, TS] = annotate_localite(time)
LocalitePairedPulseTolerance = 200; % there is a fair bit of delay; used for annotation; 200 ms
timeToNextPulse = [diff(time) inf]; % add infinite time after last pulse
timeSinceLastPulse = [inf diff(time)]; % add infinite time before first pulse

US = (timeToNextPulse > 500) & (timeSinceLastPulse > 500); % Single pulses have substantial distance in both directions!
TS = (timeSinceLastPulse >= 0 & timeSinceLastPulse < LocalitePairedPulseTolerance & timeToNextPulse > 500);
CS  = (timeSinceLastPulse > 500 & timeToNextPulse >= 0 & timeToNextPulse < LocalitePairedPulseTolerance);
end


function [time, US, CS, TS, mepFDI, mepAPB, mepADM, ...
    preInnervationFDI, preInnervationAPB, preInnervationADM] = ...
    read_neurone(basepath, subject, session, folder, index)

% Parsing the session name into the muscles needed
hemisphereCode = extractAfter(session, '-');
hemisphereCode = hemisphereCode(1:3);
% If we are in session map-R2L, then the US and TS are delivered to the
% Left M1, and the corresponding MEP is in the RIGHT HAND!
us_handside = lower(extractBefore(hemisphereCode, '2'));
ts_handside = us_handside;
cs_handside = lower(extractAfter(hemisphereCode, '2'));

fprintf('    Processing %s as: CS @ %s   TS @ %s   US @ %s \t(hand!, opposite of coil)\n', session, cs_handside, ts_handside, us_handside)


subjData = module_read_neurone(sprintf("%s/%s/%s/%s", basepath, subject, session, folder), sessionPhaseNumber=index);
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

MEPwindow = window2index(0.02, 0.04, subjData.properties.samplingRate);
MEPwindow_cs = MEPwindow; % DO NOT! + subjData.properties.samplingRate*0.010; % Shift by 10 ms



mepFDI = nan(1, length(outMarkers));
mepFDI(US) = range(subjData.signal.(sprintf('FDI%s', us_handside)).data(triggerIdcs(US) + MEPwindow)');
mepFDI(CS) = range(subjData.signal.(sprintf('FDI%s', cs_handside)).data(triggerIdcs(CS) + MEPwindow_cs)');
mepFDI(TS) = range(subjData.signal.(sprintf('FDI%s', ts_handside)).data(triggerIdcs(TS) + MEPwindow)');

mepAPB = nan(1, length(outMarkers));
mepAPB(US) = range(subjData.signal.(sprintf('APB%s', us_handside)).data(triggerIdcs(US) + MEPwindow)');
mepAPB(CS) = range(subjData.signal.(sprintf('APB%s', cs_handside)).data(triggerIdcs(CS) + MEPwindow_cs)');
mepAPB(TS) = range(subjData.signal.(sprintf('APB%s', ts_handside)).data(triggerIdcs(TS) + MEPwindow)');

mepADM = nan(1, length(outMarkers));
mepADM(US) = range(subjData.signal.(sprintf('ADM%s', us_handside)).data(triggerIdcs(US) + MEPwindow)');
mepADM(CS) = range(subjData.signal.(sprintf('ADM%s', cs_handside)).data(triggerIdcs(CS) + MEPwindow_cs)');
mepADM(TS) = range(subjData.signal.(sprintf('ADM%s', ts_handside)).data(triggerIdcs(TS) + MEPwindow)');


preinnervation_window = window2index(-0.510, -0.010, subjData.properties.samplingRate);
preinnervation_window_ts = window2index(-0.520, -0.020, subjData.properties.samplingRate); % TS is preceded by CS 10ms beforehand -- thus shift by 10 ms

preInnervationFDI = nan(1, length(outMarkers));
preInnervationFDI(US) = range(detrend(subjData.signal.(sprintf('FDI%s', us_handside)).data(triggerIdcs(US) + preinnervation_window)', 1));
preInnervationFDI(CS) = range(detrend(subjData.signal.(sprintf('FDI%s', cs_handside)).data(triggerIdcs(CS) + preinnervation_window)', 1));
preInnervationFDI(TS) = range(detrend(subjData.signal.(sprintf('FDI%s', ts_handside)).data(triggerIdcs(TS) + preinnervation_window_ts)', 1));

preInnervationAPB = nan(1, length(outMarkers));
preInnervationAPB(US) = range(detrend(subjData.signal.(sprintf('APB%s', us_handside)).data(triggerIdcs(US) + preinnervation_window)', 1));
preInnervationAPB(CS) = range(detrend(subjData.signal.(sprintf('APB%s', cs_handside)).data(triggerIdcs(CS) + preinnervation_window)', 1));
preInnervationAPB(TS) = range(detrend(subjData.signal.(sprintf('APB%s', ts_handside)).data(triggerIdcs(TS) + preinnervation_window_ts)', 1));

preInnervationADM = nan(1, length(outMarkers));
preInnervationADM(US) = range(detrend(subjData.signal.(sprintf('ADM%s', us_handside)).data(triggerIdcs(US) + preinnervation_window)', 1));
preInnervationADM(CS) = range(detrend(subjData.signal.(sprintf('ADM%s', cs_handside)).data(triggerIdcs(CS) + preinnervation_window)', 1));
preInnervationADM(TS) = range(detrend(subjData.signal.(sprintf('ADM%s', ts_handside)).data(triggerIdcs(TS) + preinnervation_window_ts)', 1));
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

function indices = window2index(from, to, samplingrate)
indices = round(from*samplingrate):round(to*samplingrate);
end