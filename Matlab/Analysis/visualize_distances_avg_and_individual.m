%% This script produces FIGURE 3 of the manuscript
% FS-avgerage brain is plotted, and the spots marked thereon in red and
% blue
% Below, the boxplots of the distances under the different methods,
% spot-comparisons etc shown.

% This script uses the results of the "conventional evaluation" --- make sure to run
% compare_SSAM_and_projection.m before

STORAGE = "D:/HighDef-operate/HighDef";

SPLIT = "all"; % "all", "first_half", "second_half"
suffix = sprintf("_%s", SPLIT);
if strcmpi(SPLIT, "all")
    suffix = "";
end

leftHanded = ["sub-006", "sub-016", "sub-017"];

muscles = ["ADM", "FDI", "APB"];
colors = [];
colors.ADM = [0 0.63 0];
colors.FDI = [0.63 0 0];
colors.APB = [0 0 0.63];
markers =[];
markers.ADM = '+';
markers.APB = 'x';
markers.FDI = 'o';
for muscle = muscles
    for other = muscles
        condition = sprintf('%svs%s', muscle, other);
        colors.(condition) = colors.(muscle) + colors.(other);
        markers.(condition) = '^';
    end
end


%% 1 Collect all spot data available for all subjects as a joint table:
%    Subject Session Hemisphere Spot Muscle x_source y_source z_source x_coil y_coil z_coil angle_coil R2
sessionNames = []; sessionNames.L = ["map-L", "map-L2R"]; sessionNames.R = ["map-R", "map-R2L"];
templates = []; templates.hot = "CsE_%s_in_uV"; templates.cold = "SIHIscore_%s";
basepath = STORAGE;
%subjects = arrayfun(@(s) string(s.name), dir(fullfile(basepath, 'sub-*')))';
subjects = ["sub-001", "sub-002", "sub-003", "sub-004", "sub-006", "sub-007", "sub-008", "sub-009", "sub-011", "sub-012", "sub-013", "sub-014", "sub-015", "sub-016", "sub-017", "sub-019", "sub-022", "sub-023"];

% Where are the hotspot data?
% - Source space:
% basepath/sub-???/results/exp_{exp_id}/r2/mesh_mesh0/roi_midlayer_{hemisphere}/CsE_{muscle}_in_uV/sigmoid4/r2_roi_data.hdf5
% - Coil space:
% basepath/sub-???/opt/{exp_id}/CsE_{muscle}_in_uV_[pos]/mesh0/E_mag/opt_coil_pos.hdf5

Result = [];
for ID = subjects
    for hemisphere = ["R", "L"]
        for muscle = ["ADM", "FDI", "APB"]
            roi_id = sprintf("midlayer_%s", lower(hemisphere));
            for session = 1:2
                exp_id = sessionNames.(hemisphere)(session);
                if session == 2 && muscle == "FDI"
                    spots = ["hot", "cold"];
                else
                    spots = ["hot"];
                end

                for spot = spots
                    response_id = sprintf(templates.(spot), muscle);
                    if strcmpi(SPLIT, "all")
                        newrow = collectSpot(basepath, ID, exp_id, roi_id, response_id);
                    else
                        newrow = collectSpot(basepath, ID, exp_id, roi_id, response_id, SPLIT);
                    end
                    if ~isempty(newrow)
                        newrow.Session = sprintf("S%d", session); newrow.Hemisphere = hemisphere; newrow.Spot = spot; newrow.Muscle = muscle;
                        fprintf(' %s/%s  :  S%d=%s %s-%s\n', ID, hemisphere, session, exp_id, muscle, spot)
                    else
                        fprintf(' missing %s/%s  :  S%d=%s %s-%s\n', ID, hemisphere, session, exp_id, muscle, spot)
                    end
                    Result = [Result; newrow];
                end
            end
        end
    end
end



% Add new hemisphere labels: dominant, nondominant
Result.Dominant = strcmpi(Result.Hemisphere, "L");

% For the single right-handed subject, flip this
% and also mirror the x-coordinate
mirror_which = ismember(Result.Subject, leftHanded);
Result.Dominant(mirror_which) = ~Result.Dominant(mirror_which);
Result.X_source(mirror_which) = -Result.X_source(mirror_which);



%% Plot distance histograms for:
% 1: FDI hot vs FDI cold
% 2: ADM hot vs APB hot
% 3: ADM hot vs FDI hot
% 4: FDI hot vs APB hot
% each: two histograms, one for left, one for right.

hemisphere_colors = []; 
%hemisphere_colors.L = "b"; hemisphere_colors.R = "r";
hemisphere_colors.L = "k"; 
hemisphere_colors.R = "k";
hemisphere_colors.dominant = "k";
hemisphere_colors.nondominant = "k";




fprintf("\n\nTest difference of spots on individual brain for SSAM (%s)\n", SPLIT)


distances = [];
coldspot_R2 = [];
for dominant = [true, false]
    if dominant
        dominance = "dominant";
        anterolateral_axis = [-1 1 -1] ./ sqrt(3);
        compare_along = [-1 -1 -1] ./ sqrt(3); % downward posterolateral (more along the surface)
    else
        dominance = "nondominant";
        anterolateral_axis = [1 1 -1] ./ sqrt(3);
        compare_along = [1 -1 -1] ./ sqrt(3); % downward posterolateral (more along the surface)
    end

    indices_A = find(Result.Spot == "hot" & Result.Muscle == "ADM" & Result.Session == "S2" & Result.Dominant == dominant);
    indices_B = arrayfun(@(i) find(Result.Spot == "hot" & Result.Subject == Result.Subject(i) & Result.Dominant == Result.Dominant(i) & Result.Muscle == "FDI" & Result.Session == "S2"), indices_A);
    differences = [Result.X_source(indices_B) - Result.X_source(indices_A), Result.Y_source(indices_B) - Result.Y_source(indices_A), Result.Z_source(indices_B) - Result.Z_source(indices_A)];
    scoring = differences * anterolateral_axis';
    distances.(dominance).ADMvsFDI = scoring; % sqrt((Result.X_source(indices_A) - Result.X_source(indices_B)).^2 + (Result.Y_source(indices_A) - Result.Y_source(indices_B)).^2 + (Result.Z_source(indices_A) - Result.Z_source(indices_B)).^2);

    indices_A = find(Result.Spot == "hot" & Result.Muscle == "ADM" & Result.Session == "S2" & Result.Dominant == dominant);
    indices_B = arrayfun(@(i) find(Result.Spot == "hot" & Result.Subject == Result.Subject(i) & Result.Dominant == Result.Dominant(i) & Result.Muscle == "APB" & Result.Session == "S2"), indices_A);
    differences = [Result.X_source(indices_B) - Result.X_source(indices_A), Result.Y_source(indices_B) - Result.Y_source(indices_A), Result.Z_source(indices_B) - Result.Z_source(indices_A)];
    scoring = differences * anterolateral_axis';
    distances.(dominance).ADMvsAPB = scoring; % sqrt((Result.X_source(indices_A) - Result.X_source(indices_B)).^2 + (Result.Y_source(indices_A) - Result.Y_source(indices_B)).^2 + (Result.Z_source(indices_A) - Result.Z_source(indices_B)).^2);

    indices_A = find(Result.Spot == "hot" & Result.Muscle == "FDI" & Result.Session == "S2" & Result.Dominant == dominant);
    indices_B = arrayfun(@(i) find(Result.Spot == "hot" & Result.Subject == Result.Subject(i) & Result.Dominant == Result.Dominant(i) & Result.Muscle == "APB" & Result.Session == "S2"), indices_A);
    differences = [Result.X_source(indices_B) - Result.X_source(indices_A), Result.Y_source(indices_B) - Result.Y_source(indices_A), Result.Z_source(indices_B) - Result.Z_source(indices_A)];
    scoring = differences * anterolateral_axis';
    distances.(dominance).FDIvsAPB = scoring; %sqrt((Result.X_source(indices_A) - Result.X_source(indices_B)).^2 + (Result.Y_source(indices_A) - Result.Y_source(indices_B)).^2 + (Result.Z_source(indices_A) - Result.Z_source(indices_B)).^2);

    indices_A = find(Result.Spot == "hot" & Result.Muscle == "FDI" & Result.Session == "S2" & Result.Dominant == dominant);
    indices_B = arrayfun(@(i) find(Result.Spot == "cold" & Result.Subject == Result.Subject(i) & Result.Dominant == Result.Dominant(i) & Result.Muscle == Result.Muscle(i) & Result.Session == "S2"), indices_A);
    differences = [Result.X_source(indices_B) - Result.X_source(indices_A), Result.Y_source(indices_B) - Result.Y_source(indices_A), Result.Z_source(indices_B) - Result.Z_source(indices_A)];
    scoring = differences * compare_along';
    distances.(dominance).hotcold = scoring; %sqrt((Result.X_source(indices_A) - Result.X_source(indices_B)).^2 + (Result.Y_source(indices_A) - Result.Y_source(indices_B)).^2 + (Result.Z_source(indices_A) - Result.Z_source(indices_B)).^2);
    coldspot_R2.(dominance) = Result.R2(indices_B);

    connector = [Result.X_source(indices_B) - Result.X_source(indices_A), Result.Y_source(indices_B) - Result.Y_source(indices_A), Result.Z_source(indices_B) - Result.Z_source(indices_A)];


    scoring = connector * compare_along';
    distances.(dominance).hotcold_posterolateral_component = scoring;

    scoring_found_coldspot = scoring(coldspot_R2.(dominance) > 0.1);
    p_projected = signrank(scoring_found_coldspot, 0, "tail", "right");
    
    fprintf("%s hemisphere (%s trials):\t Shifted by %3.3f ± %3.3f mm posterolaterally on ind. brain    p = %3.3f\n", dominance, SPLIT, mean(scoring_found_coldspot), std(scoring_found_coldspot), p_projected)
    
end

fprintf("\n\n\n")











%%



addpath(genpath('//wsl.localhost/Ubuntu-22.04/usr/local/freesurfer/matlab'))
addpath(genpath('B:/Projects/2023-01 HighDef/libraries/visualization'))
fs_subject_dir = '//wsl.localhost/Ubuntu-22.04/usr/local/freesurfer/subjects'; %getenv('SUBJECTS_DIR');

exp_names = []; exp_names.l = 'map-L2R'; exp_names.r = 'map-R2L';
response_names = []; response_names.hot = 'CsE_FDI_in_uV'; response_names.cold = 'SIHIscore_FDI';


[vertex_coords_lh, faces_lh] = read_surf([fs_subject_dir '/fsaverage/surf/lh.inflated']);
[vertex_coords_rh, faces_rh] = read_surf([fs_subject_dir '/fsaverage/surf/rh.inflated']);

curv_lh = read_curv([fs_subject_dir '/fsaverage/surf/lh.curv']);
curv_rh = read_curv([fs_subject_dir '/fsaverage/surf/rh.curv']);

vertex_coords = []; vertex_coords.l = vertex_coords_lh; vertex_coords.r = vertex_coords_rh;
faces = []; faces.l = faces_lh; faces.r = faces_rh;
curvature.l = curv_lh; curvature.r = curv_rh;
curvature.l(curvature.l > 0) = 1;
curvature.l(curvature.l <= 0) = -1;
curvature.r(curvature.r > 0) = 1;
curvature.r(curvature.r <= 0) = -1;

results = [];


subjects_sublist = 1:length(subjects);

distances.dominant.avg_hotcold = nan(length(distances.dominant.hotcold), 1);
distances.nondominant.avg_hotcold = nan(length(distances.nondominant.hotcold), 1);





for iSubject = subjects_sublist
    subject = subjects(iSubject);
    for hemisphere = ["l", "r"]
        for spot = ["hot", "cold"]
            try
                mgh_data = MRIread(sprintf('%s/%s/projected_%s_%s_%s%s_vertmax.mgh', STORAGE, subject, subject, exp_names.(hemisphere), response_names.(spot), suffix));
                [maxR2, where_max] = max(mgh_data.vol);
                results.(hemisphere).(spot).R2(iSubject) = maxR2;
    
                max_coords = vertex_coords.(hemisphere)(where_max,:);
                results.(hemisphere).(spot).x(iSubject) = max_coords(1);
                results.(hemisphere).(spot).y(iSubject) = max_coords(2);
                results.(hemisphere).(spot).z(iSubject) = max_coords(3);
            catch 
                fprintf('Missing backprojected data for subject %s/%s/%s!\n', subject, hemisphere, spot)
                results.(hemisphere).(spot).x(iSubject) = nan; results.(hemisphere).(spot).y(iSubject) = nan; results.(hemisphere).(spot).z(iSubject) = nan; results.(hemisphere).(spot).R2(iSubject) = nan;
            end
        end
        dif = nan(3,1);
        dimnames = ["x", "y", "z"];
        for i = 1:3
            dif(i) = results.(hemisphere).cold.(dimnames(i))(iSubject) - results.(hemisphere).hot.(dimnames(i))(iSubject);
        end

        if ismember(subject, leftHanded)
            if hemisphere == "l"
                compare_along = [-1 -1 -1] ./ sqrt(3);
                distances.nondominant.avg_hotcold(iSubject) = compare_along * dif; %sqrt(sum(dif.^2));
            else
                compare_along = [1 -1 -1] ./ sqrt(3);
                distances.dominant.avg_hotcold(iSubject) = compare_along * dif; % sqrt(sum(dif.^2));
            end
        else
            if hemisphere == "l"
                compare_along = [-1 -1 -1] ./ sqrt(3);
                distances.dominant.avg_hotcold(iSubject) = compare_along * dif; % sqrt(sum(dif.^2));
            else
                compare_along = [1 -1 -1] ./ sqrt(3);
                distances.nondominant.avg_hotcold(iSubject) = compare_along * dif; % sqrt(sum(dif.^2));
            end
        end
    end
    fprintf('Retrieved data for %s\n', subject)
end



R2_cutoff = 0.1;


%% Compute the posterolateral component of the distance and test if significant on fs_avg

fprintf("\n\nTest difference of spots on fs average brain for SSAM (%s)\n", SPLIT)
for hemisphere = ["dominant", "nondominant"]
    mirror_which = ismember(subjects, leftHanded);
    if hemisphere == "dominant"
        this = "l";
        other = "r";
        compare_along = [-1 -1 -1] ./ sqrt(3); % downward posterolateral (more along the surface)
    else
        this = "r";
        other = "l";
        compare_along = [1 -1 -1] ./ sqrt(3);
    end

    hotspot_locations  = [results.(this).hot.x;  results.(this).hot.y;  results.(this).hot.z]';
    coldspot_locations = [results.(this).cold.x; results.(this).cold.y; results.(this).cold.z]';
    coldspot_R2s = results.(this).cold.R2;

    hotspot_locations(mirror_which,:)  = [-results.(other).hot.x(mirror_which)',  results.(other).hot.y(mirror_which)',  results.(other).hot.z(mirror_which)'];
    coldspot_locations(mirror_which,:) = [-results.(other).cold.x(mirror_which)', results.(other).cold.y(mirror_which)', results.(other).cold.z(mirror_which)'];
    coldspot_R2s(mirror_which) = results.(other).cold.R2(mirror_which);

    differences = coldspot_locations - hotspot_locations;
    % project onto postero-lateral unit vector:
    scoring = differences * compare_along';
    scoring_found_coldspot = scoring(coldspot_R2s > 0.1);
    p_projected = signrank(scoring_found_coldspot, 0, "tail", "right");
    fprintf("%s hemisphere (%s trials):\t Shifted by %3.3f ± %3.3f mm posterolaterally on avg. brain    p = %3.3f\n", hemisphere, SPLIT, mean(scoring_found_coldspot), std(scoring_found_coldspot), p_projected)
end

fprintf("\n\n\n")

%%
all_max = max([distances.dominant.hotcold; distances.nondominant.hotcold; distances.dominant.ADMvsFDI; distances.dominant.avg_hotcold; distances.nondominant.ADMvsFDI; distances.dominant.ADMvsAPB; distances.nondominant.ADMvsAPB; distances.dominant.FDIvsAPB; distances.nondominant.FDIvsAPB; distances.nondominant.avg_hotcold]);
all_min = min([distances.dominant.hotcold; distances.nondominant.hotcold; distances.dominant.ADMvsFDI; distances.dominant.avg_hotcold; distances.nondominant.ADMvsFDI; distances.dominant.ADMvsAPB; distances.nondominant.ADMvsAPB; distances.dominant.FDIvsAPB; distances.nondominant.FDIvsAPB; distances.nondominant.avg_hotcold]);

x0 = 0.2;
dotpos = []; dotpos.nondominant = x0+0.7; dotpos.dominant = x0+0.5;
boxpos = []; boxpos.nondominant = x0+1; boxpos.dominant = x0+0.2;
boxwidth = 0.4;
boxcolor = []; boxcolor.nondominant = [154, 114, 247]./255; boxcolor.dominant = [221, 247, 114]./255;

cerulean = [0, 150, 255]./255;


names = []; 
names.hotcold = ["SSAM";"ind. brain"]; 
names.ADMvsFDI = ["ADM \leftrightarrow FDI";"{\color{red}hot}spot"]; 
names.ADMvsAPB = ["ADM \leftrightarrow APB";"{\color{red}hot}spot"]; 
names.FDIvsAPB = ["FDI \leftrightarrow APB";"{\color{red}hot}spot"];
names.avg_hotcold = ["SSAM";"avg. brain"];
names.avg_simple_hotcold = ["Projection";"avg. brain"];






fig = figure(Position=[50 50 700 800]);
layout = tiledlayout(6,4, TileSpacing="compact");

ax_avg_brain = nexttile(1, [4 4]);
offset = []; offset.mag = 42; offset.l = -offset.mag; offset.r = offset.mag;
%layout = tiledlayout(1,2);
%nexttile;
for hemisphere = ["l", "r"]
    ts = trisurf(faces.(hemisphere)+1, vertex_coords.(hemisphere)(:,1) + offset.(hemisphere), vertex_coords.(hemisphere)(:,2), vertex_coords.(hemisphere)(:,3), curvature.(hemisphere)); 
    ts.EdgeColor='none'; 
    %ts.FaceColor = [0.98 1 0.995];
    %colormap([0.95 0.95 0.94;0.8 0.81 0.8]);
    colormap(0.8.*[0.355 0.35 0.37; 0.215 0.21 0.23]);
    hold on
    coldmask = results.(hemisphere).cold.R2 > R2_cutoff;
    hotmask = results.(hemisphere).hot.R2 > R2_cutoff;
    plot3([results.(hemisphere).hot.x; results.(hemisphere).cold.x] + offset.(hemisphere), [results.(hemisphere).hot.y; results.(hemisphere).cold.y], [results.(hemisphere).hot.z; results.(hemisphere).cold.z]+5, 'w.-', LineWidth=1.4, MarkerSize=8)
    plot3(results.(hemisphere).hot.x(hotmask) + offset.(hemisphere), results.(hemisphere).hot.y(hotmask), results.(hemisphere).hot.z(hotmask)+10, 'r^', MarkerFaceColor='r')
    plot3(results.(hemisphere).cold.x(coldmask) + offset.(hemisphere), results.(hemisphere).cold.y(coldmask), results.(hemisphere).cold.z(coldmask)+10, ...
        'v', MarkerFaceColor=cerulean, MarkerEdgeColor=cerulean)
end


%light(Style="infinite",Position=[0 50 80],Color=[1 0.98 0.97]);
%light(Style="infinite",Position=[0 -10 80],Color=[1 0.98 0.97]);
c = [0.15 0.3 0.8]; %[0.07 0.2 0.8];
%light(Style="infinite",Position=[-10 -10 -10],Color=0.75.*c);
%light(Style="infinite",Position=[10 -10 -10],Color=0.75.*c);

lr = max([range(xlim), range(ylim), range(zlim)]);
lr = 220;
xlim(mean(xlim) + 0.5.*[-lr lr]); ylim(mean(ylim) + 0.5.*[-lr lr]); zlim(mean(zlim) + 0.5.*[-lr lr])
axis square
axis off
view(2)
exportgraphics(ax_avg_brain, sprintf("B:/Projects/2023-01 HighDef/Results/Evaluation/spots_on_avg_brain%s.png", suffix), Resolution=900)
text(min(xlim), max(ylim), max(zlim), "A", HorizontalAlignment="left", FontSize=14, FontWeight="bold")





% Load "conventional evaluation" --- make sure to run
% compare_SSAM_and_projection.m before
ComparisonTable = readtable(sprintf("signed-distances%s.csv", suffix));
distances.dominant.avg_simple_hotcold = ComparisonTable.simple_dominant;
distances.nondominant.avg_simple_hotcold = ComparisonTable.simple_nondominant;


collectedSubjects = ["sub-001", "sub-002", "sub-003", "sub-004", "sub-006", "sub-007", "sub-008", "sub-009", "sub-011", "sub-012", "sub-013", "sub-014", "sub-015", "sub-016", "sub-017", "sub-019", "sub-022", "sub-023"];

comparison_names = fieldnames(distances.dominant)';
comparison_names = comparison_names(~strcmpi(comparison_names, "hotcold_posterolateral_component"));
nr = length(comparison_names);
subplot_indices = []; subplot_indices.avg_simple_hotcold = 0; subplot_indices.hotcold = 2; subplot_indices.avg_hotcold = 1; subplot_indices.ADMvsAPB = 2; subplot_indices.FDIvsAPB = 0; subplot_indices.ADMvsFDI = 1;

plot_spacing = 1.6;

ax_cold = nexttile(17, [2 2]);
hold on
ax_hot = nexttile(19, [2 2]);
hold on

for i = 1:nr
    comp = comparison_names{i};
    if endsWith(comp, "hotcold")
        ax = ax_cold;
    else
        ax = ax_hot;    
    end

    axes(ax);
    hold on

    for dominant = [true, false]
        if dominant
            h = "dominant";
        else
            h = "nondominant";
        end
        
        d = distances.(h).(comp);
        fprintf("%s %s:\n\tData:  %s\n", h, comp, strjoin(cellfun(@num2str, num2cell(sort(d)'), 'UniformOutput', false)))
        labels = collectedSubjects';

        x0 = subplot_indices.(comp)*plot_spacing;

        ms = 2;
        if strcmpi(comp, 'hotcold') 
            found_coldspot = coldspot_R2.(h) > 0.1;
            plot(x0+dotpos.(h), d(found_coldspot), 'o', Color=hemisphere_colors.(h), MarkerFaceColor=hemisphere_colors.(h), MarkerSize=ms)
            plot(x0+dotpos.(h), d(~found_coldspot), 'o', Color=hemisphere_colors.(h), MarkerSize=ms)
            d = d(found_coldspot);
            labels = labels(found_coldspot);
        elseif strcmpi(comp, 'avg_hotcold')
            if dominant
                found_coldspot = results.l.cold.R2 > 0.1;
                found_coldspot(subjects == "sub-006") = results.r.cold.R2(subjects == "sub-006") > 0.1;
                found_coldspot(subjects == "sub-016") = results.r.cold.R2(subjects == "sub-016") > 0.1;
                found_coldspot(subjects == "sub-017") = results.r.cold.R2(subjects == "sub-017") > 0.1;
            else
                found_coldspot = results.r.cold.R2 > 0.1;
                found_coldspot(subjects == "sub-006") = results.l.cold.R2(subjects == "sub-006") > 0.1;
                found_coldspot(subjects == "sub-016") = results.l.cold.R2(subjects == "sub-016") > 0.1;
                found_coldspot(subjects == "sub-017") = results.l.cold.R2(subjects == "sub-017") > 0.1;
            end
            plot(x0+dotpos.(h), d(found_coldspot), 'o', Color=hemisphere_colors.(h), MarkerFaceColor=hemisphere_colors.(h), MarkerSize=ms)
            plot(x0+dotpos.(h), d(~found_coldspot), 'o', Color=hemisphere_colors.(h), MarkerSize=ms)
            labels = labels(found_coldspot(~isnan(d)));
            d = d(found_coldspot);
        else
            plot(x0+dotpos.(h), d, 'o', Color=hemisphere_colors.(h), MarkerFaceColor=hemisphere_colors.(h), MarkerSize=ms)
        end
        
        xb0 = x0 + boxpos.(h);

        writetable(table(labels, d, VariableNames=["Subject", "Distance"]), sprintf("B:/Projects/2023-01 HighDef/Results/Coil-space/ssam%s-%s-%s.csv", suffix, h, comp))
        
        q1 = quantile(d, 0.25); q3 = quantile(d, 0.75);
        w1 = min(d(d >= q1-(1.5*(q3-q1))));
        w3 = max(d(d <= q3+(1.5*(q3-q1))));
        fprintf("\tQ1 = %3.3f\n", q1)
        fprintf("\tQ3 = %3.3f\n", q3)
        fprintf("\tW1 = %3.3f\n", w1)
        fprintf("\tW3 = %3.3f\n", w3)
        fprintf("\tmedian = %3.3f mm\n", median(d, "omitmissing"))
        fprintf("\tmean = %3.3f ± %3.3f mm\n", mean(d, "omitmissing"), std(d, "omitmissing"))
        patch(xb0 + boxwidth.*[-0.5 0.5 0.5 -0.5 -0.5], [q1 q1 q3 q3 q1], boxcolor.(h), FaceAlpha=1, DisplayName=h)
        plot(xb0 + boxwidth.*[-0.5 0.5], median(d, "omitmissing") * ones(1,2), '-', Color=hemisphere_colors.(h), LineWidth=2)
        plot(xb0 .* ones(1,2), [q1 w1], 'k-'); plot(xb0 + boxwidth.*[-0.25 0.25], [w1 w1], 'k-')
        plot(xb0 .* ones(1,2), [q3 w3], 'k-'); plot(xb0 + boxwidth.*[-0.25 0.25], [w3 w3], 'k-')
        plot([0 6], [0 0], "k:", LineWidth=0.2)
        %text(xb0, all_min-2, h, HorizontalAlignment="right", Color=hemisphere_colors.(h), Rotation=90)
    end

    text(x0 + 0.5*(boxpos.dominant + boxpos.nondominant), all_min-2, names.(comp), HorizontalAlignment="center", VerticalAlignment="top")
    if endsWith(comp, "hotcold")
        text(x0 + 0.5*(boxpos.dominant + boxpos.nondominant), all_max+3, sprintf("B%d", subplot_indices.(comp)+1), FontWeight="bold", HorizontalAlignment="center", Color=hemisphere_colors.(h))
    end

    set(ax, Color="none");
    ax.XAxis.Visible = "off";
    ylim([all_min all_max])
    xlim([0 3*plot_spacing])
end

% Left hemisphere is dominant (and thus above right hemisphere)
axes(ax_cold);
ylabel("Posterolateral distance [mm]", FontSize=12)
text(0, all_max+3, "B", HorizontalAlignment="center", FontSize=14, FontWeight="bold")
text(1.5*plot_spacing, all_min-12, sprintf("FDI {\\color{red}hot} \\leftrightarrow {\\color[rgb]{%3.5f,%3.5f,%3.5f}cold}", cerulean), VerticalAlignment="top", HorizontalAlignment="center", FontWeight="bold")
axes(ax_hot);
ylabel("Anterolateral distance [mm]", FontSize=12)
text(0, all_max+3, "C", HorizontalAlignment="center", FontSize=14, FontWeight="bold")
text(1.5*plot_spacing, all_min-12, "SSAM, ind. brain", VerticalAlignment="top", HorizontalAlignment="center", FontWeight="bold")


legend_labels = repmat("", 32, 0);
legend_labels(length(subjects) + 1) = "dominant";
legend_labels(8 + 2*length(subjects)) = "nondominant";
legend(legend_labels, FontSize=10)

y = layout.Position(2);
h = layout.Position(4);
layout.Position(2) = y + 0.04;
layout.Position(4) = h - 0.04;

exportgraphics(fig, sprintf('B:/Projects/2023-01 HighDef/Results/Evaluation/group_distances_and_avg_brain%s.pdf', suffix), Resolution=600)
exportgraphics(fig, sprintf('B:/Projects/2023-01 HighDef/Results/Evaluation/group_distances_and_avg_brain%s.png', suffix), Resolution=950)
















% 
% 
% %% 2 Plot the comparisons for selected comparison (hot in S1 vs S2 intra muscle, hot muscle vs muscle intra session, hot vs cold in S2 intra muscle)
% 
% % in the end, want data of size quantity x muscles x subjects
% 
% %% Reliability of hotspots within each muscle (R)
% muscleContrasts = [muscles; muscles];
% 
% hemisphere = "L";
% f1L = plot_contrasts(comparant(hemisphere, "S1", "", "hot"), comparant(hemisphere, "S2", "", "hot"), muscleContrasts, Result, colors, markers);
% f1L.Name = sprintf("Hotspot reliability within muscle %s", hemisphere);
% 
% hemisphere = "R";
% f1R = plot_contrasts(comparant(hemisphere, "S1", "", "hot"), comparant(hemisphere, "S2", "", "hot"), muscleContrasts, Result, colors, markers);
% f1R.Name = sprintf("Hotspot reliability within muscle %s", hemisphere);
% 
% %% Specificity
% muscleContrasts = ["ADM" "FDI"; "ADM" "APB"; "FDI" "APB"]';
% 
% hemisphere = "L";
% f2S1L = plot_contrasts(comparant(hemisphere, "S1", "", "hot"), comparant(hemisphere, "S1", "", "hot"), muscleContrasts, Result, colors, markers);
% f2S1L.Name = sprintf("Distance of hotspots (%s) in S1", hemisphere);
% 
% f2S2L = plot_contrasts(comparant(hemisphere, "S2", "", "hot"), comparant(hemisphere, "S2", "", "hot"), muscleContrasts, Result, colors, markers);
% f2S2L.Name = sprintf("Distance of hotspots (%s) in S2", hemisphere);
% 
% hemisphere = "R";
% f2S1R = plot_contrasts(comparant(hemisphere, "S1", "", "hot"), comparant(hemisphere, "S1", "", "hot"), muscleContrasts, Result, colors, markers);
% f2S1R.Name = sprintf("Distance of hotspots (%s) in S1", hemisphere);
% 
% f2S2R = plot_contrasts(comparant(hemisphere, "S2", "", "hot"), comparant(hemisphere, "S2", "", "hot"), muscleContrasts, Result, colors, markers);
% f2S2R.Name = sprintf("Distance of hotspots (%s) in S2", hemisphere);
% 
% 
% %% Hotspot vs. Coldspot
% muscleContrasts = [muscles; muscles];
% 
% hemisphere = "L";
% f3L = plot_contrasts(comparant(hemisphere, "S2", "", "hot"), comparant(hemisphere, "S2", "", "cold"), muscleContrasts, Result, colors, markers);
% f3L.Name = sprintf("Distance between hot and cold for each muscle (%s)", hemisphere);
% 
% hemisphere = "R";
% f3R = plot_contrasts(comparant(hemisphere, "S2", "", "hot"), comparant(hemisphere, "S2", "", "cold"), muscleContrasts, Result, colors, markers);
% f3R.Name = sprintf("Distance between hot and cold for each muscle (%s)", hemisphere);
% 
% 
% %% Coldspot muscle specificity
% muscleContrasts = ["ADM" "FDI"; "ADM" "APB"; "FDI" "APB"]';
% 
% hemisphere = "L";
% f4L = plot_contrasts(comparant(hemisphere, "S2", "", "cold"), comparant(hemisphere, "S2", "", "cold"), muscleContrasts, Result, colors, markers);
% f4L.Name = sprintf("Distance of coldspots (%s) in S2", hemisphere);
% 
% hemisphere = "R";
% f4R = plot_contrasts(comparant(hemisphere, "S2", "", "cold"), comparant(hemisphere, "S2", "", "cold"), muscleContrasts, Result, colors, markers);
% f4R.Name = sprintf("Distance of coldspots (%s) in S2", hemisphere);
% 
% 
% 










function [fig] = plot_contrasts(c1, c2, muscleContrasts, data, colors, markers)
fig = figure(Position=[50 150 420 450]);
tiledlayout(3,2)

all_c = [];
for muscle = muscleContrasts
    c1.Muscle = muscle(1);
    c2.Muscle = muscle(end);

    condition = join(unique(muscle), "vs");
    c = make_comparison_for_each_subject(data, c1, c2);
    c.Condition = repmat([condition], size(c,1), 1);
    all_c = [all_c; c];
    unreliable = c.quality < (1/8)^2; % arbitrary threshold! May also need to depend on the muscle

    nexttile(2)
    hold on
    plot(c(~unreliable,:), "coil_distance", "source_distance", Marker=markers.(condition), LineStyle='none', Color=colors.(condition), LineWidth=2)
    plot(c(unreliable,:), "coil_distance", "source_distance", Marker=markers.(condition), LineStyle='none', Color=colors.(condition), LineWidth=0.5, MarkerSize=5)
    set(gca, "Color", "none"); axis square; box on

    nexttile(4)
    hold on
    plot(c(~unreliable,:), "coil_distance", "angle_distance_deg", Marker=markers.(condition), LineStyle='none', Color=colors.(condition), LineWidth=2)
    plot(c(unreliable,:), "coil_distance", "angle_distance_deg", Marker=markers.(condition), LineStyle='none', Color=colors.(condition), LineWidth=0.5, MarkerSize=5)
    set(gca, "Color", "none"); axis square; box on

end

nexttile(1)
b = boxchart(all_c.source_distance, GroupByColor=all_c.Condition);
for ib = 1:length(b)
    b(ib).BoxFaceColor = colors.(b(ib).DisplayName);
end
xticks([])
set(gca, "Color", "none"); axis square; box on

nexttile(3)
b = boxchart(all_c.angle_distance_deg, GroupByColor=all_c.Condition);
for ib = 1:length(b)
    b(ib).BoxFaceColor = colors.(b(ib).DisplayName);
end
xticks([])
set(gca, "Color", "none"); axis square; box on

nexttile(6)
b = boxchart(all_c.coil_distance, Orientation="horizontal", GroupByColor=all_c.Condition);
for ib = 1:length(b)
    b(ib).BoxFaceColor = colors.(b(ib).DisplayName);
    %b(ib).BoxMedianLineColor = [0 0 0];
    %b(ib).BoxEdgeColor = [0 0 0];
end
yticks([])
set(gca, "Color", "none"); axis square; box on
legend(Location="westoutside")
end






function comparisons = make_comparison_for_each_subject(data, c1, c2)
result_columns = ["source_distance", "coil_distance", "angle_distance_deg", "quality"];
columns_source = ["X_source", "Y_source", "Z_source"];
coil_pos_columns = ["X_coil", "Y_coil", "Z_coil"];

mask1 = true(size(data, 1), 1);
mask2 = true(size(data, 1), 1);
for cf = fieldnames(c1)'
    f = cf{:};
    mask1 = mask1 & data.(f) == c1.(f);
    mask2 = mask2 & data.(f) == c2.(f);
end

subjects = unique(data.Subject)';
nSubjects = length(subjects);
comparisons = [table(subjects', VariableNames=["Subject"]) array2table(nan(nSubjects, length(result_columns)), VariableNames=result_columns)];
for iSubject = 1:nSubjects
    subject = subjects(iSubject);
    smask1 = data.Subject == subject & mask1;
    smask2 = data.Subject == subject & mask2;
    if any(smask1) && any(smask2)
        comparisons.source_distance(iSubject)    = sqrt(sum((table2array(data(smask1, columns_source) - data(smask2, columns_source))).^2, 2));
        comparisons.coil_distance(iSubject)      = sqrt(sum((table2array(data(smask1, coil_pos_columns) - data(smask2, coil_pos_columns))).^2, 2));
        comparisons.angle_distance_deg(iSubject) = table2array(data(smask1, "Angle_in_rad") - data(smask2, "Angle_in_rad"));
        comparisons.quality(iSubject) = data.R2(smask1) .* data.R2(smask2);
    end
end

% Wrap angle differences
over_plus_180 = logical(table2array(comparisons(:, "angle_distance_deg") > pi));
comparisons(over_plus_180, "angle_distance_deg") = comparisons(over_plus_180, "angle_distance_deg") - 2*pi;
under_minus_180 = logical(table2array(comparisons(:, "angle_distance_deg") < -pi));
comparisons(under_minus_180, "angle_distance_deg") = comparisons(under_minus_180, "angle_distance_deg") + 2*pi;
comparisons.angle_distance_deg = rad2deg(table2array(comparisons(:, "angle_distance_deg")));
end


function cfg = comparant(hemisphere, session, muscle, spot)
cfg = [];
cfg.Hemisphere = hemisphere;
cfg.Session = session;
cfg.Muscle = muscle;
cfg.Spot = spot;
end





