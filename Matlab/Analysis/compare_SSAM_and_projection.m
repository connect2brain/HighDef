% [1] Retrieve Projection distances on the average brain
conventional_project_to_surface
close all
clc

STORAGE = "D:/HighDef-operate/HighDef";
leftHanded = ["sub-006", "sub-016", "sub-017"];

%%
%compare_along = [-1 -1 -1] ./ sqrt(3); % downward posterolateral (more along the surface)
direction.dominant = [-1 -1 -1];
direction.dominant = direction.dominant ./ norm(direction.dominant);
direction.nondominant = [1 -1 -1];
direction.nondominant = direction.nondominant ./ norm(direction.nondominant);


ComparisonTable = table();
ComparisonTable.Subject = subjects;

fprintf("Result for conventional:\n")

for hemisphere = ["dominant", "nondominant"]
    compare_along = direction.(hemisphere);

    scoring = offsets.(hemisphere) * compare_along';

    ComparisonTable.(sprintf("simple_%s", hemisphere)) = scoring;
    p_projected = signrank(scoring, 0, "tail", "right");
    fprintf("%s hemisphere:\t Shifted by %3.3f ± %3.3f mm posterolaterally     p = %3.3f\n", hemisphere, mean(scoring), std(scoring), p_projected)
end

%% [2] Retrieve SSAM distances on the average brain

SPLIT = "all"; % "all", "first_half", "second_half"
suffix = sprintf("_%s", SPLIT);
if strcmpi(SPLIT, "all")
    suffix = "";
end


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


names.avg_hotcold = ["FDI Hotspot"; "vs"; "FDI Coldspot (avg. brain)"];


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
    end
    fprintf('Retrieved data for %s\n', subject)
end



R2_cutoff = 0.1;

fprintf("SSAM on fs_avg:\n")
% Compute the posterolateral component of the distance and test if significant on fs_avg
for hemisphere = ["dominant", "nondominant"]
    mirror_which = ismember(subjects, leftHanded);
    if hemisphere == "dominant"
        this = "l";
        other = "r";
    else
        this = "r";
        other = "l";
    end

    compare_along = direction.(hemisphere);

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

    ComparisonTable.(sprintf("SSAM_%s", hemisphere)) = scoring;
    ComparisonTable.(sprintf("SSAM_%s", hemisphere))(coldspot_R2s < 0.1) = nan; % "signrank treats NaNs in x and y as missing values and ignores them."
end






%% [3] Compare (only, where R² > 0.1 for coldspot under SSAM)
bins = -30:5:30;
subplot(2,1,1)
histogram(ComparisonTable.SSAM_dominant, BinEdges=bins); hold on; histogram(ComparisonTable.simple_dominant, BinEdges=bins)
title("Dominant")

subplot(2,1,2)
histogram(ComparisonTable.SSAM_nondominant, BinEdges=bins); hold on; histogram(ComparisonTable.simple_nondominant, BinEdges=bins)
title("Non-Dominant")


p_dominant    = signrank(ComparisonTable.SSAM_dominant,    ComparisonTable.simple_dominant);
p_nondominant = signrank(ComparisonTable.SSAM_nondominant, ComparisonTable.simple_nondominant);

fprintf("\n\nComparison of SSAM (%s) and simple projection approach\n - dominant hemisphere:    p = %3.3f\n - nondominant hemisphere: p = %3.3f\n", SPLIT, p_dominant, p_nondominant)

writetable(ComparisonTable, sprintf("signed-distances%s.csv", suffix))

