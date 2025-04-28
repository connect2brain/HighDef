addpath(genpath('//wsl.localhost/Ubuntu-22.04/usr/local/freesurfer/7.4.1/matlab'))
addpath(genpath('B:/Projects/2023-01 HighDef/libraries/visualization'))
fs_subject_dir = '//wsl.localhost/Ubuntu-22.04/usr/local/freesurfer/7.4.1/subjects'; %getenv('SUBJECTS_DIR');

exp_names = []; exp_names.l = 'map-L2R'; exp_names.r = 'map-R2L';
response_names = []; response_names.hot = 'CsE_FDI_in_uV'; response_names.cold = 'SIHIscore_FDI';


[vertex_coords_lh, faces_lh] = read_surf([fs_subject_dir '/fsaverage/surf/lh.inflated']);
[vertex_coords_rh, faces_rh] = read_surf([fs_subject_dir '/fsaverage/surf/rh.inflated']);

vertex_coords = []; vertex_coords.l = vertex_coords_lh; vertex_coords.r = vertex_coords_rh;
faces = []; faces.l = faces_lh; faces.r = faces_rh;

results = [];


subjects = ["sub-001", "sub-003", "sub-004", "sub-006", "sub-007", "sub-008", "sub-009"];

for iSubject = 1:length(subjects)
    subject = subjects(iSubject);
    for hemisphere = ["l", "r"]
        for spot = ["hot", "cold"]
            mgh_data = MRIread(sprintf('U:/home/bnplab-admin/TMS_localization/HighDef/%s/projected_%s_%s_%s_vertmax.mgh', subject, subject, exp_names.(hemisphere), response_names.(spot)));
            [maxR2, where_max] = max(mgh_data.vol);
            results.(hemisphere).(spot).R2(iSubject) = maxR2;

            max_coords = vertex_coords.(hemisphere)(where_max,:);
            results.(hemisphere).(spot).x(iSubject) = max_coords(1);
            results.(hemisphere).(spot).y(iSubject) = max_coords(2);
            results.(hemisphere).(spot).z(iSubject) = max_coords(3);
        end
    end
    fprintf('Retrieved data for %s\n', subject)
end








%%

R2_cutoff = 0.05;

fig = figure(Position=[50 50 700 670]);
offset = []; offset.mag = 42; offset.l = -offset.mag; offset.r = offset.mag;
%layout = tiledlayout(1,2);
%nexttile;
for hemisphere = ["l", "r"]
    ts = trisurf(faces.(hemisphere)+1, vertex_coords.(hemisphere)(:,1) + offset.(hemisphere), vertex_coords.(hemisphere)(:,2), vertex_coords.(hemisphere)(:,3)); 
    ts.EdgeColor='none'; 
    ts.FaceColor = [0.98 1 0.995];
    hold on
    coldmask = results.(hemisphere).cold.R2 > R2_cutoff;
    hotmask = results.(hemisphere).hot.R2 > R2_cutoff;
    plot3([results.(hemisphere).hot.x; results.(hemisphere).cold.x] + offset.(hemisphere), [results.(hemisphere).hot.y; results.(hemisphere).cold.y], [results.(hemisphere).hot.z; results.(hemisphere).cold.z]+5, 'k.-')
    plot3(results.(hemisphere).hot.x(hotmask) + offset.(hemisphere), results.(hemisphere).hot.y(hotmask), results.(hemisphere).hot.z(hotmask)+5, 'r^', MarkerFaceColor='r')
    plot3(results.(hemisphere).cold.x(coldmask) + offset.(hemisphere), results.(hemisphere).cold.y(coldmask), results.(hemisphere).cold.z(coldmask)+5, 'bv', MarkerFaceColor='b')
end


light(Style="infinite",Position=[0 50 80],Color=[1 0.98 0.97]);
light(Style="infinite",Position=[0 -10 80],Color=[1 0.98 0.97]);
c = [0.15 0.3 0.8]; %[0.07 0.2 0.8];
light(Style="infinite",Position=[-10 -10 -10],Color=0.75.*c);
light(Style="infinite",Position=[10 -10 -10],Color=0.75.*c);

lr = max([range(xlim), range(ylim), range(zlim)]);
xlim(mean(xlim) + 0.5.*[-lr lr]); ylim(mean(ylim) + 0.5.*[-lr lr]); zlim(mean(zlim) + 0.5.*[-lr lr])
axis square
axis off
view(2)

exportgraphics(fig, 'B:/Projects/2023-01 HighDef/Results/AllSpots-Figures/avg_brain.png', Resolution=300)