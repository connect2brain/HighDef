
addpath(genpath('//wsl.localhost/Ubuntu-22.04/usr/local/freesurfer/7.4.1/matlab'))
addpath(genpath('B:/Projects/2023-01 HighDef/libraries/visualization'))
fs_subject_dir = '//wsl.localhost/Ubuntu-22.04/usr/local/freesurfer/7.4.1/subjects'; %getenv('SUBJECTS_DIR');

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

%%
offset = []; offset.mag = 42; offset.l = -offset.mag; offset.r = offset.mag;
for hemisphere = ["l", "r"]
    ts = trisurf(faces.(hemisphere)+1, vertex_coords.(hemisphere)(:,1) + offset.(hemisphere), vertex_coords.(hemisphere)(:,2), vertex_coords.(hemisphere)(:,3), curvature.(hemisphere)); 
    ts.EdgeColor='none'; 
    %ts.FaceColor = [0.98 1 0.995];
    %colormap([0.83 0.84 0.85; 0.98 1 0.995]);
    colormap([0.7 0.7 0.7; 0.85 0.85 0.85]);
    %ts.AmbientStrength = 0.8;
    hold on
end

%light(Style="infinite",Position=[-600 400 500],Color=0.2.*[1 0.98 0.97]);
%light(Style="infinite",Position=[600 400 500],Color=0.2.*[1 0.98 0.97]);
c = [0.15 0.3 0.8]; %[0.07 0.2 0.8];
%light(Style="infinite",Position=[-10 -10 -10],Color=0.75.*c);
%light(Style="infinite",Position=[10 -10 -10],Color=0.75.*c);

lr = max([range(xlim), range(ylim), range(zlim)]);
lr = 220;
xlim(mean(xlim) + 0.5.*[-lr lr]); ylim(mean(ylim) + 0.5.*[-lr lr]); zlim(mean(zlim) + 0.5.*[-lr lr])
axis square
axis off
%box on
view(2)