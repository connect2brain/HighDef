ID = 'sub-001';

roi_geo = sprintf('//wsl.localhost/Ubuntu-22.04/home/bnplab-admin/TMS_localization/TMS_loc_results/%s/mesh/charm/roi/midlayer_m1s1pmd/geo.hdf5', ID);
all_geo = sprintf('//wsl.localhost/Ubuntu-22.04/home/bnplab-admin/TMS_localization/TMS_loc_results/%s/mesh/charm/m2m_%s/%s.hdf5', ID, ID, ID);

roi_coordinates = h5read(roi_geo, '/mesh/nodes/node_coord');
roi_triangles = h5read(roi_geo, '/mesh/elm/triangle_number_list')' + 1;

all_coordinates = h5read(all_geo, '/mesh/nodes/node_coord');
all_triangles = h5read(all_geo, '/mesh/elm/triangle_number_list')' + 1;
all_tissue_type = h5read(all_geo, '/mesh/elm/tri_tissue_type')';

wm_triangles = all_triangles(all_tissue_type == 1001,:);



%%
ts = trisurf(wm_triangles, all_coordinates(1,:), all_coordinates(2,:), all_coordinates(3,:));
ts.EdgeColor = 'none';
ts.FaceColor = 0.2 .* ones(1,3);
hold on
ts_roi = trisurf(roi_triangles, roi_coordinates(1,:), roi_coordinates(2,:), roi_coordinates(3,:));
ts_roi.FaceColor = [0.2 0.6 1];
ts_roi.EdgeColor = 'none';
axis vis3d
light
lighting gouraud
material dull


%% Golden brain:
gm_triangles = all_triangles(all_tissue_type == 1005,:);
ts = trisurf(gm_triangles, all_coordinates(1,:), all_coordinates(2,:), all_coordinates(3,:));
ts.EdgeColor = 'none';
ts.FaceColor = [0.98 0.82 0.5];
axis equal
axis vis3d
light(Style="infinite",Position=[-10 0 10],Color="white");
light(Style="infinite",Position=[0 0 -10],Color=[0.1 0.4 0.8]);
material shiny
set(gca, 'Color', 'none')
grid off
axis off
