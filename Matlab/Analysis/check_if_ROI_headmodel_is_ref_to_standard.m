for ID = ["sub-001", "sub-002", "sub-003", "sub-004", "sub-006", "sub-007", "sub-008", "sub-009"]

    exp_id = "map-R2L";
    hemisphere = extractAfter(extractBefore(exp_id, "2"), "-");
    mesh_id = 'mesh0';
    roi = sprintf('midlayer_%s', lower(hemisphere));

    geo_file = sprintf('//wsl.localhost/Ubuntu-22.04/home/bnplab-admin/TMS_localization/HighDef/%s/mesh/roi/%s/geo.hdf5', ID, roi);
    coordinates = h5read(geo_file, '/mesh/nodes/node_coord');
    triangles = h5read(geo_file, '/mesh/elm/triangle_number_list')' + 1;
    fprintf("%s has %5.f triangles\n", ID, size(triangles,1))
end

% Conclusion: The ROIs do not match even in size!


%%
ID = "sub-001";
exp_id = "map-L2R";
roi_id = "small_l";
mesh_id = "mesh0";

geo_file = sprintf('//wsl.localhost/Ubuntu-22.04/home/bnplab-admin/TMS_localization/HighDef/%s/results/exp_%s/r2/mesh_%s/roi_%s/CsE_FDI_in_uV/sigmoid4/r2_geo.hdf5', ID, exp_id, mesh_id, roi_id);
coordinates = h5read(geo_file, '/mesh/nodes/node_coord');
triangles = h5read(geo_file, '/mesh/elm/triangle_number_list')' + 1;
tissue_type = h5read(geo_file, '/mesh/elm/tri_tissue_type')';

data_file = sprintf('//wsl.localhost/Ubuntu-22.04/home/bnplab-admin/TMS_localization/HighDef/%s/results/exp_%s/r2/mesh_%s/roi_%s/CsE_FDI_in_uV/sigmoid4/r2_data.hdf5', ID, exp_id, mesh_id, roi_id);
data = h5read(data_file, "/data/tris/c_E_mag");

trisurf(triangles, coordinates(1,:), coordinates(2,:), coordinates(3,:), data, EdgeColor="none")
axis square