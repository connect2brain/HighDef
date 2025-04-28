geo_file = '//wsl.localhost/Ubuntu-22.04/home/bnplab-admin/TMS_localization/TMS_loc_results/sub-001/results/exp_main/r2/mesh0/roi_midlayer_m1s1pmd/CsE_in_uV/sigmoid4/r2_roi_geo.hdf5';
simulation_result = '//wsl.localhost/Ubuntu-22.04/home/bnplab-admin/TMS_localization/TMS_loc_results/sub-001/results/exp_main/electric_field/mesh0/roi_midlayer_m1s1pmd/e.hdf5';


h5disp(simulation_result)

%%

coordinates = h5read(geo_file, '/mesh/nodes/node_coord');
triangles = h5read(geo_file, '/mesh/elm/triangle_number_list')' + 1;
E_mag = h5read(simulation_result, '/E_mag');

maxE = max(E_mag(:));

%%
for i_trial = 1:size(E_mag, 2)
    l = 80;
    ts = trisurf(triangles, coordinates(1,:), coordinates(2,:), coordinates(3,:), E_mag(:,i_trial));
    xlim([-l l]); ylim([-l l]); zlim([-l l])
    ts.EdgeColor = 'none';
    colormap(bone)
    colorbar;
    clim([0 maxE])
    axis vis3d
    view(2)
    pause(0.08)
end

