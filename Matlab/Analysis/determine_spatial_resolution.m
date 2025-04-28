ID = 'sub-014';
exp_id = 'R2L';
hemisphere = exp_id(1);
mesh_id = 'mesh0';
muscle = 'FDI';
roi = sprintf('midlayer_%s', lower(hemisphere));
ROOT = sprintf('//wsl.localhost/Ubuntu-22.04/home/bnplab-admin/TMS_localization/HighDef/%s/results/exp_map-%s', ID, exp_id);
R2_PATH = sprintf('%s/r2/mesh_%s/roi_%s', ROOT, mesh_id, roi);
simulation_result = sprintf('%s/electric_field/mesh_%s/roi_%s/e.hdf5', ROOT, mesh_id, roi);

excitationName = sprintf('CsE_%s_in_uV', muscle);
inhibitionName = sprintf('SIHIscore_%s', muscle);
CR_name = sprintf('inhibited_mep_%s_in_uV', muscle);


IN_X = sprintf('%s/%s/sigmoid4', R2_PATH, excitationName);
IN_I = sprintf('%s/%s/sigmoid4', R2_PATH, inhibitionName);

infix = '_roi';
geo_file = [IN_I '/r2' infix '_geo.hdf5'];
coordinates = h5read(geo_file, '/mesh/nodes/node_coord'); % in mm
triangles = h5read(geo_file, '/mesh/elm/triangle_number_list')' + 1;

nTriangles = size(triangles, 1);
centers = nan(3, nTriangles);
areas = nan(nTriangles, 1);
for iTriangle = 1:nTriangles
    node_coordinates = coordinates(:, triangles(iTriangle,:));
    centers(:,iTriangle) = mean(node_coordinates, 2);
    areas(iTriangle) = norm(cross(node_coordinates(:,2) - node_coordinates(:,1), node_coordinates(:,3) - node_coordinates(:,1))) / 2; % cross product's magnitude is area of parallelogram with vectors for sides -- need triangle: half that
    % Unit of area: mm² ( = 10^-6 m²)
end

%%
% Determine the area of each triangle
% Determine the smallest distance to any neighbor for each triangle
minDistances = nan(nTriangles, 1);
for iTriangle = 1:nTriangles
    minDistances(iTriangle) = min(vecnorm(centers(:,setdiff(1:nTriangles, iTriangle)) - centers(:,iTriangle)));
end

%%
n_vertices_in_surface = length(unique(triangles(:)));
fprintf("Average vertex density = %0.5f / mm²\n", n_vertices_in_surface / sum(areas))

%%
tiledlayout(2,1);
nexttile()
histogram(areas);
xlabel("Area [mm²]")
title("ROI triangle areas")

nexttile()
histogram(minDistances)
xlabel("Distance [mm]")
title("Distance from center of one triangle to nearest center of another triangle")