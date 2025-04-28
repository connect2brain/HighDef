ID = "sub-003";
exp_id = "L2R";
hemisphere = extractBefore(exp_id, "2");
mesh_id = 'mesh0';

roi = sprintf('midlayer_%s', lower(hemisphere));
ROOT = sprintf('//wsl.localhost/Ubuntu-22.04/home/bnplab-admin/TMS_localization/HighDef/%s/results/exp_map-%s', ID, exp_id);
R2_PATH = sprintf('%s/r2/mesh_%s/roi_%s', ROOT, mesh_id, roi);

templates = [];
templates.hot = 'CsE_%s_in_uV';
templates.cold = 'SIHIscore_%s';

geo_file = sprintf('//wsl.localhost/Ubuntu-22.04/home/bnplab-admin/TMS_localization/HighDef/%s/mesh/roi/%s/geo.hdf5', ID, roi);
coordinates = h5read(geo_file, '/mesh/nodes/node_coord');
triangles = h5read(geo_file, '/mesh/elm/triangle_number_list')' + 1;

nTriangles = size(triangles,1);
centers = nan(3, nTriangles);
normals = nan(3, nTriangles);
for iTriangle = 1:nTriangles
    node_coordinates = coordinates(:, triangles(iTriangle,:));
    centers(:,iTriangle) = mean(node_coordinates, 2);
    normals(:,iTriangle) = cross(node_coordinates(:,2) - node_coordinates(:,1), node_coordinates(:,3) - node_coordinates(:,1));
    normals(:,iTriangle) = normals(:,iTriangle) ./ norm(normals(:,iTriangle)); % unit length
end


raw_excentricity = sqrt(sum((coordinates - [0; 0; 10]) .^ 2, 1));

%%
fig = figure(Position=[50 50 600 600]);
hold on
ts = trisurf(triangles, coordinates(1,:), coordinates(2,:), coordinates(3,:), mean(raw_excentricity(triangles), 2));
ts.EdgeColor="none";
colormap("turbo")
hold on

nStreamlines = 50;
nSteps = 30;
constantforce = [-1;0;1];
constantforce = 0.67 .* constantforce ./ norm(constantforce);
for iStreamline = 1:nStreamlines
    p = nan(3,nSteps);
    p(:,1) = centers(:, randi(nTriangles));
    for iStep = 2:nSteps
        weights = exp(-vecnorm(centers - p(:, iStep-1)));
        force = sum(weights .* normals, 2) / sum(weights);
        p(:, iStep) = p(:, iStep-1) + force + constantforce;
    end
    plot3(p(1,:), p(2,:), p(3,:))
end


%plot3(centers(1,:), centers(2,:), centers(3,:), 'wo', MarkerSize=3)
%plot3(centers(1,:) + [zeros(1,nTriangles); normals(1,:)], centers(2,:) + [zeros(1,nTriangles); normals(2,:)], centers(3,:) + [zeros(1,nTriangles); normals(3,:)], 'w')

axis square; box on;
set(gca, Color="none");
view(2)

