function [ts, max_loc, whichTriangle, maxR2, corners] = plot_R2_map(in, only_roi, metric)
arguments
    in
    only_roi (1,1) {mustBeA(only_roi, "logical")} = true;
    metric = 'mag';
end

if only_roi
    infix = '_roi';
else
    infix = '';
end

if ~isempty(metric)
    metric = sprintf('_%s', metric);
end

geo_file = [in '/r2' infix '_geo.hdf5'];
data_file = [in '/r2' infix '_data.hdf5'];
coordinates = h5read(geo_file, '/mesh/nodes/node_coord');
triangles = h5read(geo_file, '/mesh/elm/triangle_number_list')' + 1;
data = h5read(data_file, sprintf('/data/tris/c_E%s', metric));
data(isnan(data)) = 0;

data = reject_R2_outliers_if_needed(data);

[maxR2, whichTriangle] = max(data);
corners = coordinates(:,triangles(whichTriangle,:));
max_loc = mean(corners, 2);

ts = trisurf(triangles, coordinates(1,:), coordinates(2,:), coordinates(3,:), data);
ts.EdgeColor = 'none';
end