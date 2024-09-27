ID = 'sub-014';
hemisphere = 'R';
roi = sprintf('midlayer_%s', lower(hemisphere));
suffix = '';
mesh_id = 'mesh0'; % sprintf('mesh_%s_%s', hemisphere, suffix);
ROOT = sprintf('//wsl.localhost/Ubuntu-22.04/home/bnplab-admin/TMS_localization/HighDef/%s/results/exp_map-%s', ID, hemisphere);
R2_PATH = sprintf('%s/r2/mesh_%s/roi_%s', ROOT, mesh_id, roi);
IN_APB = sprintf('%s/CsE_APB_in_uV/sigmoid4', R2_PATH);
IN_FDI = sprintf('%s/CsE_FDI_in_uV/sigmoid4', R2_PATH);
IN_ADM = sprintf('%s/CsE_ADM_in_uV/sigmoid4', R2_PATH);

responses = readtable(sprintf('//wsl.localhost/Ubuntu-22.04/home/bnplab-admin/TMS_localization/HighDef/%s/%s_%s_raw.csv', ID, ID, hemisphere));


l = 80;

metric = 'mag';

fig = figure(Position=[50 50 2000 600]);

ax1 = subplot(1,3,1);
[ts_APB, max_loc_APB, bestMatch_APB, maxR2_APB] = visualize(IN_APB, true, metric);
hold on
format_axis()
title('Hotspot APB')

ax2 = subplot(1,3,2);
[ts_FDI, max_loc_FDI, bestMatch_FDI, maxR2_FDI] = visualize(IN_FDI, true, metric);
hold on
format_axis()
title('Hotspot FDI')

ax3 = subplot(1,3,3);
[ts_ADM, max_loc_ADM, bestMatch_ADM, maxR2_ADM] = visualize(IN_ADM, true, metric);
hold on
format_axis()
title('Hotspot ADM')

fontsize(fig, scale=2)



plot3(ax1, max_loc_APB(1), max_loc_APB(2), max_loc_APB(3)+1, 'wx', MarkerSize=10, LineWidth=2)
plot3(ax1, max_loc_FDI(1), max_loc_FDI(2), max_loc_FDI(3)+1, 'w.', MarkerSize=10, LineWidth=3)
plot3(ax1, max_loc_ADM(1), max_loc_ADM(2), max_loc_ADM(3)+1, 'w.', MarkerSize=10, LineWidth=3)

plot3(ax2, max_loc_APB(1), max_loc_APB(2), max_loc_APB(3)+1, 'w.', MarkerSize=10, LineWidth=3)
plot3(ax2, max_loc_FDI(1), max_loc_FDI(2), max_loc_FDI(3)+1, 'wx', MarkerSize=10, LineWidth=2)
plot3(ax2, max_loc_ADM(1), max_loc_ADM(2), max_loc_ADM(3)+1, 'w.', MarkerSize=10, LineWidth=3)

plot3(ax3, max_loc_APB(1), max_loc_APB(2), max_loc_APB(3)+1, 'w.', MarkerSize=10, LineWidth=3)
plot3(ax3, max_loc_FDI(1), max_loc_FDI(2), max_loc_FDI(3)+1, 'w.', MarkerSize=10, LineWidth=3)
plot3(ax3, max_loc_ADM(1), max_loc_ADM(2), max_loc_ADM(3)+1, 'wx', MarkerSize=10, LineWidth=2)

exportgraphics(fig, sprintf('B:/Projects/2023-01 HighDef/Results/R2-Figures/%s_%s_%s_localization%s.png', ID, hemisphere, metric, suffix))

% Note spots in overview table
update_spots(ID, 'ini', hemisphere, 'APB', 'hot', max_loc_APB(1), max_loc_APB(2), max_loc_APB(3), maxR2_APB)
update_spots(ID, 'ini', hemisphere, 'FDI', 'hot', max_loc_FDI(1), max_loc_FDI(2), max_loc_FDI(3), maxR2_FDI)
update_spots(ID, 'ini', hemisphere, 'ADM', 'hot', max_loc_ADM(1), max_loc_ADM(2), max_loc_ADM(3), maxR2_ADM)



%%

simulation_result = sprintf('%s/electric_field/mesh_%s/roi_%s/e.hdf5', ROOT, mesh_id, roi);
E_mag = h5read(simulation_result, '/E_mag');
figure;
scatter(E_mag(bestMatch_FDI,:), responses.CsE_FDI_in_uV, 'k.')



function format_axis()
%t = -70:10:0;
xL = xlim; yL = ylim; zL = zlim;
rX = range(xL); rY = range(yL); rZ = range(zL);
r = max([rX rY rZ]);
xlim(mean(xL) + 0.5.*[-r r]); ylim(mean(yL) + 0.5.*[-r r]); zlim(mean(zL) + 0.5.*[-r r])
%xticks(t)
xlabel('\leftarrow Left [mm]')
ylabel('Front [mm] \rightarrow')
cb = colorbar; ylabel(cb, 'R^2', Rotation=0)
colormap(hotspotColormap())
set(gca(), 'Color', 'none')
view(2)
axis vis3d
end


function [ts, max_loc, whichTriangle, maxR2] = visualize(in, only_roi, metric)
if only_roi
    infix = '_roi';
else
    infix = '';
end

geo_file = [in '/r2' infix '_geo.hdf5'];
data_file = [in '/r2' infix '_data.hdf5'];
if isfile(geo_file) && isfile(data_file)
    coordinates = h5read(geo_file, '/mesh/nodes/node_coord');
    triangles = h5read(geo_file, '/mesh/elm/triangle_number_list')' + 1;
    data = h5read(data_file, sprintf('/data/tris/c_E_%s', metric));
    data(isnan(data)) = 0;
    [maxR2, whichTriangle] = max(data);
    max_loc = mean(coordinates(:,triangles(whichTriangle,:)), 2);
    
    ts = trisurf(triangles, coordinates(1,:), coordinates(2,:), coordinates(3,:), data);
    ts.EdgeColor = 'none';
else
    text(0.5,0.5, "Not acquired", HorizontalAlignment="center")
    max_loc = nan(1,3);
    ts = [];
    whichTriangle = nan;
    maxR2 = nan;
end
end

function [map] = hotspotColormap()
breaks = [0 0.2 0.4 0.75 0.875 1];
values = [-80 0.1 0; 
          -60 0.65 0.25; 
          -40 0.75 0.5; 
           0 0.7 0.97; 
           15 0.67 1;
           30 0.63 1];

hotspotColormapHsv = interp1(breaks, values, linspace(0,1,100));
hotspotColormapHsv(:,1) = mod(hotspotColormapHsv(:,1), 360) ./ 360;
map = hsv2rgb(hotspotColormapHsv);
end