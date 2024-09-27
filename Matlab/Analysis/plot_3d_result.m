ID = 'sub-014';
exp_id = 'R2L';
hemisphere = exp_id(1);
%mesh_id = sprintf('mesh%s_ur', hemisphere);
mesh_id = 'mesh0';

muscle = 'FDI';

roi = sprintf('midlayer_%s', lower(hemisphere));
%roi = sprintf('small_%s', lower(hemisphere));
ROOT = sprintf('//wsl.localhost/Ubuntu-22.04/home/bnplab-admin/TMS_localization/HighDef/%s/results/exp_map-%s', ID, exp_id);
R2_PATH = sprintf('%s/r2/mesh_%s/roi_%s', ROOT, mesh_id, roi);
simulation_result = sprintf('%s/electric_field/mesh_%s/roi_%s/e.hdf5', ROOT, mesh_id, roi);

excitationName = sprintf('CsE_%s_in_uV', muscle);
inhibitionName = sprintf('SIHIscore_%s', muscle);
CR_name = sprintf('inhibited_mep_%s_in_uV', muscle);


IN_X = sprintf('%s/%s/sigmoid4', R2_PATH, excitationName);
IN_I = sprintf('%s/%s/sigmoid4', R2_PATH, inhibitionName);

if strcmpi(exp_id, 'L2R-plus-pi')
    responses = readtable(sprintf('//wsl.localhost/Ubuntu-22.04/home/bnplab-admin/TMS_localization/HighDef/%s/%s_map-%s_raw.csv', ID, ID, exp_id));
else
    responses = readtable(sprintf('//wsl.localhost/Ubuntu-22.04/home/bnplab-admin/TMS_localization/HighDef/%s/%s_%s_raw.csv', ID, ID, exp_id));
end
%



metric = 'mag';

% x -60 -10   /\=50 -> -70, 0
% y -50  20   /\=70 -> -50, 20
% z  10  50   /\=40 -> -5, 65



fig = figure(Position=[50 50 1600 600]);

ax1 = subplot(1,2,1);
[ts_excitability, max_loc_I, bestMatch_I, maxR2_I] = plot_R2_map(IN_I, true, metric);
hold on
format_axis()
title(sprintf('Coldspot localization: R² = %0.3f', maxR2_I))


ax2 = subplot(1,2,2);
[ts_inhibition, max_loc_X, bestMatch_X, maxR2_X] = plot_R2_map(IN_X, true, metric);
hold on
format_axis()
title(sprintf('Hotspot localization: R² = %0.3f', maxR2_X))

fontsize(fig, scale=2)
fprintf('Best R² (hot)  = %0.5f\nBest R² (cold) = %0.5f\n\n', maxR2_X, maxR2_I)


plot3(ax1, max_loc_X(1), max_loc_X(2), max_loc_X(3)+1, 'w^', MarkerSize=10, LineWidth=1)
plot3(ax1, max_loc_I(1), max_loc_I(2), max_loc_I(3)+1, 'wv', MarkerSize=10, LineWidth=3)

plot3(ax2, max_loc_X(1), max_loc_X(2), max_loc_X(3)+1, 'w^', MarkerSize=10, LineWidth=3)
plot3(ax2, max_loc_I(1), max_loc_I(2), max_loc_I(3)+1, 'wv', MarkerSize=10, LineWidth=1)

exportgraphics(fig, sprintf('B:/Projects/2023-01 HighDef/Results/R2-Figures/%s_%s_%s_localization_%s.png', ID, exp_id, metric, muscle))
fprintf(' Distance between coldspot and hotspot: %.3f mm\n\n', sqrt(sum((max_loc_I - max_loc_X).^2)))


% Note spots in overview table
update_spots(ID, exp_id, hemisphere, muscle, 'hot', max_loc_X(1), max_loc_X(2), max_loc_X(3), maxR2_X)
update_spots(ID, exp_id, hemisphere, muscle, 'cold', max_loc_I(1), max_loc_I(2), max_loc_I(3), maxR2_I)





%%
E_mag = h5read(simulation_result, '/E_mag');
% We put dI/dt in as A/µs; thus, we get V/µm (1e-6 V/m)
unit_correction = 1e6;
usedIntensities = unique(responses.Intensity_percentMSO);
intensityGroup = responses.Intensity_percentMSO ~= max(usedIntensities);

figio = figure(Position=[150 400 1200 500]);

CR_too_low = responses.(CR_name) < 40;

%figio = figure(Position=[150 0 1200 1200]);
%subplot(2,2,1)
subplot(1,2,1)
scatter_and_movmean(unit_correction.*E_mag(bestMatch_I,~CR_too_low),responses.(inhibitionName)(~CR_too_low), intensityGroup(~CR_too_low))
title('I/O-curve at the triangle best matching inhibition', FontSize=8)
ylabel('SIHIscore')
xlabel('|E| [V/m]')
set(gca, 'Color', 'none')

% subplot(2,2,2)
% scatter_and_movmean(unit_correction.*E_mag(bestMatch_I,:),responses.(excitationName), intensityGroup)
% title('I/O-curve at the triangle best matching inhibition', FontSize=8)
% ylabel('MEP amplitude [µV]')
% xlabel('|E| [V/m]')
% set(gca, 'Color', 'none')

% subplot(2,2,3)
% scatter_and_movmean(unit_correction.*E_mag(bestMatch_X,:),responses.(inhibitionName), intensityGroup)
% title('I/O-curve at the triangle best matching excitation', FontSize=8)
% ylabel('SIHIscore')
% xlabel('|E| [V/m]')
% set(gca, 'Color', 'none')

% subplot(2,2,4)
subplot(1,2,2)
scatter_and_movmean(unit_correction.*E_mag(bestMatch_X,:), responses.(excitationName), intensityGroup)
title('I/O-curve at the triangle best matching excitation', FontSize=8)
ylabel('MEP amplitude [µV]')
xlabel('|E| [V/m]')
set(gca, 'Color', 'none')



fontsize(figio, scale=1.5)
exportgraphics(figio, sprintf('B:/Projects/2023-01 HighDef/Results/R2-Figures/%s_%s_%s_bestmatch_IO_%s.png', ID, exp_id, metric, muscle))




%figure;
%scatter_and_movmean(unit_correction.*E_mag(bestMatch_I,:),responses.inhibited_mep_in_uV, ~intensityGroup)


function scatter_and_movmean(x,y, grouping)
xmin = min(x);
xmax = max(x);
nWindows = 100;
smoothness = range(x) / 10;
xq = linspace(xmin, xmax, nWindows);
yq = nan(1, nWindows);
vq = nan(1, nWindows);
for iWindow = 1:nWindows
   center = xq(iWindow);
   weights = exp(-0.5.*((x-center)./smoothness).^2)./(smoothness*sqrt(2*pi));
   denom = sum(weights);
   yq(iWindow) = sum(weights'.*y) / denom;
   vq(iWindow) = sum(weights' .* (y - yq(iWindow)).^2) / denom;
end
fill([xq fliplr(xq)], [yq - sqrt(vq), fliplr(yq + sqrt(vq))], 'k', EdgeColor='none', FaceAlpha=0.1)
hold on
plot(xq, yq, 'k-')
scatter(x(~grouping),y(~grouping),'kx')
scatter(x(grouping),y(grouping),'b.')
end




function format_axis()
%t = -70:10:0;
limit_range = max([range(xlim) range(ylim) range(zlim)]);
window = [-limit_range limit_range] ./ 2;

xlim(mean(xlim) + window); ylim(mean(ylim) + window); zlim(mean(zlim) + window)
%xticks(t)
xlabel('\leftarrow Left [mm]')
ylabel('Front [mm] \rightarrow')
cb = colorbar; ylabel(cb, 'R^2', Rotation=0)
colormap(hotspotColormap())
set(gca(), 'Color', 'none')
view(2)
axis square
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