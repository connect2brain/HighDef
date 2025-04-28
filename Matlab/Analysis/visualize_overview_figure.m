addpath(genpath('B:/Projects/2023-01 HighDef/libraries/visualization'))

ID = 'sub-003';
exp_id = 'L2R';
hemisphere = exp_id(1);
mesh_id = 'mesh0';

%%
% Get some example MEPs


addpath('B:\Projects\2020-12 VETTERPHD Project\libraries\neurone_tools_for_matlab_1.1.3.11_mod')
addpath(genpath('B:/Projects/2023-01 HighDef/libraries/vetter'))


ROOT = 'B:/Experimental Data/2023-01 HighDef/Conventional';
tableFile = sprintf('%s/Sessions.xlsx', ROOT);
T = readtable(tableFile);
subjects = unique(T.Subject);
fprintf('Listed subjects: %s\n\n', sprintf(' %s', subjects{:}))

IN = sprintf('%s/%s', ROOT, ID);
ROWS = T(strcmpi(T.Subject, ID),:);
ROWS = ROWS(startsWith(ROWS.Session, 'map') & endsWith(ROWS.Session, exp_id) & ismember(ROWS.Condition, {'High', 'Low'}) & ~endsWith(ROWS.Session, 'old'),:);

row = ROWS(1,:);
[neurone_times, US, CS, TS, mepsFDI, mepsAPB, mepsADM] = read_neurone(ROOT, ID, sprintf("map-%s", exp_id), row.NeurOneFolder{:}, row.NeurOneIndex);

%%

selectedTrial = 47; % 1, 14, 42, 43, 45, 47
selectedCS = find(CS); selectedCS = selectedCS(selectedTrial);
selectedTS = find(TS); selectedTS = selectedTS(selectedTrial);

fig = figure(Position=[50 50 700 150]);
tiledlayout(1,4)

allaxes = [];
allaxes(1) = nexttile(2);
plot(detrend(mepsADM(:,selectedCS), 0), 'k-')
axis square; axis off


allaxes(2) = nexttile(3);
plot(detrend(mepsFDI(:,selectedCS), 0), 'k-')
axis square; axis off


allaxes(3) = nexttile(4);
plot(detrend(mepsAPB(:,selectedCS), 0), 'k-')
axis square; axis off


allaxes(4) = nexttile(1);
meanUSFDI = mean(detrend(mepsFDI(:,US), 0), 2);
plot(meanUSFDI, 'k:')
hold on
plot(detrend(mepsFDI(:,selectedTS), 0), 'k-')
axis square; axis off


linkaxes(allaxes, "y")
l = max(abs(ylim));
ylim([-l l])



exportgraphics(fig, sprintf('B:/Projects/2023-01 HighDef/Results/Evaluation/illustration_meps.pdf'))
exportgraphics(fig, sprintf('B:/Projects/2023-01 HighDef/Results/Evaluation/illustration_meps.png'))

%%


roi = sprintf('midlayer_%s', lower(hemisphere));
ROOT = sprintf('D:/HighDef-operate/HighDef/%s/results/exp_map-%s', ID, exp_id);
R2_PATH = sprintf('%s/r2/mesh_%s/roi_%s', ROOT, mesh_id, roi);

templates = [];
templates.hot = 'CsE_%s_in_uV';
templates.cold = 'SIHIscore_%s';

tileIndices = [];
% tileIndices.map.ADM.hot = 3;    tileIndices.io.ADM.hot = 4;
% tileIndices.map.FDI.hot = 7;    tileIndices.io.FDI.hot = 8;
% tileIndices.map.APB.hot = 11;   tileIndices.io.APB.hot = 12;
% tileIndices.map.FDI.cold = 6;   tileIndices.io.FDI.cold = 5;
tileIndices.map.ADM.hot  = 8;   tileIndices.io.ADM.hot  = 12;
tileIndices.map.FDI.hot  = 6;   tileIndices.io.FDI.hot  = 10;
tileIndices.map.APB.hot  = 7;   tileIndices.io.APB.hot  = 11;
tileIndices.map.FDI.cold = 5;   tileIndices.io.FDI.cold = 9;


fig = figure(Position=[50 500 675 385]);
%layout = tiledlayout(3, 4, TileSpacing="tight");
layout = tiledlayout(3, 4, TileSpacing="tight");

allaxes = [];
allaxes(1) = nexttile(4);
plot(detrend(mepsADM(:,selectedCS), 0), 'k-')
axis square; axis off


allaxes(2) = nexttile(2);
plot(detrend(mepsFDI(:,selectedCS), 0), 'k-')
axis square; axis off


allaxes(3) = nexttile(3);
plot(detrend(mepsAPB(:,selectedCS), 0), 'k-')
axis square; axis off


allaxes(4) = nexttile(1);
meanUSFDI = mean(detrend(mepsFDI(:,US), 0), 2);
plot(meanUSFDI, 'k:')
hold on
plot(detrend(mepsFDI(:,selectedTS), 0), 'k-')
axis square; axis off


linkaxes(allaxes, "y")
l = max(abs(ylim));
ylim([-l l])



for muscle = ["ADM", "FDI", "APB"]
    if muscle == "FDI"
        spots = ["hot", "cold"];
    else
        spots = ["hot"];
    end
    for spot = spots
        % Map:
        nexttile(tileIndices.map.(muscle).(spot))
        [ts, max_loc, bestMatch, maxR2, corners] = plot_R2_map(sprintf('%s/%s/sigmoid4', R2_PATH, sprintf(templates.(spot), muscle)), true, 'mag');
        hold on
        plot3(max_loc(1), max_loc(2), max_loc(3)+5, "ko")
        colormap(warm_colormap())
        view(2)
        box off; axis off; axis square
        limrange = max([range(xlim) range(ylim) range(zlim)]);
        window = [-1 1] .* limrange ./ 2;
        xlim(mean(xlim) + window); ylim(mean(ylim) + window); zlim(mean(zlim) + window)

        % I/O-curve
        E_mag = h5read(sprintf('%s/electric_field/mesh_%s/roi_%s/e.hdf5', ROOT, mesh_id, roi), '/E_mag');
        responses = readtable(sprintf('D:/HighDef-operate/HighDef/%s/%s_%s_raw.csv', ID, ID, exp_id));
        sigmoid4 = 'a + (b - a)/(1+exp(-c*(x - d)))';
        x = E_mag(bestMatch,:)' .* 1e6; % bc. the SI is given in A/us (i.e. 1e6 A/s), the E-field magnitude is given in uV/m (i.e. 1e-6 V/m)
        y = responses.(sprintf(templates.(spot), muscle));
        k0 = range(y)/range(x);
        x0 = mean(x);
        % This stuff is only for visualization! The actual regression
        % follows some complicated code by Weise, Numssen et al.
        if spot == "hot"
            if muscle == "FDI"
                k0 = 0.1e-5 * k0;
                x0 = 1.8*x0;
                bestMatchFDIhot = bestMatch;
            else
                k0 = 1e-2 * k0;
            end
        end
        [f1, gof, ~] = fit(x,y,sigmoid4,'Start', [quantile(y, 0.25) quantile(y, 0.75) k0 x0]);
        
        nexttile(tileIndices.io.(muscle).(spot))
        plot(x, y, 'k.', MarkerSize=1)
        hold on
        if spot == "hot"
            fitplot = plot(f1, "r");
        else
            fitplot = plot(f1, "b");
        end
        fitplot.LineWidth = 2;
        legend off
        axis square
        box on
        set(gca, Color="none");
        xlabel("|E| [V/m]")
        if spot == "hot"
            ylabel("MEP amplitude [\muV]")
            %set(gca, YAxisLocation="right")
        else
            ylabel("SIHI score")
        end
        xticks([0 50])

    end
end



exportgraphics(fig, sprintf('B:/Projects/2023-01 HighDef/Results/Evaluation/illustration_all_spots_and_io.pdf'), Resolution=450)
exportgraphics(fig, sprintf('B:/Projects/2023-01 HighDef/Results/Evaluation/illustration_all_spots_and_io.png'))

%%

fig_IO_only_poster = figure(Position=[50 150 136 122]);
muscle = "FDI";
spot = "cold";
cerulean = [0, 150, 255]./255;

% I/O-curve
E_mag = h5read(sprintf('%s/electric_field/mesh_%s/roi_%s/e.hdf5', ROOT, mesh_id, roi), '/E_mag');
responses = readtable(sprintf('D:/HighDef-operate/HighDef/%s/%s_%s_raw.csv', ID, ID, exp_id));
sigmoid4 = 'a + (b - a)/(1+exp(-c*(x - d)))';
x = E_mag(bestMatch,:)' .* 1e6; % bc. the SI is given in A/us (i.e. 1e6 A/s), the E-field magnitude is given in uV/m (i.e. 1e-6 V/m)
y = responses.(sprintf(templates.(spot), muscle));
k0 = range(y)/range(x);
x0 = mean(x);
k0 = 1e-2 * k0;
[f1, gof, ~] = fit(x,y,sigmoid4,'Start', [quantile(y, 0.25) quantile(y, 0.75) k0 x0]);

nexttile(tileIndices.io.(muscle).(spot))
plot(x, y, 'k.', MarkerSize=1)
hold on
fitplot = plot(f1);
fitplot.LineWidth = 2;
fitplot.Color = cerulean;
legend off
axis square
box on
set(gca, Color="none");
xlabel("|E| [V/m]")
ylabel("SIHI score")
xticks([0 50])


exportgraphics(fig_IO_only_poster, sprintf('B:/Projects/2023-01 HighDef/Results/Evaluation/illustration_mini_io.pdf'), BackgroundColor="none")

%% Visualize the surface per se, and then with some response variable (FDI MEP)
otherTriangle = 10500;

fig = figure;
geo_file = sprintf('D:/HighDef-operate/HighDef/%s/mesh/roi/%s/geo.hdf5', ID, roi);
coordinates = h5read(geo_file, '/mesh/nodes/node_coord');
triangles = h5read(geo_file, '/mesh/elm/triangle_number_list')' + 1;
ts = trisurf(triangles, coordinates(1,:), coordinates(2,:), coordinates(3,:)); %, 1:size(triangles,1));
ts.EdgeColor = 'none';
ts.FaceAlpha = 1;
colormap("gray")
%colorbar
hold on
s2 = trisurf(triangles([otherTriangle bestMatchFDIhot],:), coordinates(1,:), coordinates(2,:), coordinates(3,:));
s2.FaceColor="none";
view(2)
box off; axis off; axis square
limrange = max([range(xlim) range(ylim) range(zlim)]);
window = [-1 1] .* limrange ./ 2;
xlim(mean(xlim) + window); ylim(mean(ylim) + window); zlim(mean(zlim) + window)
exportgraphics(fig, sprintf('B:/Projects/2023-01 HighDef/Results/Evaluation/empty_roi.pdf'), Resolution=450)
exportgraphics(fig, sprintf('B:/Projects/2023-01 HighDef/Results/Evaluation/empty_roi.png'))

%%

fig = figure;
[ts, max_loc, bestMatch, maxR2, corners] = plot_R2_map(sprintf('%s/%s/sigmoid4', R2_PATH, sprintf(templates.hot, "FDI")), true, 'mag');
colormap(warm_colormap())
hold on
s2 = trisurf(triangles([otherTriangle bestMatchFDIhot],:), coordinates(1,:), coordinates(2,:), coordinates(3,:), 0);
s2.FaceColor="none";
view(2)
box off; axis off; axis square
limrange = max([range(xlim) range(ylim) range(zlim)]);
window = [-1 1] .* limrange ./ 2;
xlim(mean(xlim) + window); ylim(mean(ylim) + window); zlim(mean(zlim) + window)
cb = colorbar(XTick=[0 0.2 0.4 0.6]);
a = get(cb,'YTickLabel');
set(cb,'YTickLabel',a,'fontsize',20)
ylabel(cb, "R²", rotation=0, FontSize=20)

exportgraphics(fig, sprintf('B:/Projects/2023-01 HighDef/Results/Evaluation/colored_example_map.pdf'), Resolution=450)
exportgraphics(fig, sprintf('B:/Projects/2023-01 HighDef/Results/Evaluation/colored_example_map.png'))
%%
fig = figure(Position=[50 50 510 260]);
tiledlayout(2,2);
% I/O-curve
E_mag = h5read(sprintf('%s/electric_field/mesh_%s/roi_%s/e.hdf5', ROOT, mesh_id, roi), '/E_mag');
responses = readtable(sprintf('D:/HighDef-operate/HighDef/%s/%s_%s_raw.csv', ID, ID, exp_id));
sigmoid4 = 'a + (b - a)/(1+exp(-c*(x - d)))';



x = E_mag(bestMatch,:)' .* 1e6; % bc. the SI is given in A/us (i.e. 1e6 A/s), the E-field magnitude is given in uV/m (i.e. 1e-6 V/m)
y = responses.(sprintf(templates.(spot), muscle));
k0 = range(y)/range(x);
x0 = mean(x);
% This stuff is only for visualization! The actual regression
% follows some complicated code by Weise, Numssen et al.
%k0 = 0.1e-5 * k0;
k0 = 1e-6 * k0;
x0 = 1.8*x0;

[f1, gof, ~] = fit(x,y,sigmoid4,'Start', [quantile(y, 0.25) quantile(y, 0.75) k0 x0]);
fprintf("R² = %.4f\n", gof.rsquare)

nexttile(3)
plot(x, y, 'k.', MarkerSize=1)
axis square
box on
set(gca, Color="none");
xlabel("|E| [V/m]")
ylabel("Response")
xticks([0 50])
yticks([])
ylim([0 quantile(y, 0.99)])

nexttile(4)
plot(x, y, 'k.', MarkerSize=1)
hold on
plot(f1, "r")
legend off
axis square
box on
set(gca, Color="none");
xlabel("|E| [V/m]")
ylabel("Response")
xticks([0 50])
yticks([])
ylim([0 quantile(y, 0.99)])
text(max(xlim), f1(max(xlim)), sprintf("R² = %.4f", gof.rsquare), HorizontalAlignment="left", Color="r")

x = E_mag(otherTriangle,:)' .* 1e6; % bc. the SI is given in A/us (i.e. 1e6 A/s), the E-field magnitude is given in uV/m (i.e. 1e-6 V/m)
y = responses.(sprintf(templates.(spot), muscle));
k0 = range(y)/range(x);
x0 = mean(x);
% This stuff is only for visualization! The actual regression
% follows some complicated code by Weise, Numssen et al.
%k0 = 0.1e-5 * k0;
k0 = 1e-6 * k0;
x0 = 1.8*x0;

[f1, gof, ~] = fit(x,y,sigmoid4,'Start', [quantile(y, 0.25) quantile(y, 0.75) k0 x0]);
fprintf("R² = %.4f\n", gof.rsquare)

nexttile(1)
plot(x, y, 'k.', MarkerSize=1)
axis square
box on
set(gca, Color="none");
xlabel("|E| [V/m]")
ylabel("Response")
xticks([0 50])
yticks([])
ylim([0 quantile(y, 0.99)])

nexttile(2)
plot(x, y, 'k.', MarkerSize=1)
hold on
plot(f1, "r")
legend off
axis square
box on
set(gca, Color="none");
xlabel("|E| [V/m]")
ylabel("Response")
xticks([0 50])
yticks([])
ylim([0 quantile(y, 0.99)])
text(max(xlim), f1(max(xlim)), sprintf("R² = %.4f", gof.rsquare), HorizontalAlignment="left", Color="r")


exportgraphics(fig, sprintf('B:/Projects/2023-01 HighDef/Results/Evaluation/fits.pdf'), Resolution=450)
exportgraphics(fig, sprintf('B:/Projects/2023-01 HighDef/Results/Evaluation/fits.png'), Resolution=900)














function [time, US, CS, TS, mepFDI, mepAPB, mepADM] = read_neurone(basepath, subject, session, folder, index)
% Parsing the session name into the muscles needed
hemisphereCode = extractAfter(session, '-');
% If we are in session map-R2L, then the US and TS are delivered to the
% Left M1, and the corresponding MEP is in the RIGHT HAND!
us_handside = lower(extractBefore(hemisphereCode, '2'));
ts_handside = us_handside;
cs_handside = lower(extractAfter(hemisphereCode, '2'));

fprintf('    Processing %s as: CS @ %s   TS @ %s   US @ %s\n', session, cs_handside, ts_handside, us_handside)


subjData = module_read_neurone(sprintf("%s/%s/%s/%s", basepath, subject, session, folder), sessionPhaseNumber=index);
markerLabels = subjData.markers.type;

if ismember('144', unique(markerLabels))
    markerLabels = markerLabels(~strcmpi(markerLabels, '144'));
    markerLabels(~isnan(str2double(markerLabels))) = num2cell(cellfun(@(s) num2str(str2double(s) - 144), markerLabels(~isnan(str2double(markerLabels)))));
end

codeSingleRight = '4';
codePairedRight = '1';
codePairedLeft  = '8';

outMarkers = find(strcmpi(markerLabels, 'Out'));
truncatedOutMarkerNext = arrayfun(@(x) min([x length(markerLabels)]), outMarkers+1);

CS = (strcmpi(markerLabels(outMarkers-1), codePairedLeft) | strcmpi(markerLabels(truncatedOutMarkerNext), codePairedLeft))';
TS = (strcmpi(markerLabels(outMarkers-1), codePairedRight) | strcmpi(markerLabels(truncatedOutMarkerNext), codePairedRight))';
US = (strcmpi(markerLabels(outMarkers-1), codeSingleRight) | strcmpi(markerLabels(truncatedOutMarkerNext), codeSingleRight))';
triggerIdcs = subjData.markers.index(strcmpi(subjData.markers.type, 'Out')); % This is NOT the same as outMarkers! markerLabels that leads to outMarkers is edited.
time = subjData.markers.time(strcmpi(subjData.markers.type, 'Out'))';

MEPwindow_in_s = [0.02, 0.045];
MEPwindow = MEPwindow_in_s .* subjData.properties.samplingRate;
MEPwindow = round(MEPwindow(1)):round(MEPwindow(2));
MEPwindow_cs = MEPwindow; % DO NOT! + subjData.properties.samplingRate*0.010; % Shift by 10 ms
MEPwindow_shifted_us = MEPwindow - subjData.properties.samplingRate*0.010; % DO NOT! + subjData.properties.samplingRate*0.010; % Shift by 10 ms

mepFDI = nan(length(MEPwindow), length(outMarkers));
mepFDI(:,US) = subjData.signal.(sprintf('FDI%s', us_handside)).data(triggerIdcs(US) + MEPwindow_shifted_us)';
mepFDI(:,CS) = subjData.signal.(sprintf('FDI%s', cs_handside)).data(triggerIdcs(CS) + MEPwindow_cs)';
mepFDI(:,TS) = subjData.signal.(sprintf('FDI%s', ts_handside)).data(triggerIdcs(CS) + MEPwindow)';

mepAPB = nan(length(MEPwindow), length(outMarkers));
mepAPB(:,US) = subjData.signal.(sprintf('APB%s', us_handside)).data(triggerIdcs(US) + MEPwindow_shifted_us)';
mepAPB(:,CS) = subjData.signal.(sprintf('APB%s', cs_handside)).data(triggerIdcs(CS) + MEPwindow_cs)';
mepAPB(:,TS) = subjData.signal.(sprintf('APB%s', ts_handside)).data(triggerIdcs(CS) + MEPwindow)';

mepADM = nan(length(MEPwindow), length(outMarkers));
mepADM(:,US) = subjData.signal.(sprintf('ADM%s', us_handside)).data(triggerIdcs(US) + MEPwindow_shifted_us)';
mepADM(:,CS) = subjData.signal.(sprintf('ADM%s', cs_handside)).data(triggerIdcs(CS) + MEPwindow_cs)';
mepADM(:,TS) = subjData.signal.(sprintf('ADM%s', ts_handside)).data(triggerIdcs(CS) + MEPwindow)';
end

