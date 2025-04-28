muscles = ["ADM", "FDI", "APB"];
muscleColor = [];
muscleColor.ADM = [0 0.5 1];
muscleColor.FDI = [1 0 0];
muscleColor.APB = [0.2 0 1];
sessionNames = []; sessionNames.L = ["map-L", "map-L2R"]; sessionNames.R = ["map-R", "map-R2L"];

addpath(genpath('B:/Projects/2023-01 HighDef/libraries/visualization'))

basepath = '//wsl.localhost/Ubuntu-22.04/home/bnplab-admin/TMS_localization/HighDef';
subjects = arrayfun(@(s) string(s.name), dir(fullfile(basepath, 'sub-*')))';

%% R E T R I E V E  D A T A
% For each subject:
% Subject, from_session, from_muscle, triangle_index, to_session, to_muscle, R2

% For all of the following 24 combinations:
% S1 ADM S1 ADM > QOF
% S1 ADM S2 ADM - Reliability
% S1 ADM S1 FDI : Specificity
% S1 ADM S1 APB : Specificity
% S1 FDI S1 FDI > QOF
% S1 FDI S2 FDI - Reliability
% S1 FDI S1 ADM : Specificity
% S1 FDI S1 APB : Specificity
% S1 APB S1 APB > QOF
% S1 APB S2 APB - Reliability
% S1 APB S1 ADM : Specificity
% S1 APB S1 FDI : Specificity
% S2 ADM S2 ADM > QOF
% S2 ADM S1 ADM - Reliability
% S2 ADM S2 FDI : Specificity
% S2 ADM S2 APB : Specificity
% S2 FDI S2 FDI > QOF
% S2 FDI S1 FDI - Reliability
% S2 FDI S2 ADM : Specificity
% S2 FDI S2 APB : Specificity
% S2 APB S2 APB > QOF
% S2 APB S1 APB - Reliability
% S2 APB S2 ADM : Specificity
% S2 APB S2 FDI : Specificity
%
% i.e. eventually: nSubjects * 2 * 24 = nSubjects * 48 rows

ResultTable = table([], [], [], [], [], [], [], [], [], [], VariableNames={'Subject', 'Hemisphere', 'from_spot', 'to_spot', 'from_session', 'from_muscle', 'triangle_index', 'to_session', 'to_muscle', 'R2'});

for subject = subjects
    for hemisphere = ["L", "R"]
        for iSession = 1:2
            exp_id = sessionNames.(hemisphere)(iSession);
            for iMuscle = 1:3
                jSession = 3-iSession; % 1->2; 2->1

                fn_iSession = sprintf('%s/%s/results/exp_%s/r2/mesh_mesh0/roi_midlayer_%s/CsE_%s_in_uV/sigmoid4/r2_roi_data.hdf5', ...
                    basepath, subject, exp_id, lower(hemisphere), muscles(iMuscle));
                if ~isfile(fn_iSession)
                    fprintf(' did not find %s\n', fn_iSession)
                    continue
                end
                fn_jSession = sprintf('%s/%s/results/exp_%s/r2/mesh_mesh0/roi_midlayer_%s/CsE_%s_in_uV/sigmoid4/r2_roi_data.hdf5', ...
                    basepath, subject, sessionNames.(hemisphere)(jSession), lower(hemisphere), muscles(iMuscle));
                if ~isfile(fn_jSession)
                    fprintf(' did not find %s\n', fn_jSession)
                    continue
                end

                % Quantify QOF
                data = reject_R2_outliers_if_needed(h5read(fn_iSession, '/data/tris/c_E_mag'));
                data(isnan(data)) = 0;
                [selfR2, bestIndex] = max(data);
                fprintf(' %s/%s  >  S%d %s S%d %s\n', subject, hemisphere, iSession, muscles(iMuscle), iSession, muscles(iMuscle))
                ResultTable = [ResultTable; {subject, hemisphere, 'hot', 'hot', sprintf("S%d", iSession), muscles(iMuscle), bestIndex, sprintf("S%d", iSession), muscles(iMuscle), selfR2}];

                % Quantify reliability
                data = reject_R2_outliers_if_needed(h5read(fn_jSession, '/data/tris/c_E_mag'));
                data(isnan(data)) = 0;
                crossR2 = data(bestIndex);
                maxR2 = max(data);
                fprintf(' %s/%s  -  S%d %s S%d %s\n', subject, hemisphere, iSession, muscles(iMuscle), jSession, muscles(iMuscle))
                ResultTable = [ResultTable; {subject, hemisphere, 'hot', 'hot', sprintf("S%d", iSession), muscles(iMuscle), bestIndex, sprintf("S%d", jSession), muscles(iMuscle), crossR2/maxR2}];

                % Quantify specificity
                for jMuscle = setdiff(1:3, iMuscle)
                    % Get triangle for iSession, iMuscle (=bestIndex), and retrieve R2 from
                    % iSession, jMuscle for same triangle
                    fn_jMuscle = sprintf('%s/%s/results/exp_%s/r2/mesh_mesh0/roi_midlayer_%s/CsE_%s_in_uV/sigmoid4/r2_roi_data.hdf5', ...
                        basepath, subject, sessionNames.(hemisphere)(iSession), lower(hemisphere), muscles(jMuscle));
                    if ~isfile(fn_jMuscle)
                        fprintf(' did not find %s\n', fn_jMuscle)
                        continue
                    end

                    data = reject_R2_outliers_if_needed(h5read(fn_jMuscle, '/data/tris/c_E_mag'));
                    data(isnan(data)) = 0;
                    jMuscleR2 = data(bestIndex);
                    fprintf(' %s/%s  :  S%d %s S%d %s\n', subject, hemisphere, iSession, muscles(iMuscle), iSession, muscles(jMuscle))
                    ResultTable = [ResultTable; {subject, hemisphere, 'hot', 'hot', sprintf("S%d", iSession), muscles(iMuscle), bestIndex, sprintf("S%d", iSession), muscles(jMuscle), jMuscleR2}];
                end
            end
        end
    end
end


responseName = []; responseName.hot = "CsE_%s_in_uV"; responseName.cold = "SIHIscore_%s";

% Quantify hotspot-coldspot levels
for subject = subjects
    for hemisphere = ["L", "R"]
        iSession = 2;
        exp_id = sessionNames.(hemisphere)(iSession);
        for fromSpot = ["hot", "cold"]
            fn_fromSpot = sprintf('%s/%s/results/exp_%s/r2/mesh_mesh0/roi_midlayer_%s/%s/sigmoid4/r2_roi_data.hdf5', ...
                basepath, subject, exp_id, lower(hemisphere), sprintf(responseName.(fromSpot), "FDI"));
            if ~isfile(fn_fromSpot)
                fprintf(' did not find %s\n', fn_fromSpot)
                continue
            end
            data = reject_R2_outliers_if_needed(h5read(fn_fromSpot, '/data/tris/c_E_mag'));
            data(isnan(data)) = 0;
            [selfR2, bestIndex] = max(data);

            for toSpot = ["hot", "cold"]
                if fromSpot == "hot" && toSpot == "hot"
                    continue
                end
                fn_toSpot = sprintf('%s/%s/results/exp_%s/r2/mesh_mesh0/roi_midlayer_%s/%s/sigmoid4/r2_roi_data.hdf5', ...
                    basepath, subject, exp_id, lower(hemisphere), sprintf(responseName.(toSpot), "FDI"));
                if ~isfile(fn_toSpot)
                    fprintf(' did not find %s\n', fn_toSpot)
                    continue
                end
                data = reject_R2_outliers_if_needed(h5read(fn_toSpot, '/data/tris/c_E_mag'));
                data(isnan(data)) = 0;
                pickedR2 = data(bestIndex);
                fprintf(' %s/%s  :  S%d %s S%d %s\n', subject, hemisphere, iSession, "FDI", iSession, muscles(jMuscle))
                ResultTable = [ResultTable; {subject, hemisphere, fromSpot, toSpot, sprintf("S%d", iSession), "FDI", bestIndex, sprintf("S%d", iSession), "FDI", pickedR2}];
            end
        end
    end
end



%% P L O T
hemisphere = "L";
cfg = [];
cfg.dotpos.S1 = 1.4; cfg.dotpos.S2 = 1.6;
cfg.medpos.S1 = 0.65; cfg.medpos.S2 = 2.35;
cfg.medianColor = 'r';
fig = figure(Position=[50 300 1100 500]);
tiledlayout(3,3)

% Direct results: How good an R² do we actually achieve, and does this
% improve in S2?
% ax = nexttile;
% mask = ResultTable.Hemisphere == hemisphere & ResultTable.from_session == ResultTable.to_session & ResultTable.from_muscle == "ADM" & ResultTable.to_muscle == ResultTable.from_muscle;
% comparisonBoxplot(ResultTable.R2(mask), ResultTable.from_session(mask), ResultTable.Subject(mask), cfg)
% title('ADM MEP', FontSize=12)
% ax.YAxis.Visible = 'on'; ax.YAxis.Color ='none';
% ylabel('$\triangleright_{\bullet}^{*}$', Interpreter='latex', FontSize=30, Rotation=0, Color='k');
%
% nexttile
% mask = ResultTable.Hemisphere == hemisphere & ResultTable.from_session == ResultTable.to_session & ResultTable.from_muscle == "FDI" & ResultTable.to_muscle == ResultTable.from_muscle;
% comparisonBoxplot(ResultTable.R2(mask), ResultTable.from_session(mask), ResultTable.Subject(mask), cfg)
% title('FDI MEP', FontSize=12)
%
% nexttile
% mask = ResultTable.Hemisphere == hemisphere & ResultTable.from_session == ResultTable.to_session & ResultTable.from_muscle == "APB" & ResultTable.to_muscle == ResultTable.from_muscle;
% comparisonBoxplot(ResultTable.R2(mask), ResultTable.from_session(mask), ResultTable.Subject(mask), cfg)
% title('APB MEP', FontSize=12)

for fromMuscle = muscles
    for toMuscle = muscles
        nexttile

        if fromMuscle ~= toMuscle
            cfg.medianColor = 'k';
        else
            cfg.medianColor = muscleColor.(fromMuscle); % = toMuscle
        end

        mask = ResultTable.from_spot == "hot" & ResultTable.to_spot == "hot" & ResultTable.Hemisphere == hemisphere & ResultTable.from_session == ResultTable.to_session & ResultTable.from_muscle == fromMuscle & ResultTable.to_muscle == toMuscle;
        comparisonBoxplot(ResultTable.R2(mask), ResultTable.from_session(mask), ResultTable.Subject(mask), cfg)

        % This shows the following:
        % How well does a triangle selected to predict MEPs in one muscle
        % predict MEPs in other muscles?
        % -> Generally, should be worse
        % -> However, if the MEPs in the original muscle were poor, the other
        % muscle may actually give a better spot
        m = fromMuscle;
        mask = ResultTable.from_spot == "hot" & ResultTable.to_spot == "hot" & ResultTable.Hemisphere == hemisphere & ResultTable.from_session == ResultTable.to_session & ResultTable.from_muscle == m & ResultTable.to_muscle == m;
        med = grpstats(ResultTable.R2(mask), ResultTable.from_session(mask), @median);
        plot(med(1), cfg.medpos.S1, 'v', MarkerFaceColor=muscleColor.(m), MarkerEdgeColor=muscleColor.(m), MarkerSize=3)
        plot(med(2), cfg.medpos.S2, '^', MarkerFaceColor=muscleColor.(m), MarkerEdgeColor=muscleColor.(m), MarkerSize=3)

        % This shows the following:
        m = toMuscle;
        mask = ResultTable.from_spot == "hot" & ResultTable.to_spot == "hot" & ResultTable.Hemisphere == hemisphere & ResultTable.from_session == ResultTable.to_session & ResultTable.from_muscle == m & ResultTable.to_muscle == m;
        med = grpstats(ResultTable.R2(mask), ResultTable.from_session(mask), @median);
        plot(med(1), cfg.medpos.S1, 'v', MarkerFaceColor=muscleColor.(m), MarkerEdgeColor=muscleColor.(m), MarkerSize=3)
        plot(med(2), cfg.medpos.S2, '^', MarkerFaceColor=muscleColor.(m), MarkerEdgeColor=muscleColor.(m), MarkerSize=3)
    end
end

title(nexttile(1), 'ADM MEP', FontSize=12)
title(nexttile(2), 'FDI MEP', FontSize=12)
title(nexttile(3), 'APB MEP', FontSize=12)

ax = nexttile(1);
ax.YAxis.Visible = 'on'; ax.YAxis.Color ='none';
ylabel('$\triangleright_{\stackrel{ADM}{}}^{*}$', Interpreter='latex', FontSize=30, Rotation=0, Color='k');
ax = nexttile(4);
ax.YAxis.Visible = 'on'; ax.YAxis.Color ='none';
ylabel('$\triangleright_{\stackrel{FDI}{}}^{*}$', Interpreter='latex', FontSize=30, Rotation=0, Color='k');
ax = nexttile(7);
ax.YAxis.Visible = 'on'; ax.YAxis.Color ='none';
ylabel('$\triangleright_{\stackrel{APB}{}}^{*}$', Interpreter='latex', FontSize=30, Rotation=0, Color='k');

xlabel(nexttile(7), '$R^2$', Interpreter='latex', FontSize=15)
xlabel(nexttile(8), '$R^2$', Interpreter='latex', FontSize=15)
xlabel(nexttile(9), '$R^2$', Interpreter='latex', FontSize=15)


%% PLOT 2
diagColor = "k"; %"#a0a0a0";
dotSize = 10;
fig = figure(Position=[50 300 860 790]);
hemisphere = "L";
selectedSession = "S2";
layout = tiledlayout(4,4);

for iFromMuscle = 1:3
    for iToMuscle = 1:3
        toMuscle = muscles(iToMuscle);
        fromMuscle = muscles(iFromMuscle);
        layoutIndex = 4*(iFromMuscle-1) + iToMuscle;
        cax = nexttile(layoutIndex);
        mask = ResultTable.Hemisphere == hemisphere &  ResultTable.from_spot == "hot" & ResultTable.to_spot == "hot" & ResultTable.from_session == selectedSession & ResultTable.to_session == selectedSession & ResultTable.from_muscle == fromMuscle & ResultTable.to_muscle == toMuscle;
        data = ResultTable(mask,:);
        if fromMuscle == toMuscle
            % Do a boxplot
            b = boxchart(categorical(data.from_session), double(data.R2), Orientation="horizontal", BoxFaceColor='none', BoxEdgeColor='k', BoxMedianLineColor="k", BoxWidth=0.3);
            hold on
            plot(data.R2, 0.4+ones("like", data.R2), 'k.', MarkerSize=dotSize)
            yticks([]); box off; cax.YAxis.Color = "none"; xlim([0 1])
            axis square; set(cax, "Color", "none")

            hax = axes(layout);
            h = histogram(data.R2, BinEdges=0:0.05:1, Normalization="count", FaceColor="none");
            hax.Layout.Tile = layoutIndex;
            axis off;
            set(hax,"Color", "none"); axis square; xlim([0 1]); ylim([0 4*max(h.Values)]);
        else
            mask_to_to = ResultTable.Hemisphere == hemisphere & ResultTable.from_spot == "hot" & ResultTable.to_spot == "hot" & ResultTable.from_session == selectedSession & ResultTable.to_session == selectedSession & ResultTable.from_muscle == toMuscle & ResultTable.to_muscle == toMuscle;
            % Do a scatter-plot
            box on
            plot([0 1], [0 1], Color=diagColor)
            hold on
            x = ResultTable.R2(mask_to_to);
            y = ResultTable.R2(mask);
            [~, sortingX] = sort(ResultTable.Subject(mask_to_to));
            [~, sortingY] = sort(ResultTable.Subject(mask));
            plot(x(sortingX), y(sortingY), 'k.', MarkerSize=dotSize)
            text(0.76, -0.03, "$R^2$", Interpreter="latex",FontSize=12, VerticalAlignment="top", HorizontalAlignment="center")
            ylim([0 1]);
            xlim([0 1])
            axis square; set(cax, "Color", "none")
        end
        
    end
end

cax = nexttile(16);
mask = ResultTable.Hemisphere == hemisphere & ResultTable.from_spot == "cold" & ResultTable.to_spot == "cold" & ResultTable.from_session == selectedSession & ResultTable.to_session == selectedSession & ResultTable.from_muscle == "FDI" & ResultTable.to_muscle == "FDI";
data = ResultTable(mask,:);
b = boxchart(categorical(data.from_session), double(data.R2), Orientation="horizontal", BoxFaceColor='none', BoxEdgeColor='k', BoxMedianLineColor="k", BoxWidth=0.3);
hold on
plot(data.R2, 0.4+ones("like", data.R2), 'k.', MarkerSize=dotSize)
yticks([]); box off; cax.YAxis.Color = "none"; xlim([0 1]); axis square; set(cax, "Color", "none")
hax = axes(layout);
h = histogram(data.R2, BinEdges=0:0.05:1, Normalization="count", FaceColor="none");
hax.Layout.Tile = 16;
axis off;
set(hax,"Color", "none"); axis square; xlim([0 1]); ylim([0 4*max(h.Values)]);

cax = nexttile(8);
mask2 = ResultTable.Hemisphere == hemisphere & ResultTable.from_spot == "hot" & ResultTable.to_spot == "cold" & ResultTable.from_session == selectedSession & ResultTable.to_session == selectedSession & ResultTable.from_muscle == "FDI" & ResultTable.to_muscle == "FDI";
x = ResultTable.R2(mask);
y = ResultTable.R2(mask2);
[~, sortingX] = sort(ResultTable.Subject(mask));
[~, sortingY] = sort(ResultTable.Subject(mask2));
box on
plot([0 1], [0 1], Color=diagColor)
hold on
plot(x(sortingX), y(sortingY), 'k.', MarkerSize=dotSize)
ylim([0 1]); xlim([0 1]); axis square; set(cax, "Color", "none")

cax = nexttile(14);
mask = ResultTable.Hemisphere == hemisphere & ResultTable.from_spot == "hot" & ResultTable.to_spot == "hot" & ResultTable.from_session == selectedSession & ResultTable.to_session == selectedSession & ResultTable.from_muscle == "FDI" & ResultTable.to_muscle == "FDI";
mask2 = ResultTable.Hemisphere == hemisphere & ResultTable.from_spot == "cold" & ResultTable.to_spot == "hot" & ResultTable.from_session == selectedSession & ResultTable.to_session == selectedSession & ResultTable.from_muscle == "FDI" & ResultTable.to_muscle == "FDI";
x = ResultTable.R2(mask);
y = ResultTable.R2(mask2);
[~, sortingX] = sort(ResultTable.Subject(mask));
[~, sortingY] = sort(ResultTable.Subject(mask2));
box on
plot([0 1], [0 1], Color=diagColor)
hold on
plot(x(sortingX), y(sortingY), 'k.', MarkerSize=dotSize)
ylim([0 1]); xlim([0 1]); axis square; set(cax, "Color", "none")


template = "\\bfMEP %s";
for i=1:3
    title(nexttile(i), sprintf(template, muscles(i)), FontSize=15)
end
cax = nexttile(4);
axis off;
title(cax, "\bf SIHI FDI", FontSize=15)

template = "$\\triangleright_{\\stackrel{%s}{}}^{+}$";
for i=1:3
    ax = nexttile(4*(i-1)+1);
    if i == 1
        y = 1;
    else
        y = 0.5;
    end
    text(ax, -0.5, y, sprintf(template, muscles(i)), FontSize=30, Interpreter='latex', Rotation=0, Color='k', HorizontalAlignment="center", VerticalAlignment="baseline");
end

cax = nexttile(13);
axis off;
text(cax, -0.5, 0.5, "$\triangleright_{\stackrel{FDI}{}}^{-}$", FontSize=30, Interpreter='latex', Rotation=0, Color='k', HorizontalAlignment="center", VerticalAlignment="baseline");


exportgraphics(fig, sprintf('B:/Projects/2023-01 HighDef/Results/Evaluation/group_cross_R2-%s.pdf', hemisphere))
exportgraphics(fig, sprintf('B:/Projects/2023-01 HighDef/Results/Evaluation/group_cross_R2-%s.png', hemisphere))






%% PLOT 3
dotSize = 10;


hemisphere_colors = []; hemisphere_colors.L = "b"; hemisphere_colors.R = "r";
dotpos = []; dotpos.R = 0.9; dotpos.L = 1.1;
boxpos = []; boxpos.R = 0.5; boxpos.L = 1.5;
boxwidth = 0.4;
names = []; 
names.hotcold = ["Distance between Hotspot"; "and Coldspot of FDI [mm]"]; 
names.ADMvsFDI = ["Distance between Hotspots"; "of ADM and FDI [mm]"]; 
names.ADMvsAPB = ["Distance between Hotspots"; "of ADM and APB [mm]"]; 
names.FDIvsAPB = ["Distance between Hotspots"; "of FDI and APB [mm]"];

diagColor = "k"; %"#a0a0a0";
dotSize = 10;
fig = figure(Position=[50 50 950 1150]);


selectedSession = "S2";
layout = tiledlayout(5,4);

for iFromMuscle = 1:3
    for iToMuscle = 1:3
        toMuscle = muscles(iToMuscle);
        fromMuscle = muscles(iFromMuscle);
        layoutIndex = 4*(iFromMuscle) + iToMuscle;
        cax = nexttile(layoutIndex);
        mask = ResultTable.from_spot == "hot" & ResultTable.to_spot == "hot" & ResultTable.from_session == selectedSession & ResultTable.to_session == selectedSession & ResultTable.from_muscle == fromMuscle & ResultTable.to_muscle == toMuscle;
        data = ResultTable(mask,:);
        if fromMuscle == toMuscle
            % Do a boxplot
            %b = boxchart(categorical(data.from_session), double(data.R2), Orientation="horizontal", BoxFaceColor='none', BoxEdgeColor='k', BoxMedianLineColor="k", BoxWidth=0.3);
            hold on
            %plot(data.R2, 0.4+ones("like", data.R2), 'k.', MarkerSize=dotSize)
            bihemisphere_boxplot(data, dotpos, dotSize, hemisphere_colors, boxpos, boxwidth)
            yticks([]); box off; cax.YAxis.Color = "none"; xlim([0 1])
            axis square; set(cax, "Color", "none")

            % hax = axes(layout);
            % h = histogram(data.R2, BinEdges=0:0.05:1, Normalization="probability", FaceColor="none");
            % hax.Layout.Tile = layoutIndex;
            % axis off;
            % set(hax,"Color", "none"); axis square; xlim([0 1]); ylim([0 4*max(h.Values)]);
        else
            mask_to_to = ResultTable.from_spot == "hot" & ResultTable.to_spot == "hot" & ResultTable.from_session == selectedSession & ResultTable.to_session == selectedSession & ResultTable.from_muscle == toMuscle & ResultTable.to_muscle == toMuscle;
            % Do a scatter-plot
            box on
            plot([0 1], [0 1], Color=diagColor)
            hold on
            % x = ResultTable.R2(mask_to_to);
            % y = ResultTable.R2(mask);
            % [~, sortingX] = sort(ResultTable.Subject(mask_to_to));
            % [~, sortingY] = sort(ResultTable.Subject(mask));
            % plot(x(sortingX), y(sortingY), 'k.', MarkerSize=dotSize)
            hold on
            for h = ["R", "L"]
                x = ResultTable.R2(mask_to_to & ResultTable.Hemisphere == h);
                y = ResultTable.R2(mask & ResultTable.Hemisphere == h);
                [~, sortingX] = sort(ResultTable.Subject(mask_to_to & ResultTable.Hemisphere == h));
                [~, sortingY] = sort(ResultTable.Subject(mask & ResultTable.Hemisphere == h));
                plot(x(sortingX), y(sortingY), '.', MarkerSize=dotSize, Color=hemisphere_colors.(h))
            end
            %text(0.76, -0.03, "$R^2$", Interpreter="latex",FontSize=12, VerticalAlignment="top", HorizontalAlignment="center")
            ylim([0 1]);
            xlim([0 1])
            axis square; set(cax, "Color", "none")
            xticks([0 0.5 1]); yticks([0 0.5 1])
            %xlabel("R²");
        end
        
    end
end

cax = nexttile(20);
mask = ResultTable.from_spot == "cold" & ResultTable.to_spot == "cold" & ResultTable.from_session == selectedSession & ResultTable.to_session == selectedSession & ResultTable.from_muscle == "FDI" & ResultTable.to_muscle == "FDI";
data = ResultTable(mask,:);
hold on
bihemisphere_boxplot(data, dotpos, dotSize, hemisphere_colors, boxpos, boxwidth)
yticks([]); box off; cax.YAxis.Color = "none"; xlim([0 1]); axis square; set(cax, "Color", "none")

cax = nexttile(12);
plot([0 1], [0 1], Color=diagColor)
hold on
mask = ResultTable.from_spot == "cold" & ResultTable.to_spot == "cold" & ResultTable.from_session == selectedSession & ResultTable.to_session == selectedSession & ResultTable.from_muscle == "FDI" & ResultTable.to_muscle == "FDI";
mask2 = ResultTable.from_spot == "hot" & ResultTable.to_spot == "cold" & ResultTable.from_session == selectedSession & ResultTable.to_session == selectedSession & ResultTable.from_muscle == "FDI" & ResultTable.to_muscle == "FDI";
for h = ["R", "L"]
    x = ResultTable.R2(mask & ResultTable.Hemisphere == h);
    y = ResultTable.R2(mask2 & ResultTable.Hemisphere == h);
    [~, sortingX] = sort(ResultTable.Subject(mask & ResultTable.Hemisphere == h));
    [~, sortingY] = sort(ResultTable.Subject(mask & ResultTable.Hemisphere == h));
    plot(x(sortingX), y(sortingY), '.', MarkerSize=dotSize, Color=hemisphere_colors.(h))
end
box on
ylim([0 1]); xlim([0 1]); axis square; set(cax, "Color", "none")
xticks([0 0.5 1]); yticks([0 0.5 1])

cax = nexttile(18);
plot([0 1], [0 1], Color=diagColor)
hold on
mask = ResultTable.from_spot == "hot" & ResultTable.to_spot == "hot" & ResultTable.from_session == selectedSession & ResultTable.to_session == selectedSession & ResultTable.from_muscle == "FDI" & ResultTable.to_muscle == "FDI";
mask2 = ResultTable.from_spot == "cold" & ResultTable.to_spot == "hot" & ResultTable.from_session == selectedSession & ResultTable.to_session == selectedSession & ResultTable.from_muscle == "FDI" & ResultTable.to_muscle == "FDI";
for h = ["R", "L"]
    x = ResultTable.R2(mask & ResultTable.Hemisphere == h);
    y = ResultTable.R2(mask2 & ResultTable.Hemisphere == h);
    [~, sortingX] = sort(ResultTable.Subject(mask & ResultTable.Hemisphere == h));
    [~, sortingY] = sort(ResultTable.Subject(mask & ResultTable.Hemisphere == h));
    plot(x(sortingX), y(sortingY), '.', MarkerSize=dotSize, Color=hemisphere_colors.(h))
end
box on
ylim([0 1]); xlim([0 1]); axis square; set(cax, "Color", "none")
xticks([0 0.5 1]); yticks([0 0.5 1])


template = "\\bfMEP %s";
for i=1:3
    title(nexttile(i+4), sprintf(template, muscles(i)), FontSize=15)
end
cax = nexttile(8);
axis off;
title(cax, "\bf SIHI FDI", FontSize=15)

template = "$\\triangleright_{\\stackrel{%s}{}}^{+}$";
for i=2:3
    ax = nexttile(4*i+1);
    text(ax, -0.7, 0.5, sprintf(template, muscles(i)), FontSize=30, Interpreter='latex', Rotation=0, Color='k', HorizontalAlignment="left", VerticalAlignment="baseline");
end

ax = axes(layout);
axis square
xlim([0 1]); ylim([0 1])
axis off;
ax.Layout.Tile = 5;
text(ax, -0.7, 0.5, sprintf(template, muscles(1)), FontSize=30, Interpreter='latex', Rotation=0, Color='k', HorizontalAlignment="left", VerticalAlignment="baseline");



cax = nexttile(17);
axis off;
text(cax, -0.7, 0.5, "$\triangleright_{\stackrel{FDI}{}}^{-}$", FontSize=30, Interpreter='latex', Rotation=0, Color='k', HorizontalAlignment="left", VerticalAlignment="baseline");
axis square;



hemisphere = "L";
ID_example = 'sub-004';

% TODO: add example maps in tiles 1:4
exp_id = sessionNames.(hemisphere)(end);
mesh_id = 'mesh0';
roi = sprintf('midlayer_%s', lower(hemisphere));
zoom_on = [];

for i = [2 1 3 4]
    if i == 4
        response = "SIHIscore_FDI";
    else
        response = sprintf("CsE_%s_in_uV", muscles(i));
    end
    ax = nexttile(i);
    [ts, max_loc, whichTriangle, maxR2, corners] = plot_R2_map(sprintf('%s/%s/results/exp_%s/r2/mesh_%s/roi_%s/%s/sigmoid4', basepath, ID_example, exp_id, mesh_id, roi, response), true, 'mag');
    if i == 2
        zoom_on = max_loc;
    end
    hold on
    if hemisphere == "R"
        az = 90;
        el = 45;
    else
        az = 0;
        el = 90;
    end
    format_axis(zoom_on, 3, az, el)
    trisurf([1 2 3], corners(1,:), corners(2,:), corners(3,:), 0, EdgeColor="k", FaceColor="none", LineWidth=1)
    axis off
end




exportgraphics(fig, sprintf('B:/Projects/2023-01 HighDef/Results/Evaluation/group_cross_R2.pdf'))
exportgraphics(fig, sprintf('B:/Projects/2023-01 HighDef/Results/Evaluation/group_cross_R2.png'))








% OLD: %% P L O T
% hemisphere = "R";
% cfg = [];
% cfg.dotpos.S1 = 1.4; cfg.dotpos.S2 = 1.6;
% fig = figure(Position=[50 300 1100 900]);
% tiledlayout(4,3)
%
% % Direct results: How good an R² do we actually achieve, and does this
% % improve in S2?
% ax = nexttile;
% mask = ResultTable.Hemisphere == hemisphere & ResultTable.from_session == ResultTable.to_session & ResultTable.from_muscle == "ADM" & ResultTable.to_muscle == ResultTable.from_muscle;
% comparisonBoxplot(ResultTable.R2(mask), ResultTable.from_session(mask), ResultTable.Subject(mask), cfg)
% title('ADM MEP', FontSize=12)
% ax.YAxis.Visible = 'on'; ax.YAxis.Color ='none';
% ylabel('$\triangleright_{\bullet}^{*}$', Interpreter='latex', FontSize=30, Rotation=0, Color='k');
%
% nexttile
% mask = ResultTable.Hemisphere == hemisphere & ResultTable.from_session == ResultTable.to_session & ResultTable.from_muscle == "FDI" & ResultTable.to_muscle == ResultTable.from_muscle;
% comparisonBoxplot(ResultTable.R2(mask), ResultTable.from_session(mask), ResultTable.Subject(mask), cfg)
% title('FDI MEP', FontSize=12)
%
% nexttile
% mask = ResultTable.Hemisphere == hemisphere & ResultTable.from_session == ResultTable.to_session & ResultTable.from_muscle == "APB" & ResultTable.to_muscle == ResultTable.from_muscle;
% comparisonBoxplot(ResultTable.R2(mask), ResultTable.from_session(mask), ResultTable.Subject(mask), cfg)
% title('APB MEP', FontSize=12)
%
% for fromMuscle = muscles
%     for toMuscle = muscles
%         nexttile
%         if fromMuscle == toMuscle
%             mask = ResultTable.Hemisphere == hemisphere & ResultTable.from_muscle == fromMuscle & ResultTable.to_muscle == toMuscle & ResultTable.from_session ~= ResultTable.to_session;
%             comparisonBoxplot(ResultTable.R2(mask), ResultTable.from_session(mask), ResultTable.Subject(mask), cfg)
%             text(ones(1,2), categorical(["S1", "S2"]), ["on S2", "on S1"], HorizontalAlignment="left");
%             xlabel('$R^2_\triangleright / \max R^2$', Interpreter='latex')
%         else
%             mask = ResultTable.Hemisphere == hemisphere & ResultTable.from_session == ResultTable.to_session & ResultTable.from_muscle == fromMuscle & ResultTable.to_muscle == toMuscle;
%             comparisonBoxplot(ResultTable.R2(mask), ResultTable.from_session(mask), ResultTable.Subject(mask), cfg)
%         end
%     end
% end
%
% ax = nexttile(4);
% ax.YAxis.Visible = 'on'; ax.YAxis.Color ='none';
% ylabel('$\triangleright_{\stackrel{ADM}{}}^{*}$', Interpreter='latex', FontSize=30, Rotation=0, Color='k');
% ax = nexttile(7);
% ax.YAxis.Visible = 'on'; ax.YAxis.Color ='none';
% ylabel('$\triangleright_{\stackrel{FDI}{}}^{*}$', Interpreter='latex', FontSize=30, Rotation=0, Color='k');
% ax = nexttile(10);
% ax.YAxis.Visible = 'on'; ax.YAxis.Color ='none';
% ylabel('$\triangleright_{\stackrel{APB}{}}^{*}$', Interpreter='latex', FontSize=30, Rotation=0, Color='k');



function bihemisphere_boxplot(data, dotpos, dotSize, hemisphere_colors, boxpos, boxwidth)
for h = ["L", "R"]
    d = data.R2(data.Hemisphere == h);
    plot(d, dotpos.(h), '.', Color=hemisphere_colors.(h), MarkerSize=dotSize)
    plot(median(d) * ones(1,2), boxpos.(h) + boxwidth.*[-0.5 0.5], '-', Color=hemisphere_colors.(h), LineWidth=2)
    q1 = quantile(d, 0.25); q3 = quantile(d, 0.75);
    w1 = min(d(d >= q1-(1.5*(q3-q1))));
    w3 = max(d(d <= q3+(1.5*(q3-q1))));
    plot([q1 q1 q3 q3 q1], boxpos.(h) + boxwidth.*[-0.5 0.5 0.5 -0.5 -0.5], 'k-')
    plot([q1 w1], boxpos.(h) .* ones(1,2), 'k-'); plot([w1 w1], boxpos.(h) + boxwidth.*[-0.25 0.25], 'k-')
    plot([q3 w3], boxpos.(h) .* ones(1,2), 'k-'); plot([w3 w3], boxpos.(h) + boxwidth.*[-0.25 0.25], 'k-')
    text(-0.05, boxpos.(h), h, Color=hemisphere_colors.(h), HorizontalAlignment="center", FontSize=12)
end
end



function comparisonBoxplot(data, grouping, subjects, cfg)
boxchart(categorical(grouping), data, Orientation="horizontal", BoxFaceColor='none', BoxEdgeColor='k', BoxMedianLineColor=cfg.medianColor)
hold on
plot(data, arrayfun(@(s) cfg.dotpos.(s), grouping), 'k.')
for subject = unique(subjects)'
    x = data(subjects == subject);
    y = arrayfun(@(s) cfg.dotpos.(s), grouping(subjects == subject));
    plot(x, y, 'k-')
end
set(gca, 'Color', 'none')
box off
ax = gca;
ax.YAxis.Visible = 'off';
ax.YAxis.Direction = 'reverse';
text(zeros(1,2), categorical(["S1", "S2"]), ["S1", "S2"], HorizontalAlignment="right");
xlim([0 1])
end

function format_axis(center, zoom, az, el)
arguments
    center = [];
    zoom = 1;
    az = 0;
    el = 90;
end
limit_range = max([range(xlim) range(ylim) range(zlim)]);
window = [-limit_range limit_range] ./ 2;
window = window ./ zoom;

if isempty(center)
    xlim(mean(xlim) + window); ylim(mean(ylim) + window); zlim(max(zlim) + window - max(window))
else
    xlim(center(1) + window); ylim(center(2) + window); zlim(center(3) + window)
end
xlabel('\leftarrow Left [mm]')
ylabel('Front [mm] \rightarrow')
%cb = colorbar(Location="north"); 
%ylabel(cb, 'R^2', Rotation=0)
colormap(warm_colormap())
set(gca(), 'Color', 'none')
view(az, el)
axis square
end
