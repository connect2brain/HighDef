experiment_names = [];
experiment_names.L = ["map-L", "map-L2R"];
experiment_names.R = ["map-R", "map-R2L"];

basepath = '//wsl.localhost/Ubuntu-22.04/home/bnplab-admin/TMS_localization/HighDef';
subjects = arrayfun(@(s) string(s.name), dir(fullfile(basepath, 'sub-*')))';

allRows = [];
for subject = subjects
    for hemisphere = ["L", "R"]
        roi_id = sprintf("midlayer_%s", lower(hemisphere));
        for muscle = ["ADM", "FDI", "APB"]
            response_id = sprintf("CsE_%s_in_uV", muscle);
            % get S1 map -> get best_triangle_S1 and best_R2_S1
            [S1, S1_map] = collectSpot(basepath, subject, experiment_names.(hemisphere)(1), roi_id, response_id);
            % get S2 map -> get best_triangle_S2 and best_R2_S2
            [S2, S2_map] = collectSpot(basepath, subject, experiment_names.(hemisphere)(2), roi_id, response_id);
            
            if ~isempty(S1) && ~isempty(S2)
                [~, S1_bestTriangle] = max(S1_map);
                [~, S2_bestTriangle] = max(S2_map);

                S1.Session = "S1"; S1.Muscle = muscle; S1.Hemisphere = hemisphere; S1.R2_on_other = S2_map(S1_bestTriangle);
                S2.Session = "S2"; S2.Muscle = muscle; S2.Hemisphere = hemisphere; S2.R2_on_other = S1_map(S2_bestTriangle);
                allRows = [allRows; S1; S2];
            else
                fprintf('%s %s %s not complete\n', subject, hemisphere, muscle)
            end
        end
    end
end


%%
binEdges = linspace(0, 1, 20);
fig = figure(Position=[50 500 1200 600]);
%t = tiledlayout(4,8);

% (A) Distances
indices_S1 = find(allRows.Session == "S1");
indices_S2 = arrayfun(@(i) find(allRows.Session ~= allRows.Session(i) & allRows.Muscle == allRows.Muscle(i) & allRows.Hemisphere == allRows.Hemisphere(i) & allRows.Subject == allRows.Subject(i)), indices_S1);

d_coil = sqrt((allRows.X_coil(indices_S1) - allRows.X_coil(indices_S2)).^2 + (allRows.Y_coil(indices_S1) - allRows.Y_coil(indices_S2)).^2 + (allRows.Z_coil(indices_S1) - allRows.Z_coil(indices_S2)).^2);
d_source = sqrt((allRows.X_source(indices_S1) - allRows.X_source(indices_S2)).^2 + (allRows.Y_source(indices_S1) - allRows.Y_source(indices_S2)).^2 + (allRows.Z_source(indices_S1) - allRows.Z_source(indices_S2)).^2);

[R, P] = corrcoef(d_coil, d_source);
fprintf("Distance in source space varies between %.0f mm and %.0f mm\n", min(d_source), ceil(max(d_source)))
fprintf("Distance in   coil space varies between %.0f mm and %.0f mm\n", min(d_coil), ceil(max(d_coil)))
fprintf("Correlation of source and coil space distances: R² = %0.2f, p = %0.4f\n\n", R(1,2), P(1,2))


ax_scatter = subplot(Position=[0.1 0.1 0.35 0.7]);
scatter(d_source, d_coil, 70, 'kx')
box on; set(ax_scatter, Color="none")
xl_scatter = xlim;
yl_scatter = ylim;
xlabel("Distance on the brain [mm]", FontSize=14)
ylabel("Distance of coil positions [mm]", FontSize=14)

ax_marg_x = subplot(Position=[0.1 0.8 0.35 0.2]);
histogram(d_source, BinEdges=linspace(xl_scatter(1), xl_scatter(2), 20), FaceColor="none")
xlim(xl_scatter)
set(ax_marg_x, Color="none")
ax_marg_x.YAxis.Visible = "off";
ax_marg_x.XAxis.Visible = "off";
box off
linkaxes([ax_scatter, ax_marg_x], "x")

ax_marg_y = subplot(Position=[0.45 0.1 0.1 0.7]);
histogram(d_coil, BinEdges=linspace(yl_scatter(1), yl_scatter(2), 20), Orientation="horizontal", FaceColor="none")
ylim(yl_scatter)
set(ax_marg_y, Color="none")
ax_marg_y.YAxis.Visible = "off";
ax_marg_y.XAxis.Visible = "off";
box off
linkaxes([ax_scatter, ax_marg_y], "y")



% (B) cross R²
ax1 = subplot(Position=[0.65 0.5 0.3 0.3]);
% S1 self R²
h1self = max(histogram(allRows.R2(allRows.Session == "S1"), BinEdges=binEdges, FaceColor="b", Normalization="count", FaceAlpha=0.5).Values);
hold on
% S1 on other
h1on2 = max(histogram(allRows.R2_on_other(allRows.Session == "S1"), BinEdges=binEdges, FaceColor="r", Normalization="count", FaceAlpha=0.5).Values);
box off; set(ax1, Color="none"); ax1.YAxis.Visible = "off"; xlabel("R^2", FontSize=12)
legend(["self", "other"], Location="northwest", FontSize=10);
fprintf("Mann-Whitney-U-test for S1 (H1: median of self > median of other): p = %0.4f\n", ranksum(allRows.R2(allRows.Session == "S1"), allRows.R2_on_other(allRows.Session == "S1"), tail="right"))

ax2 = subplot(Position=[0.65 0.1 0.3 0.3]);
% S2 self R²
h2self = max(histogram(allRows.R2(allRows.Session == "S2"), BinEdges=binEdges, FaceColor="b", Normalization="count", FaceAlpha=0.5).Values);
hold on
% S2 on other
h2on1 = max(histogram(allRows.R2_on_other(allRows.Session == "S2"), BinEdges=binEdges, FaceColor="r", Normalization="count", FaceAlpha=0.5).Values);
box off; set(ax2, Color="none"); ax2.YAxis.Visible = "off"; xlabel("R^2", FontSize=12)
legend(["self", "other"], Location="northwest", FontSize=10);
fprintf("Mann-Whitney-U-test for S2 (H1: median of self > median of other): p = %0.4f\n", ranksum(allRows.R2(allRows.Session == "S2"), allRows.R2_on_other(allRows.Session == "S2"), tail="right"))


ym = max([h1self h2self h1on2 h2on1]);
ylim(ax1, [0 ym]); ylim(ax2, [0 ym])

text(ax1, -0.05, ym/2, "\bf S1", HorizontalAlignment="right", FontSize=14)
text(ax2, -0.05, ym/2, "\bf S2", HorizontalAlignment="right", FontSize=14)

text(ax1, -1.90, ym+0.02, "\bf A", HorizontalAlignment="right", VerticalAlignment="bottom", FontSize=20)
text(ax1, -0.15, ym+0.02, "\bf B", HorizontalAlignment="right", VerticalAlignment="bottom", FontSize=20)

xlim(ax1, [0 1])
xlim(ax2, [0 1])

exportgraphics(fig, sprintf('B:/Projects/2023-01 HighDef/Results/Evaluation/group_hotspot_stability.pdf'))
exportgraphics(fig, sprintf('B:/Projects/2023-01 HighDef/Results/Evaluation/group_hotspot_stability.png'))


