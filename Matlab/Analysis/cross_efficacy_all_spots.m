

data = readtable("B:/Projects/2023-01 HighDef/Results/Evaluation/group_level_spot_cross_relevance_dominant.csv");

left_handed = ["sub-006", "sub-016", "sub-017"];

individual_data = readtable("D:/HighDef-operate/Figures/contrasts.csv", TextType="string");
individual_data.Dominant = repmat("nondominant", size(individual_data, 1), 1);
individual_data.Dominant((individual_data.Hemisphere == "L" & ~ismember(individual_data.Subject, left_handed)) | (individual_data.Hemisphere == "R" & ismember(individual_data.Subject, left_handed))) = "dominant";

%%
contrasts = [                       "hFDI_vs_hADM", "hAPB_vs_hADM", ...
    "cFDI_vs_hFDI", "hADM_vs_hFDI",                 "hAPB_vs_hFDI", ...
    "hADM_vs_hAPB", "hFDI_vs_hAPB"];

significance = 0.05;
% Multiple testing correction: Benjamini Hochberg:

p_values = [];
for hemisphere = ["dominant", "nondominant"]
    p_values.(hemisphere) = arrayfun(@(c) data.p_value(strcmpi(data.Contrast, c) & strcmpi(data.Coil, extractAfter(c, "_vs_")) & strcmpi(data.Hemisphere, hemisphere)), contrasts);
end

all_ps = [p_values.dominant p_values.nondominant];



[p_s, indexing] = sort(all_ps);
% Bring along labeling of these p values
contrasts_resorted = [contrasts contrasts];
contrasts_resorted = contrasts_resorted(indexing);
hemispheres_resorted = [repmat("L", 1, length(contrasts)) repmat("R", 1, length(contrasts))];
hemispheres_resorted = hemispheres_resorted(indexing);


m = length(p_s);
thresholds = significance .* ((1:m) ./ m);
%disp(p_s)
%disp(thresholds)
accepted = cumsum(p_s > thresholds) < 1;

statistics = [];
for hemisphere = ["dominant", "nondominant"]
    for iContrast = 1:length(contrasts)
        contrast = contrasts(iContrast);
        statistics.(hemisphere).(contrast).p = data.p_value(strcmpi(data.Contrast, contrast) & strcmpi(data.Coil, extractAfter(contrast, "_vs_")) & strcmpi(data.Hemisphere, hemisphere));
        statistics.(hemisphere).(contrast).t = data.t_value(strcmpi(data.Contrast, contrast) & strcmpi(data.Coil, extractAfter(contrast, "_vs_")) & strcmpi(data.Hemisphere, hemisphere));
        statistics.(hemisphere).(contrast).accepted = accepted(strcmpi(contrasts_resorted, contrast) & strcmpi(hemispheres_resorted, hemisphere));
    end
end


%%


fig = figure(Position=[50 50 820 900]);
tileIndices = [3,4,5,6,8,10,11];



barwidth = 0.5;
whiskerendwidth = 0.1;
layout = tiledlayout(6,4, TileSpacing="tight");
hemispheres = ["dominant", "nondominant"];
for iHemisphere = 1:2
    tile_offset = (iHemisphere-1)*12;
    hemisphere = hemispheres(iHemisphere);

    for iContrast = 1:length(contrasts)
        contrast = contrasts(iContrast);
        nexttile(layout, tile_offset + tileIndices(iContrast));

        contrastants = split(contrast, "_vs_");
        hold on

        transform = @(x) x.^4;
        basevalue = 0;
        if strcmpi(contrast, "cFDI_vs_hFDI")
            transform = @(x) exp(-x);
            basevalue = 1;
        end

        mask_ref = strcmpi(data.Contrast, contrast) & strcmpi(data.Hemisphere, hemisphere) & strcmpi(data.Coil, contrastants(1));
        mask_other = strcmpi(data.Contrast, contrast) & strcmpi(data.Hemisphere, hemisphere) & strcmpi(data.Coil, contrastants(2));

        %plot(1, transform(data.Estimate(mask_ref)), 'ko', MarkerFaceColor='k');
        barheight_ref = transform(data.Estimate(mask_ref));
        %yline(barheight_ref)
        %plot(2, transform(data.Estimate(mask_ref) + data.Estimate(mask_other)), 'ko', MarkerFaceColor='k');
        barheight = transform(data.Estimate(mask_ref) + data.Estimate(mask_other));
        if startsWith(contrast, "c")
            fill(1+0.5.*barwidth.*[-1 -1 1 1], [basevalue, barheight_ref, barheight_ref, basevalue], "b", EdgeColor="k", FaceAlpha=0.7)
            fill(2+0.5.*barwidth.*[-1 -1 1 1], [basevalue, barheight, barheight, basevalue], "b", EdgeColor="k", FaceAlpha=0.7)
        else
            plot(1+0.5.*barwidth.*[-1 -1 1 1], [basevalue, barheight_ref, barheight_ref, basevalue], "k")
            plot(2+0.5.*barwidth.*[-1 -1 1 1], [basevalue, barheight, barheight, basevalue], "k")
        end

        %plot([1 1], transform([data.Lower(mask_ref) data.Upper(mask_ref)]), 'k')
        %^ this is not very meaningful!
        plot([2 2], transform(data.Estimate(mask_ref) + [data.Lower(mask_other) data.Upper(mask_other)]), 'k')
        plot(2 + 0.5.*whiskerendwidth.*[-1 1], transform(data.Estimate(mask_ref) + data.Lower(mask_other)) * ones(1,2), 'k')
        plot(2 + 0.5.*whiskerendwidth.*[-1 1], transform(data.Estimate(mask_ref) + data.Upper(mask_other)) * ones(1,2), 'k')



        [p_string, p_stars] = encode_p(statistics.(hemisphere).(contrast).p);
    
        xticks([1 2]); xticklabels(contrastants)
        xlim([0, 3])
        if strcmpi(contrast, "cFDI_vs_hFDI")
            text(2, transform(data.Estimate(mask_ref)) - 0.0, p_stars, HorizontalAlignment="center", VerticalAlignment="top", FontSize=15)
            text(2, transform(data.Estimate(mask_ref)) - 0.15, p_string, HorizontalAlignment="center", VerticalAlignment="top", FontSize=9)
            ylim([0 1])
        else
            text(2, 1.07*transform(data.Estimate(mask_ref)), p_stars, HorizontalAlignment="center", FontSize=15)
            text(2, 1.3*transform(data.Estimate(mask_ref)), p_string, HorizontalAlignment="center", FontSize=9)
            ylim([0 1.5*transform(data.Estimate(mask_ref))])
        end

        if startsWith(contrast, "h")
            ylabel(sprintf("MEP in %s [\\muV]", extractAfter(contrastants(1), 1)))
        else
            ylabel("SIHI")
        end

        box on
        axis square
        set(gca, Color="none")
    end
end


exportgraphics(fig, "B:/Projects/2023-01 HighDef/Results/Evaluation/cross_efficacy_population.pdf")
exportgraphics(fig, "B:/Projects/2023-01 HighDef/Results/Evaluation/cross_efficacy_population.png", Resolution=600)






%% Using linear model (from R):

for hemisphere = ["dominant", "nondominant"]
    fig = figure(Position=[50 50 190 200]);
    hold on
    contrast = "cFDI_vs_hFDI";
    contrastants = split(contrast, "_vs_");

    transform = @(x) exp(-x);
    basevalue = 1;

    mask_ref = strcmpi(data.Contrast, contrast) & strcmpi(data.Hemisphere, hemisphere) & strcmpi(data.Coil, contrastants(1));
    mask_other = strcmpi(data.Contrast, contrast) & strcmpi(data.Hemisphere, hemisphere) & strcmpi(data.Coil, contrastants(2));
    barheight_ref = transform(data.Estimate(mask_ref));
    fprintf("SIHI %s on coldspot = %0.3f\n", hemisphere, barheight_ref)
    %yline(barheight_ref)
    barheight = transform(data.Estimate(mask_ref) + data.Estimate(mask_other));
    fprintf("SIHI %s on hotspot  = %0.3f\n", hemisphere, barheight)
    fprintf("        Difference = %0.3f\n", barheight - barheight_ref)

    fill(1+0.5.*barwidth.*[-1 -1 1 1], [basevalue, barheight_ref, barheight_ref, basevalue], "b", EdgeColor="k", FaceAlpha=0.7)
    fill(2+0.5.*barwidth.*[-1 -1 1 1], [basevalue, barheight, barheight, basevalue], "b", EdgeColor="k", FaceAlpha=0.7)

    plot([2 2], transform(data.Estimate(mask_ref) + [data.Lower(mask_other) data.Upper(mask_other)]), 'k')
    plot(2 + 0.5.*whiskerendwidth.*[-1 1], transform(data.Estimate(mask_ref) + data.Lower(mask_other)) * ones(1,2), 'k')
    plot(2 + 0.5.*whiskerendwidth.*[-1 1], transform(data.Estimate(mask_ref) + data.Upper(mask_other)) * ones(1,2), 'k')


    % For review: Add individual data points to the figure
    % Individual data:
    SIHI_rows = individual_data(individual_data.Receiver == "-FDI" & individual_data.Dominant == hemisphere,:);
    SIHI_rows.Response = exp(-SIHI_rows.Response);
    summary_table = groupsummary(SIHI_rows, ["Coil", "Subject"], "mean", "Response");

    %plot(0.5, summary_table.mean_Response(summary_table.Coil == "-FDI"), 'k.')
    %plot(2.5, summary_table.mean_Response(summary_table.Coil == "+FDI"), 'k.')
    % Instead: Show the full individual data in supplementary figure!



    [p_string, p_stars] = encode_p(statistics.(hemisphere).(contrast).p);
    fprintf("%s Significant? %s -- with %s, t=%0.2f\n\n", hemisphere, p_stars, p_string, statistics.(hemisphere).(contrast).t)
    text(2, transform(data.Estimate(mask_ref)) - 0.0, p_stars, HorizontalAlignment="center", VerticalAlignment="top", FontSize=15)
    text(2, transform(data.Estimate(mask_ref)) - 0.15, p_string, HorizontalAlignment="center", VerticalAlignment="top", FontSize=9)

    xticks([1 2]);
    xticklabels(["coldspot", "hotspot"])
    xlim([0, 3])
    ylim([0 1])
    ylabel("SIHI")
    xlabel("Coil on")
    box on
    axis square
    set(gca, Color="none")
    set(gca,xaxisLocation='top')

    exportgraphics(fig, sprintf("B:/Projects/2023-01 HighDef/Results/Evaluation/cross_efficacy_hotncold_%s.pdf", hemisphere))
    exportgraphics(fig, sprintf("B:/Projects/2023-01 HighDef/Results/Evaluation/cross_efficacy_hotncold_%s.png", hemisphere), Resolution=600)
end



%% Direct on data:

for hemisphere = ["dominant", "nondominant"]
    fig = figure(Position=[50 50 190 200]);
    hold on

    transform = @(x) exp(-x);
    basevalue = 1;

    % For review: Add individual data points to the figure
    % Individual data:
    SIHI_rows = individual_data(individual_data.Receiver == "-FDI" & individual_data.Dominant == hemisphere,:);
    SIHI_rows.Response = transform(SIHI_rows.Response);
    summary_table = groupsummary(SIHI_rows, ["Coil", "Subject"], "mean", "Response");
       
    barheight = mean(summary_table.mean_Response(summary_table.Coil == "+FDI"));
    barstd = std(summary_table.mean_Response(summary_table.Coil == "+FDI"));
    barheight_ref = mean(summary_table.mean_Response(summary_table.Coil == "-FDI"));
    barstd_ref = std(summary_table.mean_Response(summary_table.Coil == "+FDI"));


    plot(1.4, summary_table.mean_Response(summary_table.Coil == "-FDI"), 'k.')
    hold on
    plot(1.6, summary_table.mean_Response(summary_table.Coil == "+FDI"), 'k.')

    %plot([1.4 1.6], [summary_table.mean_Response(summary_table.Coil == "-FDI") summary_table.mean_Response(summary_table.Coil == "+FDI")], 'k')

    fill(1+0.5.*barwidth.*[-1 -1 1 1], [basevalue, barheight_ref, barheight_ref, basevalue], "k", EdgeColor="k", FaceAlpha=1)
    fill(2+0.5.*barwidth.*[-1 -1 1 1], [basevalue, barheight, barheight, basevalue], "k", EdgeColor="k", FaceAlpha=1)

    % Whiskers:
    plot([1 1], barheight_ref + [-barstd_ref, barstd_ref], 'k')
    plot(1 + 0.5.*whiskerendwidth.*[-1 1], (barheight_ref + barstd_ref) * ones(1,2), 'k')
    plot(1 + 0.5.*whiskerendwidth.*[-1 1], (barheight_ref - barstd_ref) * ones(1,2), 'k')

    plot([2 2], barheight + [-barstd, barstd], 'k')
    plot(2 + 0.5.*whiskerendwidth.*[-1 1], (barheight + barstd) * ones(1,2), 'k')
    plot(2 + 0.5.*whiskerendwidth.*[-1 1], (barheight - barstd) * ones(1,2), 'k')

    [p, ~, stats] = signrank(summary_table.mean_Response(summary_table.Coil == "+FDI"), summary_table.mean_Response(summary_table.Coil == "-FDI"), Tail="right"); % paired, since its sub-1 vs sub-1

    [p_string, p_stars] = encode_p(p);
    fprintf("%s Significant? %s -- with %s, z=%0.2f\n\n", hemisphere, p_stars, p_string, stats.zval)
    text(2, barheight - barstd - 0.1, p_stars, HorizontalAlignment="center", VerticalAlignment="top", FontSize=15)
    text(2, barheight - barstd - 0.2, p_string, HorizontalAlignment="center", VerticalAlignment="top", FontSize=9)


    xticks([1 2]);
    xticklabels(["coldspot", "hotspot"])
    xlim([0, 3])
    ylim([0 1])
    ylabel("SIHI")
    xlabel("Coil on")
    box on
    axis square
    set(gca, Color="none")
    set(gca,xaxisLocation='top')

    %exportgraphics(fig, sprintf("B:/Projects/2023-01 HighDef/Results/Evaluation/cross_efficacy_hotncold_%s.pdf", hemisphere))
    %exportgraphics(fig, sprintf("B:/Projects/2023-01 HighDef/Results/Evaluation/cross_efficacy_hotncold_%s.png", hemisphere), Resolution=600)
end





function [encoded, stars] = encode_p(p, digits)
arguments
    p (1,1) {mustBeNumeric};
    digits (1,1) {mustBeInteger, mustBePositive} = 3;
end
format = sprintf("%%0.%df", digits);
if p < 10^-digits
    encoded = sprintf(sprintf("p < %s", format), 10^-digits);
else
    encoded = sprintf(sprintf("p = %s", format), p);
end

if p < 0.001
    stars = "***";
elseif p < 0.01
    stars = "**";
elseif p < 0.05
    stars = "*";
elseif p < 0.1
    stars = ".";
else
    stars = "";
end


end