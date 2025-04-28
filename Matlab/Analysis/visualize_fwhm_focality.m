
basepath = 'D:/HighDef-operate/HighDef';

hemisphere = "nondominant";

names = [];
names.righthanded.dominant.exp = "map-L2R";
names.righthanded.dominant.roi = "midlayer_l";
names.righthanded.nondominant.exp = "map-R2L";
names.righthanded.nondominant.roi = "midlayer_r";
names.lefthanded.dominant.exp = "map-R2L";
names.lefthanded.dominant.roi = "midlayer_r";
names.lefthanded.nondominant.exp = "map-L2R";
names.lefthanded.nondominant.roi = "midlayer_l";


%exp_id = "map-R2L";
%roi_id = sprintf("midlayer_%s", lower(extractBefore(extractAfter(exp_id, "-"), "2")));

collected_sim_fwhm_sihi = [];
collected_sim_fwhm_cse = [];
collected_fwhm_cse = [];
collected_fwhm_sihi = [];


%%
% Missing currently: sub-003 SIHI, sub-007 SIHI
% Because even fitting the H0-model failed!

IDs        = ["sub-001", "sub-002", "sub-003", "sub-004", "sub-006", "sub-007", "sub-008",  "sub-009", "sub-011", "sub-012", "sub-013", "sub-014", "sub-015", "sub-016", "sub-017", "sub-019", "sub-022", "sub-023"]; %
lefthanded = ["sub-006", "sub-016", "sub-017"];

all_hotspots  = table();
all_coldspots = table();

for ID = IDs
    if ismember(ID, lefthanded)
        exp_id = names.lefthanded.(hemisphere).exp;
        roi_id = names.lefthanded.(hemisphere).roi;
    else
        exp_id = names.righthanded.(hemisphere).exp;
        roi_id = names.righthanded.(hemisphere).roi;
    end

    file_cse  = sprintf("%s/%s/%s_%s_CsE_FDI_in_uV_simulated_FWHMs_heteroskedastic.mat", basepath, ID, ID, exp_id);
    file_sihi = sprintf("%s/%s/%s_%s_SIHIscore_FDI_simulated_FWHMs_heteroskedastic.mat", basepath, ID, ID, exp_id);

    if isfile(file_cse)
        load(file_cse)
        FWHM_cse = FWHM_raw;
        sim_fwhm_cse = simulated_fwhms;
    else
        FWHM_cse = nan;
        sim_fwhm_cse = nan(1,100);
    end
    if isfile(file_sihi)
        load(file_sihi)
        FWHM_sihi = FWHM_raw;
        sim_fwhm_sihi = simulated_fwhms;
    else
        FWHM_sihi = nan;
        sim_fwhm_sihi = nan(1,100);
    end

    [newrow, ~] = collectSpot(basepath, ID, exp_id, roi_id, "CsE_FDI_in_uV");
    all_hotspots = [all_hotspots; newrow];
    [newrow, ~] = collectSpot(basepath, ID, exp_id, roi_id, "SIHIscore_FDI");
    all_coldspots = [all_coldspots; newrow];

    collected_fwhm_cse  = [collected_fwhm_cse; FWHM_cse];
    collected_fwhm_sihi = [collected_fwhm_sihi; FWHM_sihi];
    collected_sim_fwhm_cse  = [collected_sim_fwhm_cse; sim_fwhm_cse];
    collected_sim_fwhm_sihi = [collected_sim_fwhm_sihi; sim_fwhm_sihi];
    

    fig = figure(Position=[50 50 430 150]);
    bins = 0.0:0.005:0.22;
    histogram(sim_fwhm_sihi, BinEdges=bins, EdgeColor='k', FaceColor='b', FaceAlpha=0.5, Normalization='probability')
    hold on
    histogram(sim_fwhm_cse, BinEdges=bins, EdgeColor='k', FaceColor='r', FaceAlpha=0.5, Normalization='probability')
    xline(FWHM_sihi, Color='b', LineWidth=2)
    xline(FWHM_cse, Color='r', LineWidth=2)
    xlabel("\leftarrow Focality")
    ylabel('probability');
    title(ID);
    drawnow;
    if ~isnan(FWHM_cse)
        p_cse  = mean(sim_fwhm_cse < FWHM_cse);
        m_cse  = encode(p_cse);
    else
        p_cse = nan;
        m_sce = "(H0 model fit failed)";
    end
    if ~isnan(FWHM_sihi)
        p_sihi = mean(sim_fwhm_sihi < FWHM_sihi);
        m_sihi = encode(p_sihi);
    else
        p_sihi = nan;
        m_sce = "(H0 model fit failed)";
    end

    set(gca, Color="none")

    fprintf("%s   Hotspot:   ratio of sim. FWHM < real FWHM: p = %3.f  %s\n", ID, 100*p_cse, m_cse)
    fprintf("          Coldspot:  ratio of sim. FWHM < real FWHM: p = %3.f  %s\n", 100*p_sihi, m_sihi)

    exportgraphics(fig, sprintf('B:/Projects/2023-01 HighDef/Results/Evaluation/%s_focality.pdf', ID))
    exportgraphics(fig, sprintf('B:/Projects/2023-01 HighDef/Results/Evaluation/%s_focality.png', ID))
end


fprintf(sprintf("%s\n", repmat('_', 100,1)))
close all

%%
fig = figure(Position=[50 50 700 700]);
tiledlayout(6,3)
nIDs = length(IDs);
verdicts = []; 
verdicts.hotspot  = table(repmat("", nIDs, 1), zeros(nIDs, 1), zeros(nIDs, 1), repmat("", nIDs, 1), VariableNames=["Subject", "p_hyperfocal", "p_nonfocal", "Verdict"]); 
verdicts.coldspot = table(repmat("", nIDs, 1), zeros(nIDs, 1), zeros(nIDs, 1), repmat("", nIDs, 1), VariableNames=["Subject", "p_hyperfocal", "p_nonfocal", "Verdict"]);

for iID = 1:nIDs %
    ID = IDs(iID);

    if ismember(ID, lefthanded)
        exp_id = names.lefthanded.(hemisphere).exp;
        roi_id = names.lefthanded.(hemisphere).roi;
    else
        exp_id = names.righthanded.(hemisphere).exp;
        roi_id = names.righthanded.(hemisphere).roi;
    end
    file_cse  = sprintf("%s/%s/%s_%s_CsE_FDI_in_uV_simulated_FWHMs_heteroskedastic.mat", basepath, ID, ID, exp_id);
    file_sihi = sprintf("%s/%s/%s_%s_SIHIscore_FDI_simulated_FWHMs_heteroskedastic.mat", basepath, ID, ID, exp_id);

    if isfile(file_cse)
        load(file_cse)
        FWHM_cse = FWHM_raw;
        sim_fwhm_cse = simulated_fwhms;
    else
        FWHM_cse = nan;
        sim_fwhm_cse = nan(1,100);
    end
    if isfile(file_sihi)
        load(file_sihi)
        FWHM_sihi = FWHM_raw;
        sim_fwhm_sihi = simulated_fwhms;
    else
        FWHM_sihi = nan;
        sim_fwhm_sihi = nan(1,100);
    end

    collected_fwhm_cse  = [collected_fwhm_cse; FWHM_cse];
    collected_fwhm_sihi = [collected_fwhm_sihi; FWHM_sihi];
    collected_sim_fwhm_cse  = [collected_sim_fwhm_cse; sim_fwhm_cse];
    collected_sim_fwhm_sihi = [collected_sim_fwhm_sihi; sim_fwhm_sihi];
    
    alpha = 0.5;
    linealpha=1;
    if all_coldspots.R2(iID) < 0.1
        alpha = 0.25;
        linealpha = 0.25;
    end


    ax = nexttile();
    if strcmpi(exp_id, "map-L2R")
        bins = 0.0:0.005:0.28;
    else
        bins = 0.0:0.008:0.40;
    end
    hcold = histogram(sim_fwhm_sihi, BinEdges=bins, EdgeColor='k', FaceColor='b', FaceAlpha=alpha, Normalization='probability');
    hold on
    hhot = histogram(sim_fwhm_cse, BinEdges=bins, EdgeColor='k', FaceColor='r', FaceAlpha=alpha, Normalization='probability');
    ymax = max([hhot.Values hcold.Values]);

    pc = plot(FWHM_sihi.*[1 1], [0 1.05.*ymax], Color='b', LineWidth=2);
    pc.Color = [pc.Color linealpha];
    ph = plot(FWHM_cse.*[1 1],  [0 1.25.*ymax], Color='r', LineWidth=2);
    ph.Color = [ph.Color linealpha];
    xlabel("\leftarrow Focality")
    ylabel('probability');
    title(ID);
    drawnow;
    if ~isnan(FWHM_cse)
        p_cse  = mean(sim_fwhm_cse < FWHM_cse);
        m_cse  = encode(p_cse);
    else
        p_cse = nan;
        m_sce = "(H0 model fit failed)";
    end
    if ~isnan(FWHM_sihi)
        p_sihi = mean(sim_fwhm_sihi < FWHM_sihi);
        m_sihi = encode(p_sihi);
    else
        p_sihi = nan;
        m_sce = "(H0 model fit failed)";
    end
    
    ylim([0 1.5*ymax])

    set(gca, Color="none")
    text(FWHM_cse, 1.25*ymax, sprintf('p = %0.2f', 1-p_cse), Color="r", VerticalAlignment="bottom", HorizontalAlignment="center", FontSize=8)
    if all_coldspots.R2(iID) < 0.1
        greyout = "#808080";

        ax.YAxis.Color=greyout;
        ax.XAxis.Color=greyout;
        ax.Title.Color=greyout;
        hhot.EdgeColor=greyout;
        hcold.EdgeColor=greyout;
    else
        text(FWHM_sihi, 1.05*ymax, sprintf('p = %0.2f', 1-p_sihi), Color="b", VerticalAlignment="bottom", HorizontalAlignment="center", FontSize=8)
    end

    verdicts.hotspot(iID, :)  = {ID, 100*p_cse, 100 - 100*p_cse, m_cse};
    verdicts.coldspot(iID, :) = {ID, 100*p_sihi, 100 - 100*p_sihi, m_sihi};
end
exportgraphics(fig, sprintf('B:/Projects/2023-01 HighDef/Results/Evaluation/all_focality_%s.pdf', hemisphere))
exportgraphics(fig, sprintf('B:/Projects/2023-01 HighDef/Results/Evaluation/all_focality_%s.png', hemisphere), Resolution=400)

verdicts.coldspot(all_coldspots.R2 > 0.1,:)
verdicts.hotspot(all_hotspots.R2 > 0.1,:)
fprintf(sprintf("%s\n", repmat('_', 100,1)))
%%
placeholder = repmat("", nIDs, 1);
final_table = table(placeholder, placeholder, placeholder, placeholder, VariableNames=["Subject", "hot_p", "cold_R2", "cold_p"]);
for iID = 1:nIDs
    final_table.Subject(iID) = IDs(iID);
    hotspot_p   = verdicts.hotspot.p_nonfocal(iID) / 100;
    coldspot_R2 = all_coldspots.R2(iID);
    coldspot_p  = verdicts.coldspot.p_nonfocal(iID) / 100;
    final_table.hot_p(iID) = sprintf("%0.2f %s", hotspot_p, enstar(hotspot_p));
    final_table.cold_R2(iID) = sprintf("%0.3f", coldspot_R2);
    star_for_cold = enstar(coldspot_p);
    if coldspot_R2 < 0.1 && strlength(star_for_cold) > 0
        star_for_cold = "(" + star_for_cold + ")";
    end
    final_table.cold_p(iID) = sprintf("%0.2f %s", coldspot_p, star_for_cold);
end

fprintf(TeXify(final_table))
fprintf("\n")



function s = enstar(p)
    s = "";
    % if p < 0.001
    %     s = "***";
    % elseif p < 0.01
    %     s = "**";
    % elseif p < 0.05
    %     s = "*";
    % end
    if p < 0.05
        s = "*";
    end
end

function s = encode(p)
if p < 0.05
    s = "more focal than expected";
elseif p > 0.95
    s = "certainly nonfocal";
else
    s = "";
end
end

function s = TeXify(t)
s = join(t.Properties.VariableNames, " & ");
s = string(s{:});
nRows = size(t,1);
nColumns = size(t,2);
for iRow = 1:nRows
    rowString = savestring(table2array(t(iRow, 1)));
    for iColumn = 2:nColumns
        rowString = rowString + " & " + savestring(table2array(t(iRow, iColumn)));
    end
    s = s + "\\\\\n" + rowString;
end
end

function s = savestring(v)
    s = string(v);
    if isnumeric(v) && isnan(v)
        s = "nan";
    end
end