% In this script, check the results of the mapping experiments for all
% subjects, using "conventional/coil space" analysis
addpath(genpath('B:/Projects/2023-01 HighDef/libraries/vetter'))
addpath(genpath('B:/Projects/2023-01 HighDef/libraries/visualization'))
basepath = "//wsl.localhost/Ubuntu-22.04/home/bnplab-admin/TMS_localization/HighDef";
subjects = arrayfun(@(s) string(s.name), dir(fullfile(basepath, 'sub-*')))';
%% 1: Load all data
data.L2R = [];
data.R2L = [];
cog.L2R.hot  = [];
cog.L2R.cold = [];
cog.R2L.hot  = [];
cog.R2L.cold = [];

peak.L2R.hot  = [];
peak.L2R.cold = [];
peak.R2L.hot  = [];
peak.R2L.cold = [];

for subject = subjects
    fig = figure(Position=[100 100 1000 1000]);
    for exp_id = ["L2R", "R2L"]
        filename = sprintf("%s/%s/%s_%s_raw.csv", basepath, subject, subject, exp_id);
        if isfile(filename)
            fprintf("<< %s %s\n", subject, exp_id)
            T = readtable(filename);
            T = T(T.Intensity_percentMSO == max(T.Intensity_percentMSO),:);
            T.Subject = repmat(subject, size(T,1), 1);
            % 1b: Reject trials with overly extreme rotations?
            % Would go here

            % 2: Project trials onto a posterior-anterior x lateral-medial plane (using PCA)
            projected = project_to_standard(T, extractBefore(exp_id, "2"));
            T.Anterior = projected(:,1);
            T.Transversal = projected(:,2);

            data.(exp_id) = [data.(exp_id); T];
            inhibitionRatioWeights = exp(T.SIHIscore_FDI) ./ sum(exp(T.SIHIscore_FDI));
            excitationWeights = T.CsE_FDI_in_uV ./ sum(T.CsE_FDI_in_uV);
            
            
            [mepCOM, mepCov, mepMedian]    = getWeightedMetrics(projected,excitationWeights);
            [sihiCOM, sihiCov, sihiMedian] = getWeightedMetrics(projected,inhibitionRatioWeights);
            cog.(exp_id).hot  = [cog.(exp_id).hot;  array2table(subject, VariableNames="Subject") array2table(mepCOM,  VariableNames=["Anterior", "Transversal"])];
            cog.(exp_id).cold = [cog.(exp_id).cold; array2table(subject, VariableNames="Subject") array2table(sihiCOM, VariableNames=["Anterior", "Transversal"])];

            
            
            gridLength = 150;
            sd = 4;
            maxDistanceOfPeak = 4*sd;
            kernel = @(x,y) (1/(2*pi*sd^2)).*exp(-(x.^2 + y.^2)./(2*sd^2));
            valueThreshold = 10*kernel(1,0);
            for response = ["CsE_FDI_in_uV", "SIHIscore_FDI"]
                ax = nexttile;
                [interpolatedMap, evalX, evalY] = interpolateResponseMap([T.Transversal T.Anterior], T.(response), gridLength, sd, kernel, valueThreshold, false);
                [~, i] = max(interpolatedMap(:));
                x = evalX(:); x = x(i);
                y = evalY(:); y = y(i);

                s=surf(evalX, evalY, interpolatedMap);
                s.EdgeColor="none";
                hold on
                aboveMap = 1.01*max(interpolatedMap(:));
                if startsWith(response, "CsE")
                    plot3(mepCOM(2), mepCOM(1), aboveMap, 'k+', MarkerSize=7, LineWidth=2)
                    plot3(sihiCOM(2), sihiCOM(1), aboveMap, 'k+', MarkerSize=7, LineWidth=0.5)
                    plot3(x, y, aboveMap, "ko")
                    colormap(ax, warm_colormap())

                    peak.(exp_id).hot = [peak.(exp_id).hot; array2table(subject, VariableNames="Subject") array2table([x y], VariableNames=["Transversal", "Anterior"])];
                else
                    plot3(mepCOM(2), mepCOM(1), aboveMap, 'k+', MarkerSize=7, LineWidth=0.5)
                    plot3(sihiCOM(2), sihiCOM(1), aboveMap, 'k+', MarkerSize=7, LineWidth=2)
                    plot3(x, y, aboveMap, "ko")
                    colormap(ax, [flipud(spring_colormap()); cold_colormap()])
                    climit = max(abs(1-interpolatedMap(:)));
                    clim([1-climit 1+climit]);
                    
                    peak.(exp_id).cold = [peak.(exp_id).cold; array2table(subject, VariableNames="Subject") array2table([x y], VariableNames=["Transversal", "Anterior"])];
                end
                colorbar()
                view(2)
                title(sprintf("%s %s %s", subject, exp_id, response), Interpreter="none");
                ax.Color="none";
            end
        end
    end
    drawnow
    
    exportgraphics(fig, sprintf("B:/Projects/2023-01 HighDef/Results/Coil-space/maps-%s.png", subject), ...
            'BackgroundColor', 'none', 'ContentType', 'image', "Resolution", 400)
end



%%
fig = figure(Position=[500 500 500 530]);
tiledlayout(2,2);
for exp_id = ["L2R", "R2L"]
    nexttile;
    plot([cog.(exp_id).hot.Transversal cog.(exp_id).cold.Transversal]', [cog.(exp_id).hot.Anterior cog.(exp_id).cold.Anterior]', 'k-');
    hold on
    plot(cog.(exp_id).hot.Transversal, cog.(exp_id).hot.Anterior, 'r^', MarkerFaceColor="r"); 
    plot(cog.(exp_id).cold.Transversal, cog.(exp_id).cold.Anterior, 'bv', MarkerFaceColor="b");
    if exp_id == "L2R" 
        xlabel("lateral \rightarrow medial [mm]")
    else
        xlabel("medial \leftarrow lateral [mm]")
    end
    ylabel("posterior \rightarrow anterior [mm]")
    
    set(gca, Color="none")
    title(exp_id)
    axis square

    nexttile;
    plot(([cog.(exp_id).hot.Transversal cog.(exp_id).cold.Transversal] - cog.(exp_id).hot.Transversal)', ([cog.(exp_id).hot.Anterior cog.(exp_id).cold.Anterior] - cog.(exp_id).hot.Anterior)', 'k-');
    hold on
    plot(cog.(exp_id).cold.Transversal - cog.(exp_id).hot.Transversal, cog.(exp_id).cold.Anterior - cog.(exp_id).hot.Anterior, 'ko', MarkerFaceColor="k");
    plot(mean([cog.(exp_id).hot.Transversal cog.(exp_id).cold.Transversal] - cog.(exp_id).hot.Transversal), mean([cog.(exp_id).hot.Anterior cog.(exp_id).cold.Anterior] - cog.(exp_id).hot.Anterior), 'r-', LineWidth=2);
    set(gca, Color="none")

    meanOffset = sqrt(mean(cog.(exp_id).cold.Transversal - cog.(exp_id).hot.Transversal)^2 + mean(cog.(exp_id).cold.Anterior - cog.(exp_id).hot.Anterior)^2);
    title(sprintf("%s: mean offset = %3.3f mm", exp_id, meanOffset))
    if exp_id == "L2R" 
        xlabel("lateral \rightarrow medial [mm]")
    else
        xlabel("medial \leftarrow lateral [mm]")
    end
    ylabel("posterior \rightarrow anterior [mm]")
    axis square
end



exportgraphics(fig, "B:/Projects/2023-01 HighDef/Results/Coil-space/spot-comparisons.pdf", ...
        'BackgroundColor', 'none', 'ContentType', 'vector')



%% Mean offsets for dominant/nondominant:
flip_subject = "sub-006";
valid_subjects = subjects(~ismember(subjects, ["sub-005", "sub-010"]));
mask_right_handed = valid_subjects ~= flip_subject;
mask_left_handed = valid_subjects == flip_subject;

distances.cog.dominant     = nan(length(valid_subjects), 1);
distances.cog.nondominant  = nan(length(valid_subjects), 1);
distances.peak.dominant    = nan(length(valid_subjects), 1);
distances.peak.nondominant = nan(length(valid_subjects), 1);



for dominance = ["dominant", "nondominant"]
    offsets = nan(length(valid_subjects), 2);
    peakOffsets = nan(length(valid_subjects), 2);
    directions = ["Anterior", "Transversal"];
    if dominance == "dominant"
        e1 = "L2R";
        e2 = "R2L";
    else
        e1 = "R2L";
        e2 = "L2R";
    end

    for iDirection = 1:2
        direction = directions(iDirection);
        if direction == "Transversal"
            offsets(:,iDirection) = [cog.(e1).cold.(direction)(mask_right_handed) - cog.(e1).hot.(direction)(mask_right_handed); -(cog.(e2).cold.(direction)(mask_left_handed) - cog.(e2).hot.(direction)(mask_left_handed))];
            peakOffsets(:,iDirection) = [peak.(e1).cold.(direction)(mask_right_handed) - peak.(e1).hot.(direction)(mask_right_handed); -(peak.(e2).cold.(direction)(mask_left_handed) - peak.(e2).hot.(direction)(mask_left_handed))];
        else
            offsets(:,iDirection) = [cog.(e1).cold.(direction)(mask_right_handed) - cog.(e1).hot.(direction)(mask_right_handed); cog.(e2).cold.(direction)(mask_left_handed) - cog.(e2).hot.(direction)(mask_left_handed)];
            peakOffsets(:,iDirection) = [peak.(e1).cold.(direction)(mask_right_handed) - peak.(e1).hot.(direction)(mask_right_handed); peak.(e2).cold.(direction)(mask_left_handed) - peak.(e2).hot.(direction)(mask_left_handed)];
        end
    end

    meanLateralOffset = mean(offsets(:,2));
    meanAnteriorOffset = mean(offsets(:,1));
    meanOffset = mean(vecnorm(offsets'));
    distances.cog.(dominance)  = vecnorm(offsets');
    distances.peak.(dominance) = vecnorm(peakOffsets');
    
    fprintf("Mean lateral offset on %s hemisphere: %3.3f mm\n", dominance, meanLateralOffset)
    fprintf("Mean frontal offset on %s hemisphere: %3.3f mm\n", dominance, meanAnteriorOffset)
    fprintf("Mean offset on %s hemisphere: %3.3f mm\n\n", dominance, meanOffset)
end









