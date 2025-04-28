addpath(genpath('B:/Projects/2023-01 HighDef/libraries/visualization'))


%%
ID = 'sub-014';
direction = 'R2L';
expID = sprintf('map-%s', direction);
roiID = 'midlayer_r';
meshID = 'mesh0';
response = 'SIHIscore_FDI';
n_clusters = 4;
OUT = sprintf('B:/Projects/2023-01 HighDef/Results/R2-Figures/%s_%s_%s_%dclust_', ID, expID, response, n_clusters);

basepath = '//wsl.localhost/Ubuntu-22.04/home/bnplab-admin/TMS_localization/HighDef';

geo_file = sprintf('%s/%s/results/exp_%s/r2/mesh_%s/roi_%s/%s/sigmoid4/r2_roi_geo.hdf5', basepath, ID, expID, meshID, roiID, response);
data_file = sprintf('%s/%s/results/exp_%s/r2/mesh_%s/roi_%s/%s/sigmoid4/r2_roi_data.hdf5', basepath, ID, expID, meshID, roiID, response);
simulation_result = sprintf('%s/%s/results/exp_%s/electric_field/mesh_%s/roi_%s/e.hdf5', basepath, ID, expID, meshID, roiID);
responses = readtable(sprintf('B:/Projects/2023-01 HighDef/Results/%s_%s_raw.csv', ID, direction));

%h5disp(geo_file)
%h5disp(simulation_result)
coordinates = h5read(geo_file, '/mesh/nodes/node_coord');
triangles = h5read(geo_file, '/mesh/elm/triangle_number_list')' + 1;

%%
E_mag = h5read(simulation_result, '/E_mag');

data = h5read(data_file, '/data/tris/c_E_mag');
%data(isnan(data)) = 0;
[bestR2, bestR2Triangle] = max(data);
max_loc_R2 = mean(coordinates(:,triangles(bestR2Triangle,:)), 2);


%%
highSItrials = responses.Intensity_percentMSO == max(responses.Intensity_percentMSO);
R = responses.(response)(highSItrials);
cutoff = prctile(R, 75);
whichTrials = R > cutoff;

directMask = false(size(responses, 1), 1);
directBelowThresholdMask = false(size(responses, 1), 1);
directMask(highSItrials) = whichTrials;
directBelowThresholdMask(highSItrials) = ~whichTrials;

E_mag_thresholded_trials = E_mag(:, directMask);
%E_mag_thresholded_trials = E_mag_thresholded_trials(:, whichTrials);
R_thresholded_trials = R(whichTrials);

%%
distances_between_superthreshold_trials = pdist(E_mag_thresholded_trials', 'correlation');
E_correlations_between_trials = squareform(distances_between_superthreshold_trials);
Z = linkage(distances_between_superthreshold_trials, "average");
dendrogram(Z)
clustering_labels = cluster(Z, maxclust=n_clusters);
tabulate(clustering_labels)
tabulatedClustering = tabulate(clustering_labels);

clusters = struct();
clusters.IDs = unique(clustering_labels)';
clusters.Names = arrayfun(@(i) {sprintf('c%d', i)}, clusters.IDs);
clusters.labels = clustering_labels;
for i = clusters.IDs
    clusters.(clusters.Names{i}) = struct();
end



%%
baseColors = [
    0 0 1; 
    0 1 0;
    1 0 0; 
    1 0 1;
    0 1 1;
    1 1 0;
    ];
colors = baseColors(clusters.labels,:);



v = struct();
v.origin = zeros(1, sum(directMask));
v.p1 = responses.p1(directMask)';
v.p2 = responses.p2(directMask)';
v.p3 = responses.p3(directMask)';
v.x1 = responses.x1(directMask)';
v.x2 = responses.x2(directMask)';
v.x3 = responses.x3(directMask)';
v.y1 = responses.y1(directMask)';
v.y2 = responses.y2(directMask)';
v.y3 = responses.y3(directMask)';
v.z1 = responses.z1(directMask)';
v.z2 = responses.z2(directMask)';
v.z3 = responses.z3(directMask)';


[~, ordering] = sort(clusters.labels);
clusterEdges = cumsum(tabulatedClustering(:,2))+0.5;


fig = figure(Position=[50 50 1300 500]);
subplot(1,2,1)
imagesc(E_correlations_between_trials(ordering, ordering))
hold on
xline(clusterEdges, 'k', LineWidth=1.4)
yline(clusterEdges, 'k', LineWidth=1.4)
for i = clusters.IDs
    p = clusterEdges(i) - 0.5*tabulatedClustering(i,2);
    c = baseColors(i,:);
    plot(p, clusterEdges(1), MarkerEdgeColor=c, MarkerFaceColor=c, Marker='o')
    text(p, clusterEdges(1)-3, sprintf('%d', i), Color=c, VerticalAlignment="bottom", HorizontalAlignment="center")
end

colormap("turbo")
clim([0 2])
colorbar()
axis square

subplot(1,2,2)
ts_roi = trisurf(triangles, coordinates(1,:), coordinates(2,:), coordinates(3,:));
ts_roi.FaceColor = [0.87 0.87 0.87];
ts_roi.FaceAlpha = 1;
ts_roi.EdgeColor = 'none';
%light(Style="infinite",Position=[0 -10 1],Color=[0.1 0.3 0.7]);%[0.1 0.4 0.8]);
light(Style="infinite",Position=[0 10 10],Color=[1 1 1]);
light(Style="infinite",Position=[10 -10 10],Color=[1 1 1]);
lighting gouraud
material dull
hold on
scatter3(responses.p1(directBelowThresholdMask), responses.p2(directBelowThresholdMask), responses.p3(directBelowThresholdMask), 4, 'k', "filled")
scatter3(v.p1, v.p2, v.p3, 10, colors, "filled")
for i = clusters.IDs
    fprintf('Adding cluster %d\n', i)
    for dimension = {'y'}
        x = v.p1 + [v.origin; 3.*v.(sprintf('%s1', dimension{:}))];
        y = v.p2 + [v.origin; 3.*v.(sprintf('%s2', dimension{:}))];
        z = v.p3 + [v.origin; 3.*v.(sprintf('%s3', dimension{:}))];
        plot3(x(:,clusters.labels == i), y(:,clusters.labels == i), z(:,clusters.labels == i), Color=baseColors(i,:), LineWidth=1.4)
    end
end

grid off; box on
set(gca, 'Color', 'none')
xlabel('x'); ylabel('y'); zlabel('z')


r = max([range(xlim) range(ylim) range(zlim)]);
xlim(mean(xlim) + 0.5.*[-r r]); ylim(mean(ylim) + 0.5.*[-r r]); zlim(mean(zlim) + 0.5.*[-r r])
axis square
view(2)

exportgraphics(fig, [OUT 'correlations.png'])



%% S e a r c h  for source of these clusters
% 1: Plot the minimum (or median?) E-field magnitude of the trials in the
% cluster for each Triangle



fig = figure(Position=[50 550 1300 600]);
tiledlayout(2,3, 'TileSpacing', 'tight', Padding='tight'); 
for i = clusters.IDs
    nexttile;
    clusterMask = clusters.labels == i;
    %score = min(E_mag_thresholded_trials(:,clusterMask), [], 2);
    score = median(E_mag_thresholded_trials(:,clusterMask), 2); 
    [~, whichTriangle] = max(score);
    clusters.(clusters.Names{i}).bestTriangle = whichTriangle;
    max_loc = mean(coordinates(:,triangles(whichTriangle,:)), 2);

    ts = trisurf(triangles, coordinates(1,:), coordinates(2,:), coordinates(3,:), score);
    ts.EdgeColor = 'none';
    hold on
    scatter3(v.p1(clusterMask), v.p2(clusterMask), v.p3(clusterMask), 8, baseColors(i,:), "filled")
    scatter3(max_loc(1), max_loc(2), max_loc(3), 30, 'kv')
    scatter3(max_loc_R2(1), max_loc_R2(2), max_loc_R2(3), 20, 'wo', "filled")

    ylim(1.1.*[min([min(coordinates(2,:)) min(v.p2)]) max([max(coordinates(2,:)) max(v.p2)])])
    cube_axis(gca)
    colormap(spring_colormap())
    cb = colorbar;
    cb.Label.String = 'smallest |E| for any of the trials in this cluster';
    title(sprintf('Cluster %d', i), Color=baseColors(i,:))
    axis square; grid off; set(gca, 'Color', 'none'); view(0, 90); box on
end




exportgraphics(fig, [OUT 'cluster_locations.png'])



%%
% Weise method has sigmoid4:
% Literally: y = y_0 + \\frac{amp - y_0}{1+e^{-r(x-x_0)}}
% i.e.:      y = y_0 + (amp - y_0)/(1+exp(-r*(x - x_0)))
% So first question: do i find the same result R² for the best triangle?
sigmoid4 = 'a + (b - a)/(1+exp(-c*(x - d)))';
x = E_mag(bestR2Triangle,:)';
y = responses.(response);
[f1, gof, ~] = fit(x,y,sigmoid4,'Start', [quantile(y, 0.25) quantile(y, 0.75) range(y)/range(x) mean(x)]);


figure;
plot(x, y, 'k.', MarkerSize=1)
hold on
plot(f1)

R2_in_matlab = gof.rsquare;
fprintf('   R² in Matlab   = %0.4f\nvs R² from python = %.4f\n\n', R2_in_matlab, bestR2)
%%

% Multivariate fit:
Edata = [];
for i = clusters.IDs
    Edata = [Edata E_mag(clusters.(clusters.Names{i}).bestTriangle,:)'];
end
y = responses.(response);
A0 = quantile(y,0.75);
r0 = range(y) ./ range(Edata, 1);
x00 = mean(Edata, 1);
p0 = [quantile(y,0.25) repmat([A0 nan nan], 1, n_clusters)];
p0(4:3:end) = x00;
p0(3:3:end) = r0;

%options = optimoptions('lsqcurvefit', Algorithm='interior-point', MaxFunctionEvaluations = 5e3);
options = optimoptions('lsqcurvefit', Algorithm='levenberg-marquardt', MaxFunctionEvaluations = 5e3, ScaleProblem='jacobian');
% Try to get levenberg-marquardt to work!
[bestparams_mv, resnorm, ~, exitflag] = lsqcurvefit(@mvarsigmoid, p0, Edata, y, [], [], options);
R2_mv = 1-(resnorm / (length(y) * var(y)));
fprintf('   R² univariate   = %0.4f\n   R² multivariate = %0.4f\n\n', bestR2, R2_mv)
% TODO: continue here!



if max([bestR2 R2_in_matlab]) < R2_mv
    fprintf("! Evidence for nonsingularity\n")
end







%% Inspect data directly

fig = figure(Position=[1200 50 1000 1000]);
tiledlayout(3,3)
nexttile(2)
plot(E_mag(bestR2Triangle, directBelowThresholdMask), responses.(response)(directBelowThresholdMask), 'k.', MarkerSize=1)
hold on
for i = clusters.IDs
    plot(E_mag_thresholded_trials(bestR2Triangle, clusters.labels == i), R_thresholded_trials(clusters.labels == i), '.', Color=baseColors(i,:));
end
axis square; box on; set(gca, 'Color', 'none')
xlabel('|E|'); ylabel(response, Interpreter='none'); title('At global best R² triangle')

for i = clusters.IDs
    nexttile(i + 3)
    whichTriangle = clusters.(clusters.Names{i}).bestTriangle;
    plot(E_mag(whichTriangle, directBelowThresholdMask), responses.(response)(directBelowThresholdMask), 'k.', MarkerSize=1)
    hold on
    for j = clusters.IDs
        marker = '.';
        if i == j
            marker = 'o';
        end
        plot(E_mag_thresholded_trials(whichTriangle, clusters.labels == j), R_thresholded_trials(clusters.labels == j), marker, Color=baseColors(j,:));
    end
    x = linspace(min(E_mag(whichTriangle, :)), max(E_mag(whichTriangle, :)), 100)';
    X = zeros(length(x), n_clusters);
    X(:, i) = x;
    y = mvarsigmoid(bestparams_mv, X);
    plot(x,y,'-',Color=baseColors(i,:))
    axis square; box on; set(gca, 'Color', 'none')
    xlabel('|E|'); ylabel(response, Interpreter='none'); title(sprintf('For cluster %d', i), Color=baseColors(i,:))
end


exportgraphics(fig, [OUT 'IO-scatter.png'])




%% Multivariate fit:
if false
% Toy example:
[X,Y] = meshgrid(linspace(0,3,100), linspace(0,3,100));
E = [X(:) Y(:)];
params = [-1 2 2 1.5 1 10 1.5];
y = mvarsigmoid(params, E);
y = y + randn(size(y));
plot3(E(:,1), E(:,2), y, 'k.', MarkerSize=1)

bestparams = lsqcurvefit(@mvarsigmoid, params+randn(size(params)), E, y);

hold on
predicted = mvarsigmoid(bestparams, E);
s = surf(X, Y, reshape(predicted, [100, 100]), EdgeColor='none');

end








%% mRMR?

if false
[idx, mrmrRank] = fsrmrmr(E_mag', responses.(response), Verbose=2);
save([OUT 'mRMR.mat'], "mrmrRank", "idx")
%%
fig = figure;
ts = trisurf(triangles, coordinates(1,:), coordinates(2,:), coordinates(3,:), mrmrRank);
ts.EdgeColor = 'none';
hold on

for i = 1:5
    nodes = triangles(idx(i),:);
    x = mean(coordinates(1,nodes));
    y = mean(coordinates(2,nodes));
    z = mean(coordinates(3,nodes));
    scatter3(x,y,z, 40, [1 0 0],'o', "filled")
end


cube_axis(gca)
colormap(honey_colormap())
cb = colorbar;
cb.Label.String = 'mRMR rank';
axis square; grid off; set(gca, 'Color', 'none'); view(0, 90); box on

%exportgraphics(fig, [OUT 'mRMR.png'])

%% Multivariate fit:
y = responses.(response);
Edata = [];
nMrmrPredictors = 5;
for i = 1:nMrmrPredictors
    Edata = [Edata E_mag(idx(i),:)'];
end
A0 = quantile(y,0.9);
r0 = 1e5;
x00 = mean(Edata, 1);
p0 = [0 repmat([A0 r0 nan], 1, nMrmrPredictors)];
p0(4:3:end) = x00;

bestparams = lsqcurvefit(@mvarsigmoid, p0, Edata, y);
R2_mv = 1-(mean((y - mvarsigmoid(bestparams, Edata)).^2) / var(y));
fprintf('   R² univariate   = %0.4f\n   R² multivariate = %0.4f\n\n', R2_in_matlab, R2_mv)

% mRMR also yields only a marginal improvement from R"(uv)=0.2433 to R²(mv)
% =~ 0.25 (at most 0.01 increase)
% We can probably take this to indicate that there is likely no secondary
% source!

end




















function ydata = mvarsigmoid(params, xdata)
% params will be 1+(n_clusters*3), as: y0, A_1, r_1, x0_1, A_2, r_2, x0_2,
% A_3, etc
% xdata is (ntrials, n_clusters)
y0 = params(1);
A = params(2:3:end);
r = params(3:3:end);
x0 = params(4:3:end);
ydata = y0 + sum((A - y0) ./ (1+exp(-r.*(xdata - x0))), 2);
end






