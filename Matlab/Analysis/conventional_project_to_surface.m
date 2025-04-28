addpath(genpath('//wsl.localhost/Ubuntu-22.04/usr/local/freesurfer/matlab'))
addpath(genpath('B:/Projects/2023-01 HighDef/libraries/visualization'))
addpath(genpath('B:/Projects/2023-01 HighDef/libraries/vetter'))

STORAGE = "D:/HighDef-operate/HighDef";


fs_subject_dir = '//wsl.localhost/Ubuntu-22.04/usr/local/freesurfer/subjects';


DO_CUBIC = true;



ResultTable = [];

subjects = ["sub-001", "sub-002", "sub-003", "sub-004", "sub-006", "sub-007", "sub-008", "sub-009", "sub-011", "sub-012", "sub-013", "sub-014", "sub-015", "sub-016", "sub-017", "sub-019", "sub-022", "sub-023"];
leftHanded = ["sub-006", "sub-016", "sub-017"];

for sub_id = subjects
    %% Load and merge freesurfer hemispheres
    fs_surface = "inflated";

    for hemisphere = ["l", "r"]
        [vertex_coords, faces] = read_surf(sprintf("%s/fsaverage/surf/%sh.%s", fs_subject_dir, hemisphere, fs_surface));
        separate_cortices.(hemisphere).vertex_coords = vertex_coords;
        if fs_surface == "inflated"
            if hemisphere == "l"
                x = separate_cortices.(hemisphere).vertex_coords(:, 1);
                separate_cortices.(hemisphere).vertex_coords(:, 1) = x - quantile(x, 0.95);
            else
                x = separate_cortices.(hemisphere).vertex_coords(:, 1);
                separate_cortices.(hemisphere).vertex_coords(:, 1) = x - quantile(x, 0.05);
            end
        end

        separate_cortices.(hemisphere).triangles = faces+1;
        separate_cortices.(hemisphere).curvature = read_curv(sprintf('%s/fsaverage/surf/%sh.curv', fs_subject_dir, hemisphere));
    end


    fs_coords = [separate_cortices.l.vertex_coords; separate_cortices.r.vertex_coords];
    % fs_coords(1:N_left, :) will be the left hemispheric vertex coordinates
    % --> thus, the left-hemispheric triangle indices are still correct
    % But fs_coords(N_left+1:N_left+N_right, :) will be the right hemispheric
    % coordinates
    % --> add N_left to right triangle indices
    n_left = size(separate_cortices.l.vertex_coords, 1);
    fs_triangles = [separate_cortices.l.triangles; separate_cortices.r.triangles + n_left];

    fs_curvature = [separate_cortices.l.curvature; separate_cortices.r.curvature];

    % Need to add cerebellum
    %custom_cerebellum
    res = 100;
    [el,az] = meshgrid(linspace(0, pi, res), linspace(-pi, pi, res));
    T = delaunay(el, az);
    el = el(:);
    az = az(:);
    r = 30; a = 0.5; b = 0.8; c = 1;
    x = r .* sin(el) .* cos(az) ./ sqrt(a);
    y = r .* sin(el) .* sin(az) ./ sqrt(b);
    z = r .* cos(el) ./ sqrt(c);

    y = y - 55; z = z - 60;

    fs_triangles = [fs_triangles; T + size(fs_coords, 1)];
    fs_curvature = [fs_curvature; ones(size(x, 1), 1)];
    fs_coords = [fs_coords; x y z];


    %% Load individual brain:

    hdf5_file = sprintf("%s/%s/mesh/m2m_%s/%s.hdf5", STORAGE, sub_id, sub_id, sub_id);

    sub_coords = h5read(hdf5_file, '/mesh/nodes/node_coord')';
    sub_all_triangles = h5read(hdf5_file, '/mesh/elm/triangle_number_list')' + 1;
    sub_triangle_tissue_types = h5read(hdf5_file, '/mesh/elm/tri_tissue_type');

    sub_triangles = sub_all_triangles(sub_triangle_tissue_types == 1002,:);


    %%
    all_nodes_used = unique(sub_triangles(:));         % some list of indices [i1 i2 i3 ... iN]
    sub_coords = sub_coords(all_nodes_used,:);         % The first row will correspond to i1, etc
    % We now need to re-label the triangles -- i1 is outdated, it should be 1,
    % i2 -> 2, i3 -> 3, ... iN -> N
    all_nodes_used_indices = nan(1, max(all_nodes_used)); % [1 2 3 ... N]
    all_nodes_used_indices(all_nodes_used) = 1:length(all_nodes_used);
    sub_triangles = all_nodes_used_indices(sub_triangles);




    figure;
    subplot(1,2,1)
    ts = trisurf(fs_triangles, fs_coords(:, 1), fs_coords(:, 2), fs_coords(:, 3));
    ts.FaceColor=[0.5 0.5 0.5];
    ts.FaceAlpha = 0.5;
    ts.EdgeColor="none";

    hold on

    ts = trisurf(sub_triangles, sub_coords(:, 1), sub_coords(:, 2), sub_coords(:, 3));
    ts.FaceColor=[0.8 1 0.5];
    ts.FaceAlpha = 0.5;
    ts.EdgeColor="none";

    l = max([range(xlim) range(ylim) range(zlim)]);
    xlim(mean(xlim) + [-l l]./2); ylim(mean(ylim) + [-l l]./2); zlim(mean(zlim) + [-l l]./2)
    axis square
    view(-90, 0)


    %% Align individual and avg brains

    % prepare: Center and scale
    fs_coords = fs_coords - mean(fs_coords, 1);
    sub_coords = sub_coords - mean(sub_coords, 1);

    scale_fs  = mean(vecnorm(fs_coords'));
    scale_sub = mean(vecnorm(sub_coords'));

    scale = scale_fs / scale_sub;

    sub_coords = sub_coords .* scale;


    ptCloud_fs  = pointCloud(fs_coords);
    ptCloud_sub = pointCloud(sub_coords);

    % Use pcregistericp to align point clouds
    [tform, movingReg] = pcregistericp(ptCloud_sub, ptCloud_fs);


    subplot(1,2,2)
    ts = trisurf(fs_triangles, fs_coords(:, 1), fs_coords(:, 2), fs_coords(:, 3));
    ts.FaceColor=[0.3 0.3 0.3];
    ts.FaceAlpha = 0.1;
    ts.EdgeColor="none";

    hold on

    sub_coords_transformed = [sub_coords ones(size(sub_coords, 1), 1)] * tform.A';
    sub_coords_transformed = sub_coords_transformed(:, 1:3);
    ts = trisurf(sub_triangles, sub_coords_transformed(:, 1), sub_coords_transformed(:, 2), sub_coords_transformed(:, 3));
    ts.FaceColor=[0.8 1 0.5];
    ts.FaceAlpha = 0.5;
    ts.EdgeColor="none";

    l = max([range(xlim) range(ylim) range(zlim)]);
    xlim(mean(xlim) + [-l l]./2); ylim(mean(ylim) + [-l l]./2); zlim(mean(zlim) + [-l l]./2)
    axis square
    view(-90, 0)


    %%


    fig = figure(Position=[100 50 1700 800]);
    ax1 = nexttile();
    hold on
    axis tight
    ax2 = nexttile();
    hold on
    axis tight
    axs = [ax1, ax2];

    fig_2d = figure(Position=[150 100 800 700]);
    axs2d.dominant.CsE_FDI_in_uV    = nexttile(); axis tight
    axs2d.nondominant.CsE_FDI_in_uV = nexttile(); axis tight
    axs2d.dominant.SIHIscore_FDI    = nexttile(); axis tight
    axs2d.nondominant.SIHIscore_FDI = nexttile(); axis tight

    set(0, 'CurrentFigure', fig)

    response_names = ["CsE_FDI_in_uV", "SIHIscore_FDI"];
    for iResponse = 1:2
        response_name = response_names(iResponse);

        for hemisphere = ["dominant", "nondominant"]

            if ~ismember(sub_id, leftHanded)
                if strcmpi(hemisphere, "dominant")
                    exp_id = "L2R";
                else
                    exp_id = "R2L";
                end
            else
                if strcmpi(hemisphere, "dominant")
                    exp_id = "R2L";
                else
                    exp_id = "L2R";
                end
            end

            axes(axs(iResponse));


            data = readtable(sprintf("%s/%s/%s_%s_raw.csv", STORAGE, sub_id, sub_id, exp_id));

            rejected_by_preinnervation = data.preinnervation_FDI_in_uV > 50;
            only_high_SI_trials = data.Intensity_percentMSO == max(data.Intensity_percentMSO);
            accepted_trials = only_high_SI_trials & ~rejected_by_preinnervation;
            %accepted_trials

            locations = [data.p1 data.p2 data.p3 ones(size(data,1), 1)];
            locations = locations(accepted_trials, :);
            transformed_locations = (tform.A * locations')';
            transformed_locations = transformed_locations(:,1:3);

            [Idx, ~] = knnsearch(fs_coords, transformed_locations);
            projected = fs_coords(Idx,:);


            %% [A]: Plotting the data and COM on the average brain
            clamped_curv = ones(size(fs_curvature));
            clamped_curv(fs_curvature < 0) = -1;

            response = data.(response_name);
            response = response(accepted_trials);
            if startsWith(response_name, "CsE")
                [COM, ~, MedOM]    = getWeightedMetrics(projected, response ./ sum(response));
                colors = warm_colormap();
            else
                [COM, ~, MedOM]    = getWeightedMetrics(projected, exp(response) ./ sum(exp(response)));
                colors = cold_colormap();
            end
            projectedCOM = fs_coords(knnsearch(fs_coords, COM),:);

            ts = trisurf(fs_triangles, fs_coords(:,1), fs_coords(:,2), fs_coords(:,3), clamped_curv, EdgeColor="none");
            colormap(0.8.*[0.355 0.35 0.37; 0.215 0.21 0.23]);
            hold on


            normalized_response = (response - min(response)) / range(response);
            color_indices = 1:size(colors,1);
            colors = interp1(color_indices, colors, normalized_response * (size(colors, 1) - 1) + 1);
            scatter3(projected(:,1), projected(:,2), projected(:,3), 10, colors, "filled")

            plot3(projectedCOM(1), projectedCOM(2), projectedCOM(3)+3, "w+", LineWidth=3, MarkerSize=10)
            otheraxis = axs(3-iResponse);
            plot3(otheraxis, projectedCOM(1), projectedCOM(2), projectedCOM(3)+3, "w+", LineWidth=1, MarkerSize=10)


            xlabel("x")
            ylabel("y")
            zlabel("z")

            l = max([range(xlim) range(ylim) range(zlim)]);
            xlim(mean(xlim) + [-l l]./2); ylim(mean(ylim) + [-l l]./2); zlim(mean(zlim) + [-l l]./2)
            axis square

            title(sprintf("%s %s", sub_id, response_name), Interpreter="none")

            view(2)
            set(gca, Color="none")
            grid off; box on






            %% [B]: Flatten projected to 2D, interpolate, find highest peak and area

            axes(axs2d.(hemisphere).(response_name))

            projected_2D = project_to_standard(array2table(projected, VariableNames=["p1", "p2", "p3"]), extractBefore(exp_id, "2"));
            x = projected_2D(:,2);
            y = projected_2D(:,1);

            if startsWith(response_name, "CsE")
                v = response;
            else
                v = exp(-response);
            end

            if DO_CUBIC
                gridlength = 300;
                gridSide = linspace(min([min(x) min(y)]), max([max(x) max(y)]), gridlength);
                [xq, yq] = meshgrid(gridSide, gridSide);
                vq = griddata(x,y,v,xq,yq, "cubic");
            else
                gridLength = 150;
                sd = 4;
                maxDistanceOfPeak = 4*sd;
                kernel = @(x,y) (1/(2*pi*sd^2)).*exp(-(x.^2 + y.^2)./(2*sd^2));
                valueThreshold = 20*kernel(1,0);
                [vq, xq, yq] = interpolateResponseMap([x y], v, gridLength, sd, kernel, valueThreshold, false);
            end

            if startsWith(response_name, "CsE")
                [~, best_index] = max(vq(:));
            else
                [~, best_index] = min(vq(:));
            end
            xs = xq(:); x_peak = xs(best_index);
            ys = yq(:); y_peak = ys(best_index);
            aboveMap = 1.01*max(vq(:));

            s=surf(xq, yq, vq, EdgeColor="none");
            hold on
            plot3(x_peak, y_peak, aboveMap, "ko")

            if startsWith(response_name, "SIHI")
                clim([min(clim), 1])
                colormap(axs2d.(hemisphere).(response_name), flipud(cold_colormap(10)))
            else
                colormap(axs2d.(hemisphere).(response_name), warm_colormap(10))
            end

            axis square; grid off; box on; set(axs2d.(hemisphere).(response_name), Color="none")
            l = max([range(xlim) range(ylim)]);
            xlim(mean(xlim) + [-l l]./2); ylim(mean(ylim) + [-l l]./2);
            view(2)

            colorbar
            title(sprintf("%s %s %s", sub_id, hemisphere, response_name), Interpreter="none")

            %% Compute area:
            grid_cell_dx = mean(diff(xq'), "all"); % in mm
            grid_cell_dy = mean(diff(yq), "all");  % in mm
            grid_cell_area = grid_cell_dx * grid_cell_dy; % in mm²

            if startsWith(response_name, "CsE")
                area_superthreshold = grid_cell_area * sum(vq(:) > 50);
                area_strong = grid_cell_area * sum(vq(:) > 1000);
            else
                area_superthreshold = grid_cell_area * sum(vq(:) < 1);
                area_strong = grid_cell_area * sum(vq(:) < 2/3);
            end

            fprintf("%s %s %s\n Area of response passing threshold: %3.3f mm²\n Area of strong response: %3.3f mm²\n\n", sub_id, hemisphere, response_name, area_superthreshold, area_strong)


            %% Collect in result table

            newrow = table(sub_id, response_name, hemisphere, projectedCOM(1), projectedCOM(2), projectedCOM(3), x_peak, y_peak, area_superthreshold, area_strong, VariableNames=["Subject", "Response", "Hemisphere", "x_COM", "y_COM", "z_COM", "x_peak", "y_peak", "area_where_response_passes_threshold_in_mm2", "area_where_response_is strong_in_mm2"]);
            ResultTable = [ResultTable; newrow];
        end

    end


    exportgraphics(fig, sprintf("B:/Projects/2023-01 HighDef/Results/Coil-space/avg-brain-%s.png", sub_id), ...
        'BackgroundColor', 'none', 'ContentType', 'image', "Resolution", 400)

    if DO_CUBIC
        interpolation_type = "cubic";
    else
        interpolation_type = "kernel";
    end

    exportgraphics(fig_2d, sprintf("B:/Projects/2023-01 HighDef/Results/Coil-space/maps-%s-fs_avg-%s.png", interpolation_type, sub_id), ...
        'BackgroundColor', 'none', 'ContentType', 'image', "Resolution", 400)
end





%%


subjects = unique(ResultTable.Subject);
nSubjects = length(subjects);
areas.dominant.hot     = nan(nSubjects,1);
areas.nondominant.hot  = nan(nSubjects,1);
areas.dominant.cold    = nan(nSubjects,1);
areas.nondominant.cold = nan(nSubjects,1);

spots.dominant.hot     = nan(nSubjects, 3);
spots.dominant.cold    = nan(nSubjects, 3);
spots.nondominant.hot  = nan(nSubjects, 3);
spots.nondominant.cold = nan(nSubjects, 3);

offsets.dominant     = nan(nSubjects,3);
offsets.nondominant  = nan(nSubjects,3);

distances.dominant     = nan(nSubjects,1);
distances.nondominant  = nan(nSubjects,1);

names.hot = "CsE_FDI_in_uV";
names.cold = "SIHIscore_FDI";

for iSubject = 1:nSubjects
    for hemisphere = ["dominant", "nondominant"]
        for spot = ["hot", "cold"]
            mask = ResultTable.Subject == subjects(iSubject) & ResultTable.Response == names.(spot) & ResultTable.Hemisphere == hemisphere;
            row = ResultTable(mask,:);
            areas.(hemisphere).(spot)(iSubject) = row.area_where_response_passes_threshold_in_mm2;
            spots.(hemisphere).(spot)(iSubject,:) = [row.x_COM row.y_COM row.z_COM];
        end
        mask_hot  = ResultTable.Subject == subjects(iSubject) & ResultTable.Response == names.hot  & ResultTable.Hemisphere == hemisphere;
        mask_cold = ResultTable.Subject == subjects(iSubject) & ResultTable.Response == names.cold & ResultTable.Hemisphere == hemisphere;
        offsets.(hemisphere)(iSubject,1) = ResultTable(mask_cold,:).x_COM - ResultTable(mask_hot,:).x_COM;
        if ismember(subjects(iSubject), leftHanded)
            offsets.(hemisphere)(iSubject,1) = -offsets.(hemisphere)(iSubject,1);
        end
        offsets.(hemisphere)(iSubject,2) = ResultTable(mask_cold,:).y_COM - ResultTable(mask_hot,:).y_COM;
        offsets.(hemisphere)(iSubject,3) = ResultTable(mask_cold,:).z_COM - ResultTable(mask_hot,:).z_COM;

        distances.(hemisphere)(iSubject) = norm(offsets.(hemisphere)(iSubject,:));
    end
end



%%
fig_measures = figure;
subplot(2,2,1)
plot([zeros(nSubjects, 1) offsets.dominant(:,1)]', [zeros(nSubjects, 1) offsets.dominant(:,2)]', "k-")
hold on
plot(offsets.dominant(:,1), offsets.dominant(:,2), "bv", MarkerFaceColor="b")
plot([0 mean(offsets.dominant(:,1))], [0 mean(offsets.dominant(:,2))], "r-", LineWidth=2)
title("Dominant hemisphere (CoG)")
axis square; set(gca, Color="none")
%l = max([range(xlim) range(ylim)]); xlim(mean(xlim) + 0.5.*[-l l]); ylim(mean(ylim) + 0.5.*[-l l])
xlim([-7 7]); ylim([-7 7])
plot([-8 8], [-8 8], 'k:')

subplot(2,2,2)
plot([zeros(nSubjects, 1) offsets.nondominant(:,1)]', [zeros(nSubjects, 1) offsets.nondominant(:,2)]', "k-")
hold on
plot(offsets.nondominant(:,1), offsets.nondominant(:,2), "bv", MarkerFaceColor="b")
plot([0 mean(offsets.nondominant(:,1))], [0 mean(offsets.nondominant(:,2))], "r-", LineWidth=2)
title("Nondominant hemisphere (CoG)")
axis square; set(gca, Color="none")
%l = max([range(xlim) range(ylim)]); xlim(mean(xlim) + 0.5.*[-l l]); ylim(mean(ylim) + 0.5.*[-l l])
xlim([-7 7]); ylim([-7 7])
plot([-8 8], [8 -8], 'k:')

subplot(2,2,3)
x = areas.dominant.hot;
y = areas.dominant.cold;
plot(x, y, 'kx')
hold on
m = min([x; y]);
M = max([x; y]);
plot([m M], [m M], 'k-')
axis tight; axis square; set(gca, Color="none")
xlabel("Area hotspot [mm²]")
ylabel("Area coldspot [mm²]")
title("Dominant hemisphere")

subplot(2,2,4)
x = areas.nondominant.hot;
y = areas.nondominant.cold;
plot(x, y, 'kx')
hold on
m = min([x; y]);
M = max([x; y]);
plot([m M], [m M], 'k-')
axis tight; axis square; set(gca, Color="none")
xlabel("Area hotspot [mm²]")
ylabel("Area coldspot [mm²]")
title("Nondominant hemisphere")


p_dominant = signrank(areas.dominant.hot, areas.dominant.cold); % Test if hotspot areas in median is less than coldspot area
fprintf("Dominant hemisphere:      Wilcoxon signed rank test whether hotspot has smaller area: p = %3.3f\n", p_dominant)
fprintf("                          avg.  hotspot area = %3.3f cm²\n", mean(areas.dominant.hot) / 100)
fprintf("                          avg. coldspot area = %3.3f cm²\n", mean(areas.dominant.cold) / 100)

p_nondominant = signrank(areas.nondominant.hot, areas.nondominant.cold); % Test if hotspot areas in median is less than coldspot area
fprintf("Nondominant hemisphere:   Wilcoxon signed rank test whether hotspot has smaller area: p = %3.3f\n", p_nondominant)
fprintf("                          avg.  hotspot area = %3.3f cm²\n", mean(areas.nondominant.hot) / 100)
fprintf("                          avg. coldspot area = %3.3f cm²\n", mean(areas.nondominant.cold) / 100)


exportgraphics(fig_measures, 'B:/Projects/2023-01 HighDef/Results/Coil-space/projection-approach_spot-comparison.pdf')

%% Compute avg vector:
for hemisphere = ["dominant", "nondominant"]
    mx = mean(offsets.(hemisphere)(:,1));
    my = mean(offsets.(hemisphere)(:,2));
    sx = std(offsets.(hemisphere)(:,1));
    sy = std(offsets.(hemisphere)(:,2));

    C = cov(offsets.(hemisphere)(:,1), offsets.(hemisphere)(:,2));
    covxy = C(1,2);

    l = sqrt(mx^2 + my^2);
    term_x = mx*sx;
    term_y = my*sy;
    term_cov = 2*mx*my*covxy;
    sl = sqrt(term_x^2 + term_y^2 + term_cov) / l;
    fprintf("%s hemisphere:\t shifted by %3.3f ± %3.3f mm along mean direction\n", hemisphere, l, sl)

    if hemisphere == "dominant"
        compare_along = [-1 -1 -1] ./ sqrt(3); % downward posterolateral (more along the surface)
    else
        compare_along = [1 -1 -1] ./ sqrt(3); % downward posterolateral (more along the surface)
    end

    scoring = offsets.(hemisphere) * compare_along';
    p_projected = signrank(scoring, 0, "tail", "right");
    fprintf("%s hemisphere:\t Shifted by %3.3f ± %3.3f mm posterolaterally     p = %3.3f\n", hemisphere, mean(scoring), std(scoring), p_projected)
end




%%
fig_avg_brain = figure(Position=[50 50 900 900]);
ts = trisurf(fs_triangles, fs_coords(:,1), fs_coords(:,2), fs_coords(:,3), clamped_curv, EdgeColor="none");
colormap(0.8.*[0.355 0.35 0.37; 0.215 0.21 0.23]);
hold on
axis square; axis tight

z_offset = 10;
plot3([spots.dominant.hot(:,1) spots.dominant.cold(:,1)]', ...
    [spots.dominant.hot(:,2) spots.dominant.cold(:,2)]', ...
    [spots.dominant.hot(:,3) spots.dominant.cold(:,3)]' + z_offset - 1, "w.-", LineWidth=2)
plot3([spots.nondominant.hot(:,1) spots.nondominant.cold(:,1)]', ...
    [spots.nondominant.hot(:,2) spots.nondominant.cold(:,2)]', ...
    [spots.nondominant.hot(:,3) spots.nondominant.cold(:,3)]' + z_offset - 1, "w.-", LineWidth=2)

plot3(spots.dominant.hot(:,1), spots.dominant.hot(:,2), spots.dominant.hot(:,3) + z_offset, 'r^', MarkerFaceColor="r")
plot3(spots.nondominant.hot(:,1), spots.nondominant.hot(:,2), spots.nondominant.hot(:,3) + z_offset, 'r^')

cold_color = [0, 150, 255]./255;

plot3(spots.dominant.cold(:,1), spots.dominant.cold(:,2), spots.dominant.cold(:,3) + z_offset, 'v', MarkerFaceColor=cold_color, MarkerEdgeColor=cold_color)
plot3(spots.nondominant.cold(:,1), spots.nondominant.cold(:,2), spots.nondominant.cold(:,3) + z_offset, 'v', MarkerEdgeColor=cold_color)


xlabel("x")
ylabel("y")
zlabel("z")

view(2)
set(gca, Color="none")
grid off; box on

l = max([range(xlim) range(ylim) range(zlim)]);
xlim(mean(xlim) + [-l l]./2); ylim(mean(ylim) + [-l l]./2); zlim(max(zlim) + [-l 0])

exportgraphics(fig_avg_brain, "B:/Projects/2023-01 HighDef/Results/Coil-space/cogs_on_fs_avg.png", ...
        'BackgroundColor', 'none', 'ContentType', 'image', "Resolution", 400)


l = l * 0.6;
xlim(mean(xlim) + [-l l]./2); ylim(mean(ylim) + [-l l]./2); zlim(max(zlim) + [-l 0])

exportgraphics(fig_avg_brain, "B:/Projects/2023-01 HighDef/Results/Coil-space/cogs_on_fs_avg-zoomed.png", ...
        'BackgroundColor', 'none', 'ContentType', 'image', "Resolution", 400)

