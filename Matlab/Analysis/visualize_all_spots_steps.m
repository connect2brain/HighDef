
for ID = ["sub-001", "sub-002", "sub-003", "sub-004", "sub-006", "sub-007", "sub-008", "sub-009", "sub-011", "sub-012", "sub-013", "sub-014"]
    for exp_id = ["R2L", "L2R"]

        %ID = 'sub-001';
        %exp_id = 'R2L';
        hemisphere = extractBefore(exp_id, "2");
        mesh_id = 'mesh0';

        roi = sprintf('midlayer_%s', lower(hemisphere));
        ROOT = sprintf('//wsl.localhost/Ubuntu-22.04/home/bnplab-admin/TMS_localization/HighDef/%s/results/exp_map-%s', ID, exp_id);
        R2_PATH = sprintf('%s/r2/mesh_%s/roi_%s', ROOT, mesh_id, roi);

        templates = [];
        templates.hot = 'CsE_%s_in_uV';
        templates.cold = 'SIHIscore_%s';

        spots = struct();
        spots.hot = [];
        spots.cold = [];

        S1hotspots = struct();


        for muscle = ["ADM", "FDI", "APB"]
            for spot = ["hot", "cold"]
                figure;
                [~, max_loc, bestMatch, maxR2, corners] = plot_R2_map(sprintf('%s/%s/sigmoid4', R2_PATH, sprintf(templates.(spot), muscle)), true, 'mag');
                close all;
                spots.(spot).(muscle).max_loc = max_loc';
                spots.(spot).(muscle).max_ind = bestMatch;
                spots.(spot).(muscle).max_R2 = maxR2;
                spots.(spot).(muscle).normal = double(cross(corners(:,2)-corners(:,1), corners(:,3)-corners(:,1))');
                spots.(spot).(muscle).normal = spots.(spot).(muscle).normal ./ norm(spots.(spot).(muscle).normal);
            end

            figure;
            S1_path = sprintf('//wsl.localhost/Ubuntu-22.04/home/bnplab-admin/TMS_localization/HighDef/%s/results/exp_map-%s/r2/mesh_%s/roi_%s/%s/sigmoid4', ID, hemisphere, mesh_id, roi, sprintf(templates.("hot"), muscle));
            [~, max_loc, bestMatch, maxR2, corners] = plot_R2_map(S1_path, true, 'mag');
            close all;
            S1hotspots.(muscle).max_loc = max_loc';
            S1hotspots.(muscle).max_ind = bestMatch;
            S1hotspots.(muscle).max_R2 = maxR2;
            S1hotspots.(muscle).normal = double(cross(corners(:,2)-corners(:,1), corners(:,3)-corners(:,1))');
            S1hotspots.(muscle).normal = spots.(spot).(muscle).normal ./ norm(spots.(spot).(muscle).normal);
        end



        %% Compute avg normal and transform:


        rotate_to_face_normal = eye(3);


        %% Plotting


        colors = [];                    markers = [];
        colors.cold.ADM = "k";          markers.cold.ADM = "";
        colors.cold.FDI = "b";    markers.cold.FDI = "o";
        colors.cold.APB = "k";          markers.cold.APB = "";
        colors.hot.ADM  = "k";          markers.hot.ADM  = "^";
        colors.hot.FDI  = "r";          markers.hot.FDI  = "o";
        colors.hot.APB  = "k";          markers.hot.APB  = "square";

        geo_file = sprintf('//wsl.localhost/Ubuntu-22.04/home/bnplab-admin/TMS_localization/HighDef/%s/mesh/roi/%s/geo.hdf5', ID, roi);
        coordinates = h5read(geo_file, '/mesh/nodes/node_coord');
        triangles = h5read(geo_file, '/mesh/elm/triangle_number_list')' + 1;

        raw_excentricity = sqrt(sum((coordinates - [0; 0; 10]) .^ 2, 1));
        coordinates = rotate_to_face_normal * coordinates;

        % Pick top 5 % of triangles from hotspot FDI, plot in red; pick top 5 % of
        % triangles from coldspot FDI, color blue; color intersection purple.
        hotspot_data = h5read(sprintf('//wsl.localhost/Ubuntu-22.04/home/bnplab-admin/TMS_localization/HighDef/%s/results/exp_map-%s/r2/mesh_%s/roi_%s/CsE_FDI_in_uV/sigmoid4/r2_roi_data.hdf5', ID, exp_id, mesh_id, roi), '/data/tris/c_E_mag');
        coldspot_data = h5read(sprintf('//wsl.localhost/Ubuntu-22.04/home/bnplab-admin/TMS_localization/HighDef/%s/results/exp_map-%s/r2/mesh_%s/roi_%s/SIHIscore_FDI/sigmoid4/r2_roi_data.hdf5', ID, exp_id, mesh_id, roi), '/data/tris/c_E_mag');

        hotspot_data(isnan(hotspot_data)) = 0;
        hotspot_data = reject_R2_outliers_if_needed(hotspot_data);
        coldspot_data(isnan(coldspot_data)) = 0;
        coldspot_data = reject_R2_outliers_if_needed(coldspot_data);


        best_R2_hotspot = max(hotspot_data);
        best_R2_coldspot = max(coldspot_data);

        height_threshold = 0.5;
        triangles_hot = hotspot_data > height_threshold * best_R2_hotspot;
        triangles_cold = coldspot_data > height_threshold * best_R2_coldspot;

        fprintf('hotspot:  R² = %0.4f    FWHM = %0.4f\n', best_R2_hotspot, sum(triangles_hot))
        fprintf('coldspot: R² = %0.4f    FWHM = %0.4f\n\n', best_R2_coldspot, sum(triangles_cold))

        triangles_both = triangles_hot & triangles_cold;
        triangles_only_hot = triangles_hot & ~triangles_cold;
        triangles_only_cold = triangles_cold & ~triangles_hot;


        %% Step 1: Only brain

        fig1 = figure(Position=[50 50 600 600]);
        hold on
        ts = trisurf(triangles, coordinates(1,:), coordinates(2,:), coordinates(3,:), mean(raw_excentricity(triangles), 2));
        ts.EdgeColor = 'none';
        ts.FaceAlpha = 1;
        %ts.FaceColor = "w";
        cr = max([range(xlim) range(ylim) range(zlim)]);
        xlim(mean(xlim) + [-cr cr]./2); ylim(mean(ylim) + [-cr cr]./2); zlim(mean(zlim) + [-cr cr]./2)
        clim([min(clim) 1.05*max(clim)])
        set(gca(), 'Color', 'none')
        colormap("gray")
        grid on; grid minor; box on; axis square; view(2)
        xticks(-100:5:100); yticks(-100:5:100);
        labels = arrayfun(@(v) sprintf("%d", v), xticks); labels(mod(xticks, 10) ~= 0) = "";
        xticklabels(labels); yticklabels(labels)

        where_anterior = rotate_to_face_normal * [0 10 0]';

        %cr = max([range(xlim) range(ylim) range(zlim)]);
        shrink = 0.6;
        crs = shrink*0.5.*[-cr cr];
        xlim(mean([spots.hot.FDI.max_loc(1) spots.cold.FDI.max_loc(1)]) + crs); ylim(mean([spots.hot.FDI.max_loc(2) spots.cold.FDI.max_loc(2)]) + crs);

        zlim([min(zlim) max(zlim) + abs(where_anterior(3))])

        text(min(xlim)+2, min(ylim) + 2, max(zlim)-1, ID, HorizontalAlignment="left", VerticalAlignment="baseline", FontSize=30)


        %if where_anterior(3) < 0
        %    z = max(zlim);
        %else
        %    z = max(zlim) - where_anterior(3);
        %end
        %q = quiver3(max(xlim) - 2, min(ylim) + 2, z, where_anterior(1), where_anterior(2), where_anterior(3));
        %q.LineWidth = 2;
        %q.Color = "k";
        %q.MaxHeadSize=1;
        %q.Marker='.';
        %q.MarkerSize=20;

        fontsize(fig1, scale=1.5)
        exportgraphics(fig1, sprintf('B:/Projects/2023-01 HighDef/Results/AllSpots-Figures/empty-brain_%s_%s.pdf', ID, hemisphere))


        %% Step 2: Mark hotspot
        off = 2;

        ts_hot = trisurf(triangles(triangles_hot,:), coordinates(1,:), coordinates(2,:), coordinates(3,:));
        ts_hot.FaceColor=[1 0 0];
        ts_hot.FaceAlpha=0.5;
        ts_hot.EdgeColor="none";

        pos = spots.hot.FDI.max_loc;
        pos = rotate_to_face_normal * pos';

        if best_R2_hotspot > 0.1
            mark = 'o';
        else
            mark = 'x';
        end

        plot3(pos(1), pos(2), max(zlim) - off, MarkerSize=20, LineWidth=2, Marker=mark, MarkerFaceColor="red", MarkerEdgeColor="k")

        exportgraphics(fig1, sprintf('B:/Projects/2023-01 HighDef/Results/AllSpots-Figures/only-hotspot_%s_%s.pdf', ID, hemisphere))

        %% Step 3: Mark coldspot

        ts_hot.Visible = "off";

        ts_hot_only = trisurf(triangles(triangles_only_hot,:), coordinates(1,:), coordinates(2,:), coordinates(3,:));
        ts_hot_only.FaceColor=[1 0 0];
        ts_hot_only.FaceAlpha=0.5;
        ts_hot_only.EdgeColor="none";

        ts_cold_only = trisurf(triangles(triangles_only_cold,:), coordinates(1,:), coordinates(2,:), coordinates(3,:));
        ts_cold_only.FaceColor=[0 0 1];
        ts_cold_only.FaceAlpha=0.5;
        ts_cold_only.EdgeColor="none";

        ts_both = trisurf(triangles(triangles_both,:), coordinates(1,:), coordinates(2,:), coordinates(3,:));
        ts_both.FaceColor=[0.5 0 0.5];
        ts_both.FaceAlpha=0.5;
        ts_both.EdgeColor="none";

        pos = spots.cold.FDI.max_loc;
        pos = rotate_to_face_normal * pos';

        if best_R2_coldspot > 0.1
            mark = 'o';
        else
            mark = 'x';
        end

        plot3(pos(1), pos(2), max(zlim) - off, MarkerSize=20, LineWidth=2, Marker=mark, MarkerFaceColor="blue", MarkerEdgeColor="k")

        exportgraphics(fig1, sprintf('B:/Projects/2023-01 HighDef/Results/AllSpots-Figures/both_%s_%s.pdf', ID, hemisphere))


    end
end




