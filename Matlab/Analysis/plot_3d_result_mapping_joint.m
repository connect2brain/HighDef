DO_STREAMLINES = true;



for ID = ["sub-001", "sub-002", "sub-003", "sub-004", "sub-006", "sub-007", "sub-008", "sub-009", "sub-011", "sub-012","sub-013", "sub-014"]
    for exp_id = ["R2L", "L2R"]
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
        % lower weights to APB and ADM hotspot
        %avg_normal = [1 1 0.5 0.5] * [spots.hot.FDI.normal; spots.cold.FDI.normal; spots.hot.APB.normal; spots.hot.ADM.normal];
        %avg_normal = avg_normal + [0 0 2]; % Turn slightly upwards
        %avg_normal(2) = 0; % Force into x-z-plane

        %avg_normal = avg_normal ./ norm(avg_normal);

        % Use vector from coldspot to hotspot instead:
        connection_hot_cold = spots.hot.FDI.max_loc - spots.cold.FDI.max_loc;
        connection_hot_cold(2) = 0; % Force into x-z-plane

        if norm(connection_hot_cold) > 0
            connection_hot_cold = connection_hot_cold ./ norm(connection_hot_cold);
        else
            connection_hot_cold = [-1 0 0];
        end

        % Take average with simple upwards vector
        connection_hot_cold = 0.5.* connection_hot_cold + 0.5.*[-1 0 0];

        % FROM: https://de.mathworks.com/matlabcentral/answers/445994-how-to-calculate-a-rotation-matrix-between-two-3d-points#answer_361888
        % two random 3D vectors

        % (a) Normals-based
        %p0 = avg_normal;
        %p1 = [0 0 1];
        % (b) Connection-based
        p0 = connection_hot_cold;
        p1 = [-1 0 0];

        fprintf('Transforming coordinates from  %d, %d, %d  to  %d, %d, %d\n', p0(1), p0(2), p0(3), p1(1), p1(2), p1(3))

        % calculate cross and dot products
        C = cross(p0, p1);
        D = dot(p0, p1);
        NP0 = norm(p0); % used for scaling
        if ~all(C==0) % check for colinearity
            Z = [0 -C(3) C(2); C(3) 0 -C(1); -C(2) C(1) 0] ;
            rotate_to_face_normal = (eye(3) + Z + Z^2 * (1-D)/(norm(C)^2)) / NP0^2 ; % rotation matrix
        else
            rotate_to_face_normal = sign(D) * (norm(p1) / NP0) ; % orientation and scaling
        end

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

        nTriangles = size(triangles,1);
        centers = nan(3, nTriangles);
        normals = nan(3, nTriangles);
        for iTriangle = 1:nTriangles
            node_coordinates = coordinates(:, triangles(iTriangle,:));
            centers(:,iTriangle) = mean(node_coordinates, 2);
            normals(:,iTriangle) = cross(node_coordinates(:,2) - node_coordinates(:,1), node_coordinates(:,3) - node_coordinates(:,1));
            normals(:,iTriangle) = normals(:,iTriangle) ./ norm(normals(:,iTriangle)); % unit length
        end


        % Pick top 5 % of triangles from hotspot FDI, plot in red; pick top 5 % of
        % triangles from coldspot FDI, color blue; color intersection purple.
        hotspot_data = h5read(sprintf('//wsl.localhost/Ubuntu-22.04/home/bnplab-admin/TMS_localization/HighDef/%s/results/exp_map-%s/r2/mesh_%s/roi_%s/CsE_FDI_in_uV/sigmoid4/r2_roi_data.hdf5', ID, exp_id, mesh_id, roi), '/data/tris/c_E_mag');
        coldspot_data = h5read(sprintf('//wsl.localhost/Ubuntu-22.04/home/bnplab-admin/TMS_localization/HighDef/%s/results/exp_map-%s/r2/mesh_%s/roi_%s/SIHIscore_FDI/sigmoid4/r2_roi_data.hdf5', ID, exp_id, mesh_id, roi), '/data/tris/c_E_mag');

        hotspot_data(isnan(hotspot_data)) = 0;
        hotspot_data = reject_R2_outliers_if_needed(hotspot_data);
        coldspot_data(isnan(coldspot_data)) = 0;
        coldspot_data = reject_R2_outliers_if_needed(coldspot_data);

        top_percent = 2.5;

        %triangles_hot = hotspot_data > quantile(hotspot_data, 1 - (top_percent / 100));
        %triangles_cold = coldspot_data > quantile(coldspot_data, 1 - (top_percent / 100));
        triangles_hot  = selectTopQuantile(hotspot_data,  top_percent);
        triangles_cold = selectTopQuantile(coldspot_data, top_percent);
        triangles_hot  = selectRelativeToBestR2(hotspot_data,  0.5);
        triangles_cold = selectRelativeToBestR2(coldspot_data, 0.5);
        

        triangles_both = triangles_hot & triangles_cold;
        triangles_only_hot = triangles_hot & ~triangles_cold;
        triangles_only_cold = triangles_cold & ~triangles_hot;


        fig = figure(Position=[50 50 600 600]);
        hold on
        ts = trisurf(triangles, coordinates(1,:), coordinates(2,:), coordinates(3,:), mean(raw_excentricity(triangles), 2));
        ts.EdgeColor = 'none';
        ts.FaceAlpha = 1;
        %ts.FaceColor = "w";
        cr = max([range(xlim) range(ylim) range(zlim)]);
        xlim(mean(xlim) + [-cr cr]./2); ylim(mean(ylim) + [-cr cr]./2); zlim(mean(zlim) + [-cr cr]./2)
        clim([0.6*min(clim) 1.05*max(clim)])
        set(gca(), 'Color', 'none')
        colormap("gray")
        grid on; grid minor; box on; axis square; view(2)
        xticks(-100:5:100); yticks(-100:5:100);
        labels = arrayfun(@(v) sprintf("%d", v), xticks); labels(mod(xticks, 10) ~= 0) = "";
        xticklabels(labels); yticklabels(labels)
        %title(sprintf("%s / %s", ID, hemisphere))

        ts_hot = trisurf(triangles(triangles_only_hot,:), coordinates(1,:), coordinates(2,:), coordinates(3,:));
        ts_hot.FaceColor=[1 0 0];
        ts_hot.FaceAlpha=0.5;
        ts_hot.EdgeColor="none";

        ts_both = trisurf(triangles(triangles_both,:), coordinates(1,:), coordinates(2,:), coordinates(3,:));
        ts_both.FaceColor=[0.5 0 0.5];
        ts_both.FaceAlpha=0.5;
        ts_both.EdgeColor="none";

        ts_cold = trisurf(triangles(triangles_only_cold,:), coordinates(1,:), coordinates(2,:), coordinates(3,:));
        ts_cold.FaceColor=[0 0 1];
        ts_cold.FaceAlpha=0.5;
        ts_cold.EdgeColor="none";




        xs = []; ys = [];
        for spot = ["hot", "cold"]
            for muscle = ["APB", "ADM", "FDI"]
                if spot == "cold" && muscle ~= "FDI"
                    continue
                end
                pos = spots.(spot).(muscle).max_loc;

                pos = rotate_to_face_normal * pos';

                n = spots.(spot).(muscle).normal;
                l = 0;
                plot3(pos(1) + l.*[0 n(1)], pos(2) + l.*[0 n(2)], pos(3) + l.*[0 n(3)], Color=colors.(spot).(muscle), Marker=".")


                if DO_STREAMLINES
                    nSteps = 60;
                    %constantforce = [-1;0;1];
                    constantforce = pos;
                    constantforce = 0.67 .* constantforce ./ norm(constantforce);
                    p = nan(3,nSteps);
                    p(:,1) = pos;
                    for iStep = 2:nSteps
                        weights = exp(-vecnorm(centers - p(:, iStep-1)));
                        force = sum(weights .* normals, 2) / sum(weights);
                        p(:, iStep) = p(:, iStep-1) + force + constantforce;
                        linelength = vecnorm(diff(p,[],2));
                        linelength = sum(linelength(~isnan(linelength)));
                        if linelength > 20
                            break
                        end
                    end
                    plot3(p(1,:), p(2,:), p(3,:), colors.(spot).(muscle), LineWidth=2)
                end



                if muscle == "FDI"
                    s = 20;
                    off = 2;
                    mfc = colors.(spot).FDI;
                    mec = "w";
                    lw = 1;
                else
                    s = 14;
                    off = 1;
                    mfc = "none";
                    mec = colors.(spot).(muscle);
                    lw = 2;
                end


                if DO_STREAMLINES
                    plot3(p(1,iStep), p(2,iStep), max(zlim) - off, MarkerSize=s, LineWidth=lw, Marker=markers.(spot).(muscle), MarkerFaceColor=mfc, MarkerEdgeColor=mec, Color="none")
                    xs = [xs pos(1) p(1,iStep)]; ys = [ys pos(2) p(2, iStep)];
                else
                    plot3(pos(1), pos(2), max(zlim) - off, MarkerSize=s, LineWidth=lw, Marker=markers.(spot).(muscle), MarkerFaceColor=mfc, MarkerEdgeColor=mec, Color="none")
                    xs = [xs pos(1)]; ys = [ys pos(2)];
                end

                %if n(3) < 0
                %    plot3(pos(1), pos(2), pos(3)+5, MarkerSize=s, LineWidth=2, Marker=markers.(spot).(muscle), MarkerFaceColor=colors.(spot).(muscle), MarkerEdgeColor=colors.(spot).(muscle))
                %    plot3(pos(1), pos(2), max(zlim) - 1, MarkerSize=s, LineWidth=2, Marker=markers.(spot).(muscle), MarkerFaceColor="none", MarkerEdgeColor=colors.(spot).(muscle))
                %else
                %    plot3(pos(1), pos(2), pos(3)+5, MarkerSize=s, LineWidth=2, Marker=markers.(spot).(muscle), MarkerFaceColor=colors.(spot).(muscle), MarkerEdgeColor=colors.(spot).(muscle))
                %    %plot3(pos(1), pos(2), max(zlim) - 1, MarkerSize=s, LineWidth=2, Marker=markers.(spot).(muscle), MarkerFaceColor="none", MarkerEdgeColor=colors.(spot).(muscle))
                %end
                drawnow;
                pause(0.1);
            end
        end

        % S1-hotspot:
        pos = S1hotspots.FDI.max_loc;
        pos = rotate_to_face_normal * pos';

        xs = [xs pos(1)]; ys = [ys pos(2)];
        n = S1hotspots.FDI.normal;
        %plot3(pos(1), pos(2), max(zlim) - 2, MarkerSize=10, LineWidth=2, Marker="o", MarkerFaceColor="none", MarkerEdgeColor=colors.hot.(muscle))




        where_anterior = rotate_to_face_normal * [0 10 0]';

        %cr = max([range(xlim) range(ylim) range(zlim)]);
        shrink = 0.4;
        crs = shrink*0.5.*[-cr cr];
        xlim(mean(xs) + crs); ylim(mean(ys) + crs);

        zlim([min(zlim) max(zlim) + abs(where_anterior(3))])

        if strcmpi(exp_id, "R2L")
            text(max(xlim)-2, min(ylim) + 2, max(zlim)-1, ID, HorizontalAlignment="right", VerticalAlignment="baseline", FontSize=30)
        else
            text(min(xlim)+2, min(ylim) + 2, max(zlim)-1, ID, HorizontalAlignment="left", VerticalAlignment="baseline", FontSize=30)
        end

        if DO_STREAMLINES
            legend(["", "", "", "", "", "", "", "", "", "", "", "", sprintf("R² = %0.2f", max(hotspot_data)), "", "", sprintf("R² = %0.2f", max(coldspot_data))], FontSize=17)
        else
            legend(["", "", "", "", "", "", "", "", "", sprintf("R² = %0.2f", max(hotspot_data)), "", sprintf("R² = %0.2f", max(coldspot_data))], FontSize=17)
        end


        % if where_anterior(3) < 0
        %     z = max(zlim);
        % else
        %     z = max(zlim) - where_anterior(3);
        % end
        % q = quiver3(max(xlim) - 2, min(ylim) + 2, z, where_anterior(1), where_anterior(2), where_anterior(3));
        % q.LineWidth = 2;
        % q.Color = "k";
        % q.MaxHeadSize=1;
        % q.Marker='.';
        % q.MarkerSize=20;

        % where_lateral = rotate_to_face_normal * [-10 0 0]';
        % q = quiver3(max(xlim) - 2, min(ylim) + 2, max(zlim) - where_lateral(3), where_lateral(1), where_lateral(2), where_lateral(3));
        % q.LineWidth = 2;
        % q.Color = "#00ffff";
        % q.MaxHeadSize=1;

        fontsize(fig, scale=1.5)
        if DO_STREAMLINES
            exportgraphics(fig, sprintf('B:/Projects/2023-01 HighDef/Results/AllSpots-Figures/all_spots_streamlines_%s_%s.pdf', ID, hemisphere))
        else
            exportgraphics(fig, sprintf('B:/Projects/2023-01 HighDef/Results/AllSpots-Figures/all_spots_%s_%s.pdf', ID, hemisphere))
        end

        close all
    end
end


function mask = selectTopQuantile(data, q)
% q in % (i.e. q=2.5 means 2.5%)
mask = data > quantile(data, 1 - (q / 100));
end


function mask = selectRelativeToBestR2(data, q)
mask = data > q*max(data);
end