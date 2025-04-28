sessionNames = []; sessionNames.L = ["map-L", "map-L2R"]; sessionNames.R = ["map-R", "map-R2L"];
templates = []; templates.hot = "CsE_%s_in_uV"; templates.cold = "SIHIscore_%s";
basepath = '//wsl.localhost/Ubuntu-22.04/home/bnplab-admin/TMS_localization/HighDef';
subjects = arrayfun(@(s) string(s.name), dir(fullfile(basepath, 'sub-*')))';

Result = [];
for ID = subjects
    for hemisphere = ["R", "L"]
        for muscle = ["ADM", "FDI", "APB"]
            roi_id = sprintf("midlayer_%s", lower(hemisphere));
            for session = 1:2
                exp_id = sessionNames.(hemisphere)(session);
                if session == 2 && muscle == "FDI"
                    spots = ["hot", "cold"];
                else
                    spots = ["hot"];
                end

                for spot = spots
                    response_id = sprintf(templates.(spot), muscle);
                    newrow = collectSpot(basepath, ID, exp_id, roi_id, response_id);
                    if ~isempty(newrow)
                        newrow.Session = sprintf("S%d", session); newrow.Hemisphere = hemisphere; newrow.Spot = spot; newrow.Muscle = muscle;
                        fprintf(' %s/%s  :  S%d=%s %s-%s\n', ID, hemisphere, session, exp_id, muscle, spot)
                    end
                    Result = [Result; newrow];
                end
            end
        end
    end
end



% Add new hemisphere labels: dominant, nondominant
Result.Dominant = strcmpi(Result.Hemisphere, "L");

% For the single right-handed subject, flip this
% and also mirror the x-coordinate
mirror_which = strcmpi(Result.Subject, "sub-006");
Result.Dominant(mirror_which) = ~Result.Dominant(mirror_which);
Result.X_source(mirror_which) = -Result.X_source(mirror_which);




%% Hot-cold

fig_hotcold = figure(Position=[50 50 600 350]);
tiledlayout(1,2, TileSpacing="compact");

fig_round_dominant    = figure(Position=[100 50 300 260], Name="Dominant Round");
fig_round_nondominant = figure(Position=[100 100 300 260], Name="Nondominant Round");


for dominant = [true, false]
    if dominant
        hemisphere = "dominant";
        fig_round = fig_round_dominant;
    else
        hemisphere = "nondominant";
        fig_round = fig_round_nondominant;
    end
    indices_A = find(Result.Spot == "hot" & Result.Muscle == "FDI" & Result.Session == "S2" & Result.Dominant == dominant);
    indices_B = arrayfun(@(i) find(Result.Spot == "cold" & Result.Subject == Result.Subject(i) & Result.Dominant == Result.Dominant(i) & Result.Muscle == Result.Muscle(i) & Result.Session == "S2"), indices_A);

    mask_A_found = Result.R2(indices_A) > 0.1;
    mask_B_found = Result.R2(indices_B) > 0.1;
    mask = mask_A_found & mask_B_found;
    n = sum(mask);

    indices_A = indices_A(mask);
    indices_B = indices_B(mask);

    locations_A = [Result.X_source(indices_A) Result.Y_source(indices_A) Result.Z_source(indices_A)];
    locations_B = [Result.X_source(indices_B) Result.Y_source(indices_B) Result.Z_source(indices_B)];

    differences = locations_B - locations_A;


    % Joint rectangular figure
    set(0, 'CurrentFigure', fig_hotcold)
    ax = nexttile;
    mx = mean(differences(:,1)); my = mean(differences(:,2)); mz = mean(differences(:,3));
    plot3(mx, my, mz, 'xk', LineWidth=3, MarkerSize=10)
    text(mx,my,mz, "mean", HorizontalAlignment="left", VerticalAlignment="top")
    hold on
    plot3([zeros(n, 1) differences(:, 1)]', [zeros(n, 1) differences(:, 2)]', [zeros(n, 1) differences(:, 3)]', 'k-')
    plot3(differences(:, 1), differences(:, 2), differences(:, 3), 'ko', MarkerFaceColor="k", MarkerSize=5)
    
    %l = max([range(xlim) range(ylim) range(zlim)]);
    %xlim(mean(xlim) + 0.5.*[-l l]); ylim(mean(ylim) + 0.5.*[-l l]); zlim(mean(zlim) + 0.5.*[-l l])
    xlim([-15 15]); ylim([-20 10])
    xticks(-15:5:15); yticks(-20:5:10)
    axis square; box on; grid on
    ax.Color = "none";
    title([sprintf("%s hemisphere", hemisphere); "coldspot relative to hotspot"])
    xlabel("\leftarrow Left [mm]  /  Right [mm] \rightarrow")
    ylabel("\leftarrow Back [mm]  /  Front [mm] \rightarrow")
    view(2)
    if ~dominant
        ax.YAxisLocation = "right";
    end

    % Separate round figures
    set(0, 'CurrentFigure', fig_round)
    angles = atan2(differences(:,2), differences(:,1));
    distances = sqrt(differences(:,1).^2 + differences(:,2).^2);
    polarplot([0 atan2(my, mx)], [0 sqrt(mx^2 + my^2)], '-', LineWidth=3, Color=[0.85 0 0.95])
    hold on
    polarplot(angles, distances, 'ko', MarkerFaceColor="k", MarkerSize=5);
    
    polarplot([zeros(size(angles)) angles]', [zeros(size(distances)) distances]', 'k-');
    
    

    rticks([5 10]); rticklabels(["5 mm", "10 mm"])

    thetaticks([0 90 180 270]); 
    if dominant
        thetaticklabels(["medial", "anterior", "lateral", "posterior"])
    else
        thetaticklabels(["lateral", "anterior", "medial", "posterior"])
    end
    thetaticklabels([])
    set(gca, Color="none", GridColor=[0.05 0.05 0.05], GridAlpha=0.4)
    

end

exportgraphics(fig_hotcold, "B:/Projects/2023-01 HighDef/Results/Evaluation/hot_cold_relative.pdf")
exportgraphics(fig_round_dominant, "B:/Projects/2023-01 HighDef/Results/Evaluation/hot_cold_relative_round_dominant.pdf")
exportgraphics(fig_round_nondominant, "B:/Projects/2023-01 HighDef/Results/Evaluation/hot_cold_relative_round_nondominant.pdf")
































































%% Hotspots relative to each other: attempt 1 -- only xy coords.
figure(Name="Hotspots relative 1");
tiledlayout(1,2);

for dominant = [true, false]
    indices_A = find(Result.Spot == "hot" & Result.Muscle == "APB" & Result.Session == "S2" & Result.Dominant == dominant);
    indices_B = arrayfun(@(i) find(Result.Spot == "hot" & Result.Subject == Result.Subject(i) & Result.Dominant == Result.Dominant(i) & Result.Muscle == "FDI" & Result.Session == "S2"), indices_A);
    indices_C = arrayfun(@(i) find(Result.Spot == "hot" & Result.Subject == Result.Subject(i) & Result.Dominant == Result.Dominant(i) & Result.Muscle == "ADM" & Result.Session == "S2"), indices_A);

    mask_A_found = Result.R2(indices_A) > 0.1;
    mask_B_found = Result.R2(indices_B) > 0.1;
    mask_C_found = Result.R2(indices_C) > 0.1;
    mask = mask_A_found & mask_B_found & mask_C_found;
    n = sum(mask);

    indices_A = indices_A(mask);
    indices_B = indices_B(mask);
    indices_C = indices_C(mask);

    % For plotting, ignore z coordinate
    locations_A = [Result.X_source(indices_A) Result.Y_source(indices_A)];
    locations_B = [Result.X_source(indices_B) Result.Y_source(indices_B)];
    locations_C = [Result.X_source(indices_C) Result.Y_source(indices_C)];

    locations_A = locations_A - locations_B;
    locations_B = locations_B - locations_B;
    locations_C = locations_C - locations_B;

    angles = -atan2(locations_A(:,2), locations_A(:,1));
    rotation_matrices = [reshape(cos(angles), 1, 1, n), reshape(-sin(angles), 1, 1, n); reshape(sin(angles), 1, 1, n), reshape(cos(angles), 1, 1, n)];

    locations_A = squeeze(pagemtimes(rotation_matrices, reshape(locations_A', 2, 1, n)));
    locations_B = squeeze(pagemtimes(rotation_matrices, reshape(locations_B', 2, 1, n)));
    locations_C = squeeze(pagemtimes(rotation_matrices, reshape(locations_C', 2, 1, n)));

    ax = nexttile();
    z_B = Result.Z_source(indices_B)';
    z_A = Result.Z_source(indices_A)' - z_B;
    z_C = Result.Z_source(indices_C)' - z_B;
    z_B = z_B - z_B;
    plot3([locations_A(1,:); locations_B(1,:); locations_C(1,:)], [locations_A(2,:); locations_B(2,:); locations_C(2,:)], [z_A; z_B; z_C], '-')
    hold on
    plot3(locations_A(1,:), locations_A(2,:), z_A, 'kx', MarkerFaceColor="k")
    plot3(locations_B(1,:), locations_B(2,:), z_B, 'k.')
    plot3(locations_C(1,:), locations_C(2,:), z_C, 'ko', MarkerFaceColor="k")
    axis square; grid on;
    xlim([-50 50]); ylim([-50 50]);
    ax.Color = "none";
    view(2)
end



%% Hotspots relative to each other: attempt 2 -- change of basis
figure(Name="Hotspots relative 2");
tiledlayout(1,2);

for dominant = [true, false]
    indices_A = find(Result.Spot == "hot" & Result.Muscle == "APB" & Result.Session == "S2" & Result.Dominant == dominant);
    indices_B = arrayfun(@(i) find(Result.Spot == "hot" & Result.Subject == Result.Subject(i) & Result.Dominant == Result.Dominant(i) & Result.Muscle == "FDI" & Result.Session == "S2"), indices_A);
    indices_C = arrayfun(@(i) find(Result.Spot == "hot" & Result.Subject == Result.Subject(i) & Result.Dominant == Result.Dominant(i) & Result.Muscle == "ADM" & Result.Session == "S2"), indices_A);

    mask_A_found = Result.R2(indices_A) > 0.1;
    mask_B_found = Result.R2(indices_B) > 0.1;
    mask_C_found = Result.R2(indices_C) > 0.1;
    mask = mask_A_found & mask_B_found & mask_C_found;
    n = sum(mask);

    indices_A = indices_A(mask);
    indices_B = indices_B(mask);
    indices_C = indices_C(mask);

    % For plotting, ignore z coordinate
    locations_A = [Result.X_source(indices_A) Result.Y_source(indices_A) Result.Z_source(indices_A)]';
    locations_B = [Result.X_source(indices_B) Result.Y_source(indices_B) Result.Z_source(indices_B)]';
    locations_C = [Result.X_source(indices_C) Result.Y_source(indices_C) Result.Z_source(indices_C)]';

    connect_A = locations_A - locations_B;
    connect_C = locations_C - locations_B;

    A = nan(3,n);
    C = nan(3,n);
    for i = 1:n
        u_A = connect_A(:,i) ./ norm(connect_A(:,i));
        if isnan(u_A)
            u_A = [1;0;0];
            fprintf("APB and FDI hotspots coincide (i: %d)\n", i)
        end
        if dominant
            u_A = -u_A;
        end
        u_C = connect_C(:,i) ./ norm(connect_C(:,i));
        orth_C = u_C - (u_A' * u_C) .* u_A;
        u_orth_C = orth_C ./ norm(orth_C);
        %fprintf("Product of u_orth_C and u_A: %0.5f\n", u_orth_C' * u_A)
        u_normal = cross(u_A, u_orth_C);
        %fprintf("Norm of u_normal: %0.5f\n", norm(u_normal))
        T = [u_A u_orth_C u_normal]; % since T is orthonormal, T' = inv(T)
        A(:,i) = T \ connect_A(:,i);
        C(:,i) = T \ connect_C(:,i);
    end

    ax = nexttile();
    plot3([zeros(1,n); C(1,:)], [zeros(1,n); C(2,:)], [zeros(1,n); C(3,:)], 'k-')
    hold on
    plot3(A(1,:), A(2,:), A(3,:), 'kx', MarkerFaceColor="k")
    plot3(C(1,:), C(2,:), C(3,:), 'ko', MarkerFaceColor="k")
    plot3(mean(C(1,:)), mean(C(2,:)), mean(C(3,:)), 'go', LineWidth=3)
    axis square; grid on;
    xlim([-20 20]); ylim([-20 20]); zlim([-20 20])
    ax.Color = "none";
    view(2)
end



