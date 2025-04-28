% Output the results
N = 10; D = 3;
nMC = 1e5;


ps = nan(1, nMC);
for iMC = 1:nMC
    locations_in_cond_A = randn(N, D);
    locations_in_cond_B = randn(N, D);
    ps(iMC) = hotelling_test(locations_in_cond_A, locations_in_cond_B);
end

fprintf("p < 0.05 in %0.3f %% of simulations\n", 100 * mean(ps < 0.05))
% This seems to be quite safe!


%% Minimum covariance determinant Monte Carlo:
% See here: https://doi.org/10.1007/s001840200192
nMC = 1e3;
d = nan;
q = nan;
ps = nan(1, nMC);
for iMC = 1:nMC
    locations_in_cond_A = randn(N, D);
    locations_in_cond_B = randn(N, D);
    [p_value, ~, ~, d, q] = MCD_test(locations_in_cond_A, locations_in_cond_B, d, q);
    ps(iMC) = p_value;
    if mod(iMC, 1000) == 0
        fprintf(".")
    end
end
fprintf("\np < 0.05 in %0.3f %% of simulations\n", 100 * mean(ps < 0.05))

% Then we have: TR2 \approx d*F(p,q)
% Turn into:
% Compute the TR2 for the real data
% F = TR2 / d;
% p_value = 1 - fcdf(F, D, q)


%% Thus:
% Obtain the positions of all the spots

%% 1 Collect all spot data available for all subjects as a joint table:
%    Subject Session Hemisphere Spot Muscle x_source y_source z_source x_coil y_coil z_coil angle_coil R2
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



%%
% Compare FDI hot vs FDI cold for left and right hemisphere, each at
% alpha=2.5 %

alpha = 0.025;

tiledlayout(2,2);
i = 1;

rejected_subject = ""; %"sub-011";
% This subject had very unstable MEP amplitudes (also in the experience of
% other experimenters), and has also a fairly bad hotspot localization (R²
% < 0.5)


for dominant = [true, false]
    if dominant
        hemisphere = "dominant   ";
    else
        hemisphere = "nondominant";
    end
    indices_A = find(Result.Spot == "hot" & Result.Muscle == "FDI" & Result.Session == "S2" & Result.Dominant == dominant);
    indices_B = arrayfun(@(i) find(Result.Spot == "cold" & Result.Subject == Result.Subject(i) & Result.Dominant == Result.Dominant(i) & Result.Muscle == Result.Muscle(i) & Result.Session == "S2"), indices_A);

    mask_A_found = Result.R2(indices_A) > 0.1;
    mask_B_found = Result.R2(indices_B) > 0.1;
    mask = mask_A_found & mask_B_found & ~strcmpi(rejected_subject, Result.Subject(indices_A));
    n = sum(mask);

    indices_A = indices_A(mask);
    indices_B = indices_B(mask);

    locations_A = [Result.X_source(indices_A) Result.Y_source(indices_A) Result.Z_source(indices_A)];
    locations_B = [Result.X_source(indices_B) Result.Y_source(indices_B) Result.Z_source(indices_B)];

    [p_hotelling, T2, F] = hotelling_test(locations_A, locations_B);
    [p_mcd] = 1; %MCD_test(locations_A, locations_B);

    differences = locations_B - locations_A;
    % project onto postero-lateral unit vector:
    if dominant
        compare_along = [-1 -1 0] ./ sqrt(2); % flat posterolateral (which, due to curvature of the head means: outwards
        compare_along = [-1 -1 -1] ./ sqrt(3); % downward posterolateral (more along the surface)
    else
        compare_along = [1 -1 -1] ./ sqrt(3);
    end
    scoring = differences * compare_along';
    p_projected = signrank(scoring, 0, "tail", "right");

    %[p, T2, F] = MCD_test(locations_A, locations_B);

    fprintf("%s Hemisphere FDI hot vs cold:\n    Hotelling \t\t\tp = %0.4f (N=%d) \t %s \n         MCD \t\t\tp = %0.4f (N=%d) \t %s \n    Projection+Ranksum \tp = %0.4f (N=%d) \t %s\n", hemisphere, p_hotelling, n, enstar(p_hotelling, alpha), p_mcd, n, enstar(p_mcd, alpha), p_projected, n, enstar(p_projected, alpha))

    ax = nexttile(i);
    plot3(locations_A(:,1), locations_A(:,2), locations_A(:,3), 'r+')
    hold on
    plot3(locations_B(:,1), locations_B(:,2), locations_B(:,3), 'bx')
    plot3([locations_A(:,1) locations_B(:,1)]', [locations_A(:,2) locations_B(:,2)]', [locations_A(:,3) locations_B(:,3)]', '-k')
    l = max([range(xlim) range(ylim) range(zlim)]);
    xlim(mean(xlim) + 0.5.*[-l l]); ylim(mean(ylim) + 0.5.*[-l l]); zlim(mean(zlim) + 0.5.*[-l l])
    axis square; box on; grid on
    ax.Color = "none";
    title(sprintf("%s: absolute", hemisphere))
    view(2)

    differences = locations_B - locations_A;

    ax = nexttile(i+2);
    plot3([zeros(n, 1) differences(:, 1)]', [zeros(n, 1) differences(:, 2)]', [zeros(n, 1) differences(:, 3)]', 'ko-', MarkerFaceColor="k")
    l = max([range(xlim) range(ylim) range(zlim)]);
    xlim(mean(xlim) + 0.5.*[-l l]); ylim(mean(ylim) + 0.5.*[-l l]); zlim(mean(zlim) + 0.5.*[-l l])
    axis square; box on; grid on
    ax.Color = "none";
    title(sprintf("%s: differences", hemisphere))
    view(2)
    i = i+1;
end



%%
% Additionally, do the tests (but uncorrected) for ADM vs APB, APB vs FDI,
% ADM vs FDI
alpha = 0.05 / 3;

comparisons = ["ADM" "APB"; "ADM" "FDI"; "APB" "FDI"];

for dominant = [true, false]
    if dominant
        hemisphere = "dominant   ";
    else
        hemisphere = "nondominant";
    end
    for iComparison = 1:size(comparisons, 1)
        indices_A = find(Result.Spot == "hot" & Result.Muscle == comparisons(iComparison, 1) & Result.Session == "S2" & Result.Dominant == dominant);
        indices_B = arrayfun(@(i) find(Result.Spot == "hot" & Result.Subject == Result.Subject(i) & Result.Dominant == Result.Dominant(i) & Result.Muscle == comparisons(iComparison, 2) & Result.Session == "S2"), indices_A);

        mask_A_found = Result.R2(indices_A);
        mask_B_found = Result.R2(indices_B);
        mask = mask_A_found & mask_B_found;
        n = sum(mask);

        indices_A = indices_A(mask);
        indices_B = indices_B(mask);

        locations_A = [Result.X_source(indices_A) Result.Y_source(indices_A) Result.Z_source(indices_A)];
        locations_B = [Result.X_source(indices_B) Result.Y_source(indices_B) Result.Z_source(indices_B)];

        [p_hotelling, T2, F] = hotelling_test(locations_A, locations_B);

        fprintf("%s Hemisphere: p = %0.4f (N=%d)   %s (%s hot vs %s hot)\n", hemisphere, p_hotelling, n, enstar(p_hotelling, alpha), comparisons(iComparison, 1), comparisons(iComparison, 2))
    end
end


function [p_value, T2, F] = hotelling_test(A, B)
% A: N x P matrix
% B: N x P matrix
% Performs a Hotelling T² test to check whether the paired difference
% between A and B is nonzero.
% A(i,:) is treated as the paired comparant of B(i,:), e.g. A(i,:) is from
% condition A, and B(i,:) from condition B, and both from subject i
%
% See here: https://doi.org/10.1007/s001840200192
D = A - B;

mean_diff = mean(D);
S = cov(D);

N = size(D, 1);
dof = size(D, 2);
T2 = N * (mean_diff / S) * mean_diff';
F = ((N-dof) * T2) / (dof * (N-1));

p_value = 1 - fcdf(F, dof, N-dof);
end


function [p_value, TR2, F, d, q] = MCD_test(A, B, d, q, nMC)
% Performs MC approximation of d and q iff neither d nor q is given
% See here: https://doi.org/10.1007/s001840200192
arguments
    A {mustBeNumeric};
    B {mustBeNumeric};
    d (1,1) {mustBeNumeric} = nan;
    q (1,1) {mustBeNumeric} = nan;
    nMC (1,1) {mustBeNumeric, mustBeInteger, mustBePositive} = 1e4;
end

D = A - B;
N = size(D, 1);
dof = size(D, 2);

if isnan(d) || isnan(q)
    TR2 = nan(1, nMC);
    for iMC = 1:nMC
        X = randn(N,dof); % For given sample size and dimensionality N x D
        [C, T] = robustcov(X, Method="fmcd"); % OutlierFraction = 0.95 may also be worthwhile
        TR2(iMC) = N * (T / C) * T'; % This entails that the tested mean vector mu is zero!
    end

    % Then estimate the expected value of TR2 by mean(TR2), and its variance by
    % var(TR2):
    % This yields: https://doi.org/10.1007/s001840200192 Equations 8 and 7
    q = 4 + (dof + 2) / ((dof * var(TR2))/(mean(TR2)^2 * 2) - 1);
    d = mean(TR2) * (q - 2)/q;
end

[C, T] = robustcov(D, 'Method', 'fmcd');
TR2 = N * (T / C) * T';
F = TR2 / d;
p_value = 1 - fcdf(F, dof, q);
end



function s = enstar(p, alpha)
arguments
    p (1,1) {mustBeNumeric, mustBeNonnegative};
    alpha (1,1) {mustBeNumeric, mustBeNonnegative} = 0.05;
end
if p < alpha
    s = "*";
else
    s = "-";
end
end