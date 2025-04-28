true_x = 3;
true_y = 5;
spread = 2;
noise_level = 0.5;

x_lim = [-10, 10];
y_lim = [-10, 10];

plot(true_x, true_y, 'r+', MarkerSize=30, LineWidth=3)
hold on
xlim(x_lim)
ylim(y_lim)

n_samples = 1000;
n_mc_iterations = 2000;
x_est = nan(n_mc_iterations, 1);
y_est = nan(n_mc_iterations, 1);
for iMC = 1:n_mc_iterations
    sampled_x = range(x_lim) .* rand(n_samples, 1) + x_lim(1);
    sampled_y = range(y_lim) .* rand(n_samples, 1) + y_lim(1);
    noiseless = exp(-((sampled_x - true_x).^2 / (2*spread^2)) -((sampled_y - true_y).^2 / (2*spread^2)));
    noise = noise_level .* randn(n_samples, 1);

    response = noiseless + noise;
    cog_x = (response' * sampled_x) ./ sum(response);
    cog_y = (response' * sampled_y) ./ sum(response);
    x_est(iMC) = cog_x;
    y_est(iMC) = cog_y;
    plot(cog_x, cog_y, 'k.')
    drawnow
end

plot(true_x, true_y, 'r+', MarkerSize=30, LineWidth=3)
plot(mean(x_est), mean(y_est), 'yx', LineWidth=2)
plot(median(x_est), median(y_est), 'yo', LineWidth=2)
% mean is is slightly biased (?), in median though, COG is good