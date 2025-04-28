ID = 'sub-001';
exp_ID = 'R2L';
record = readtable(sprintf('B:/Projects/2023-01 HighDef/Results/%s_%s_raw.csv', ID, exp_ID));

nTrials = size(record, 1);

transforms = zeros(4, 4, nTrials);
transforms(4,4,:) = 1;
dimnames = {'x', 'y', 'z', 'p'};
for j = 1:4
    for i = 1:3
        transforms(i,j,:) = record.(sprintf('%s%d', dimnames{j}, i));
    end
end

positions = squeeze(transforms(1:3,4,:))';
%%
transforms = se3(transforms);

angles = eul(transforms);
roll_angles = angles(:,3);

complex_representation = exp(1i .* roll_angles);
mean_roll = angle(mean(complex_representation));

%%
figure;
polarhistogram(roll_angles, 60)
hold on
polarplot([mean_roll mean_roll], [0 100], 'r')

%%
% Project to 2D
[projectors, scores] = pca(positions);

figure;
subplot(1,2,1)
scatter3(positions(:,1),positions(:,2),positions(:,3), 'k.')
hold on
s = 10;
plot3([0 s*projectors(1,1)], [0 s*projectors(2,1)], [0 s*projectors(3,1)], 'r')
plot3([0 s*projectors(1,2)], [0 s*projectors(2,2)], [0 s*projectors(3,2)], 'g')
plot3([0 s*projectors(1,3)], [0 s*projectors(2,3)], [0 s*projectors(3,3)], 'b')
limit_range = max([range(xlim) range(ylim) range(zlim)]);
window = 0.5.*[-limit_range limit_range];
xlim(mean(xlim) + window); ylim(mean(ylim) + window); zlim(mean(zlim) + window)

axis vis3d
title('Raw locations and PCs')

subplot(1,2,2)
scatter(scores(:,1), scores(:,2), 'k.')
title('Projected to first two PCs')



%%
% Assign to each point two values:
% 1) Number of neighbors (distance-weighted)
% 2) Variability of orientation of neighbors (distance-weighted)

x = scores(:,1);
y = scores(:,2);

[X,Y,position_coverage, angle_coverage, mean_angles] = diagnose_coverage(x,y, roll_angles, 3);

figure;
subplot(2,2,1)
contourf(X,Y,position_coverage, EdgeColor='none')
colorbar
hold on
plot(x,y,'w.')
title('Position Coverage (~Density of samples)')

subplot(2,2,2)
contourf(X,Y,angle_coverage, 20, EdgeColor='none')
colorbar
hold on
plot(x,y,'w.')
title('Angle coverage (Variability of angles sampled nearby): Uniform gives 1.1-1.15')

subplot(2,2,4)
contourf(X,Y,rad2deg(mean_angles), 60, EdgeColor='none')
colorbar
hold on
plot(x,y,'w.')
title('Mean Angle sampled nearby')


%% Uniform example for comparison:
x = (rand(1e5,1) .* range(x)) + min(x);
y = (rand(1e5,1) .* range(y)) + min(y);
a = (rand(1e5,1) .* deg2rad(90)) - deg2rad(45);

[X,Y,position_coverage, angle_coverage, mean_angles] = diagnose_coverage(x,y, a, 3);

figure;
subplot(2,2,1)
contourf(X,Y,position_coverage, EdgeColor='none')
colorbar
title('Position Coverage (~Density of samples)')

subplot(2,2,2)
contourf(X,Y,angle_coverage, EdgeColor='none')
colorbar
title('Angle coverage (Variability of angles sampled nearby)')

subplot(2,2,4)
contourf(X,Y,rad2deg(mean_angles), 60, EdgeColor='none')
colorbar
title('Mean Angle sampled nearby')




function [X,Y,position_coverage, angle_coverage, mean_angles] = diagnose_coverage(x,y,a, smoothness)
xlower = floor(min(x)/10)*10;
ylower = floor(min(y)/10)*10;
xupper = ceil(max(x)/10)*10;
yupper = ceil(max(y)/10)*10;

nSamples = 100;
[X,Y] = meshgrid(linspace(xlower, xupper, nSamples), linspace(ylower, yupper, nSamples));

position_coverage = nan(nSamples, nSamples);
angle_coverage = nan(nSamples, nSamples);
mean_angles = nan(nSamples, nSamples);
for i = 1:nSamples
    for j = 1:nSamples
        distance_weights = exp(-0.5.*((x - X(i,j)).^2 + (y-Y(i,j)).^2)./(smoothness^2))./(2*pi*smoothness*smoothness);
        position_coverage(i,j) = sum(distance_weights);
        angle_coverage(i,j) = 1/abs(sum(distance_weights .* exp(1i .* a)) / sum(distance_weights));
        mean_angles(i,j) = angle(sum(distance_weights .* exp(1i .* a)) / sum(distance_weights));
    end
end
end