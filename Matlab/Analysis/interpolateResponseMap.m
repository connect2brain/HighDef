function [Z, evalX, evalY] = interpolateResponseMap(locations, response, gridLength, sd, kernel, threshold, doMedian)
arguments
    locations (:,2) {mustBeNumeric}
    response (1,:) {mustBeNumeric}
    gridLength (1,1) {mustBeNumeric,mustBePositive,mustBeInteger} = 30
    sd (1,1) {mustBeNumeric, mustBePositive} = 2;
    kernel function_handle = @(x,y) (1/(2*pi*sd^2)).*exp(-(x.^2 + y.^2)./(2*sd^2));
    threshold (1,1) {mustBeNumeric} = -Inf;
    doMedian = false;
end

axmin = min(locations(:));
axmax = max(locations(:));
[evalX, evalY] = meshgrid(linspace(axmin, axmax, gridLength), linspace(axmin, axmax, gridLength));

Z = nan(gridLength);
for i = 1:gridLength
    for j = 1:gridLength
        dx = locations(:,1) - evalX(i,j);
        dy = locations(:,2) - evalY(i,j);
        weights = kernel(dx, dy);
        normalizer = sum(weights);
        % thresholding:
        if normalizer > threshold
            if doMedian
                Z(i,j) = weightedMedian(response, weights./sum(weights));
            else
                Z(i,j) = (response * weights)/sum(weights);
            end            
        else
            Z(i,j) = nan;
        end
    end
end
end