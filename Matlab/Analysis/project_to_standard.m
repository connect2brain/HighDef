function projectedLocations = project_to_standard(table, hemisphere)
% table must have columns p1, p2, p3 encoding the 3d position of the coil
% Result is a n_samples x 2 matrix, where the first column is the location
% in posterior-anterior direction, and second column is locations in
% lateral-medial direction
% Hemisphere: either "L" or "R"

% For projection, ignore rotation part of matrix
locations = [table.p1 table.p2 table.p3];
principalComponents = pca(locations);
flatteningProjection = principalComponents(:,1:2);

projectedLocations = locations - mean(locations,1);
projectedLocations = projectedLocations * flatteningProjection;

towardsAnterior = [0 1 0] * flatteningProjection;
towardsAnterior = towardsAnterior ./ sqrt(towardsAnterior * towardsAnterior');
if strcmpi(hemisphere, "L")
    towardsMedial = [sqrt(0.5) 0 sqrt(0.5)] * flatteningProjection;
else
    towardsMedial = [sqrt(0.5) 0 -sqrt(0.5)] * flatteningProjection;
end
towardsMedial = towardsMedial ./ sqrt(towardsMedial * towardsMedial');

changeToPhysiologicalBasis = [towardsAnterior; towardsMedial];
projectedLocations = (changeToPhysiologicalBasis * projectedLocations')';

end