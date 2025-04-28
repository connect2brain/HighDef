function update_spots(ID, experiment, hemisphere, muscle, spot, x, y, z, quality)
%UPDATE_SPOTS Summary of this function goes here
%   Detailed explanation goes here
write_to = 'B:/Projects/2023-01 HighDef/Results/Evaluation/spot_locations.csv';
T = readtable(write_to);
matchingRow = strcmpi(T.Subject, ID) & strcmpi(T.Session, experiment) & strcmpi(T.Muscle, muscle) & strcmpi(T.Spot, spot) & strcmpi(T.Hemisphere, hemisphere);
newRow = table({ID}, {experiment}, {hemisphere}, {spot}, {muscle}, x, y, z, quality, VariableNames=T.Properties.VariableNames);
if any(matchingRow)
    fprintf('Overwriting %s spot for %s in %s\n', spot, muscle, ID)
    assert(sum(matchingRow) == 1)
    T(matchingRow,:) = newRow;
else
    T = [T; newRow];
end
writetable(T, write_to);
fprintf('Updated %s (%s/%s/%s -> x=%f, y=%f, z=%f, R²=%f) in %s\n', ID, muscle, spot, hemisphere, x, y, z, quality, write_to)
end

