function missing = coil_positions_missing(transform, distance_threshold)
% Exact matches to either identity matrix (rotation), or, more typically:
% to location 0,0,0 indicate that the coil was not seen.
arguments
    transform (4,4,:) {mustBeNumeric};
    distance_threshold (1,1) {mustBeNumeric} = 200; % in mm
end
    missing = squeeze(all(transform(1:3,1:3,:) == eye(3), [1 2]))' | squeeze(all(transform(1:3,4,:) == zeros(3,1), [1 2]))';
    
    too_far = sqrt(sum((squeeze(transform(1:3,4,:)) - [0 0 40]').^2, 1)) > distance_threshold;
    missing = missing | too_far;
end