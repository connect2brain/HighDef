spots = readtable('B:/Projects/2023-01 HighDef/Results/Evaluation/spot_locations.csv');
subjects = unique(spots.Subject)';

encoding = [];
encoding.APB = [0.2 0 1];
encoding.FDI = [1 0 0];
encoding.ADM = [0 0.5 0.9];
encoding.R2.MarkerSize = 20;
encoding.R2.LineWidth = 3;


hemisphereSpots = spots(strcmpi(spots.Hemisphere, 'L'),:);

hold on
for iSubject = 1:length(subjects)
    subject = subjects{iSubject};
    subjectSpots = hemisphereSpots(strcmpi(hemisphereSpots.Subject, subject),:);
    
    for muscle = {'FDI', 'APB', 'ADM'}
        fprintf('%s - %s\n', subject, muscle{:})
        muscleSpots = subjectSpots(strcmpi(subjectSpots.Muscle, muscle{:}),:);
        hotMask = strcmpi(muscleSpots.Spot, 'hot');
        coldMask = strcmpi(muscleSpots.Spot, 'cold');
        iniMask = strcmpi(muscleSpots.Session, 'ini');
        
        % Connect ini and full, if present.
        if sum(iniMask) > 0 && sum(~iniMask) > 0
            plot3([muscleSpots.X(hotMask & iniMask) muscleSpots.X(hotMask & ~iniMask)], ...
                  [muscleSpots.Y(hotMask & iniMask) muscleSpots.Y(hotMask & ~iniMask)], ...
                  [muscleSpots.Z(hotMask & iniMask) muscleSpots.Z(hotMask & ~iniMask)], 'k')
        end

        % Connect Hotspot and Coldspot
        plot3([muscleSpots.X(hotMask & ~iniMask) muscleSpots.X(coldMask & ~iniMask)], ...
              [muscleSpots.Y(hotMask & ~iniMask) muscleSpots.Y(coldMask & ~iniMask)], ...
              [muscleSpots.Z(hotMask & ~iniMask) muscleSpots.Z(coldMask & ~iniMask)], Color=encoding.(muscle{:}))

        % Scale linewidth with R²
        scatterScaleWithR2(muscleSpots(hotMask & iniMask,:), muscle{:}, '^', encoding)
        scatterScaleWithR2(muscleSpots(hotMask & ~iniMask,:), muscle{:}, '^', encoding)
        scatterScaleWithR2(muscleSpots(~hotMask & iniMask,:), muscle{:}, 'v', encoding)
        scatterScaleWithR2(muscleSpots(~hotMask & ~iniMask,:), muscle{:}, 'v', encoding)
    end

    maskHotMain = strcmpi(subjectSpots.Spot, 'hot') & ~strcmpi(subjectSpots.Session, 'ini');
    plot3(subjectSpots.X(maskHotMain), subjectSpots.Y(maskHotMain), subjectSpots.Z(maskHotMain), ':k')

end

axis square



function scatterScaleWithR2(tab, muscle, marker, encoding)
for iRow = 1:size(tab,1)
    plot3(tab(iRow,:), 'X', 'Y', 'Z', Color=encoding.(muscle), Marker=marker, LineWidth=encoding.R2.LineWidth*tab(iRow,:).R2, MarkerSize=encoding.R2.MarkerSize*(1-tab(iRow,:).R2))
end
end