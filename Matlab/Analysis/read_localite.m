function [time, transform, coordinate_space] = read_localite(path, coil, timestamp)
localiteFile = sprintf('%s/TMSTrigger/TriggerMarkers_Coil%d_%s.xml', path, coil, timestamp);
content = parseXML(localiteFile);
coordinate_space = content.Attributes(arrayfun(@(s) strcmpi(s.Name, 'coordinateSpace'), content.Attributes)).Value;
[time, transform] = extractTransformationMatrices(content);
[time, reordering] = sort(time);
transform = transform(:,:,reordering);
end