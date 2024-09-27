function out = createTimeAxis(timeWindow,samplingFrequency)
% timeWindow in ms
% samplingFrequency in Hz
%
% returns: timeaxis in ms

out = (timeWindow(1)+1000/samplingFrequency):1000/samplingFrequency:timeWindow(2); % time axis (ms)