function mep = extractMEP(emg,timeAxis,mepTimeWindow,baselineTimeWindow)

% select MEP
selMEP = mepTimeWindow(1) <= timeAxis & timeAxis <= mepTimeWindow(2);
% select baseline
selBase = baselineTimeWindow(1) <= timeAxis & timeAxis <= baselineTimeWindow(2);

% remove baseline mean
emg = emg - mean(emg(selBase));

% MEP max and min
mepmin = min(emg(selMEP));
mepmax = max(emg(selMEP));

% Time of min and max amplitudes
mepMinTime = timeAxis(find(selMEP,1)-1+find(emg(selMEP)==mepmin,1));
mepMaxTime = timeAxis(find(selMEP,1)-1+find(emg(selMEP)==mepmax,1));

mep.minmax = [min(emg(selMEP)) max(emg(selMEP))];
mep.minmaxTime = [mepMinTime mepMaxTime];

% MEP amplitude
mep.amplitude = mepmax-mepmin;

% MEP baseline amplitude range
mep.baseline = [min(emg(selBase)) max(emg(selBase))];

end