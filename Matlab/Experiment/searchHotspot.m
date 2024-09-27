function searchHotspot(stimulator, bd, port)
stimulator.arm;
stimulator.setAmplitude(65); 

fprintf('\nHotspot search. Observe EMG, change  Press Ctrl+C to abort\n')

% Initializing Hotspot search with randomly triggered pulses (1.7 - 2.7 seconds) to avoid modulation effect
while(true), bd.sendPulse(port), fprintf('.'), pause(1.7+rand()), end
end