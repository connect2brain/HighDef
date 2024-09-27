function [] = runBlock(bd, times, channel, markers, leftStimulator, rightStimulator, intensityLeft, intensityRight)
% RUNBLOCK run one stimulation block, sending pulses (specified by times, channel, markers) at
% specified intensities, using the given stimulators.
%
% Arguments:
%  bd: Bossdevice object
%  times: a column vector of the times of the stimuli to be delivered, from the start of the block
%  channel: a column vector of the channel from which to deliver the pulse at the corresponding time
%  markers: a column vector of the marker to be set in the NeurOne-recording for the pulse at the
%  corresponding time
%  leftStimulator: MAGIC-toolbox object representing and connected to the stimulator connected to
%  the coil placed on the left hemisphere
%  rightStimulator: MAGIC-toolbox object representing the other stimulator
%  intensityLeft: stimulation intensity to be used by the leftStimulator (in %MSO)
%  intensityRight: stimulation intensity to be used by the rightStimulator (in %MSO)


%arguments
%    bd;
%    times (:,1) {mustBeNumeric};
%    channel (:,1) {mustBeNumeric};
%    markers (:,1);
%    leftStimulator;
%    rightStimulator;
%    intensityLeft (1,1) int {mustBePositive, mustBeLessThanOrEqual(intensityLeft, 100)};
%    intensityRight (1,1) int {mustBePositive, mustBeLessThanOrEqual(intensityRight, 100)};
%end
    bd.configure_time_port_marker([times, channel, markers]);

    leftStimulator.arm;
    leftStimulator.setAmplitude(intensityLeft);
    rightStimulator.arm;
    rightStimulator.setAmplitude(intensityRight);

    pause(1)
    bd.manualTrigger;
end