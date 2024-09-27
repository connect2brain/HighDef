function [intensities, meps, fig, fitresult] = runIOCurve(bd, stimulator, outChannel, rmt, maxNumberUniqueIntensities, channelAPB, channelFDI) 
% RUNIOCURVE Generates an input-output (IO) curve, by running the needed measurement using the
% given stimulator.
% Arguments:
%  bd: Bossdevice object
%  stimulator: MAGIC-toolbox stimulator object
%  outChannel: Channel from which the bossdevice will send the pulses
%  rmt: Resting motor threshold estimate. Pulses will be delivered at 100%, 110%, 120%, 130% and
%  140% RMT
%  nPulsesPerIntensity: How many pulses will be delivered per intensity (100%, ..., 140%). At least
%  5 is recommended for stability of the estimate
% channelAPB: Bossdevice-internal channel-index of the ABP-muscle
% channelFDI: Bossdevice-internal channel-index of the FDI-muscle
intensities = unique(round(linspace(0.9*rmt, min([1.4*rmt, 100]), maxNumberUniqueIntensities)));
[meps, fig] = MEP_ioresponse(bd, stimulator, outChannel, intensities, 1, channelAPB, channelFDI);
fprintf('Done collecting MEP-io\n\n')


set(fig, 'Renderer','Painters') %export vectorized
set(fig, 'PaperUnits', 'centimeter', 'PaperSize', [16 16]) % set size
set(fig, 'PaperPositionMode', 'manual', 'PaperUnits', 'normalized', 'PaperPosition',[0 0 1 1]); % fill page
set(fig, 'PaperUnits', 'centimeter') % set back to something other than normalized in order to enable copy to clipboard


FDIMeps = squeeze(meps(:,:,2));
intensityPerMep = intensities;
FDIMeps = FDIMeps(:);
intensityPerMep = intensityPerMep(:);

[xData, yData] = prepareCurveData(intensityPerMep, FDIMeps);
% Set up fittype and options.
ft = fittype( 'a/(1+exp(-b*(x-s)))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [max(FDIMeps) 1 mean(intensities)];
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
scatter(xData, yData); hold on; plot(fitresult)
end