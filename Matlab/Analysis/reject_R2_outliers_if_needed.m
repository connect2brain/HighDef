function [data] = reject_R2_outliers_if_needed(data)
% Load the config file:
fname = 'config.json'; 
fid = fopen(fname); 
cfg = jsondecode(char(fread(fid,inf)')); 
fclose(fid); 

if max(data) < cfg.R2_rejection.only_apply_if_best_R2_is_lower_than
    % OUTLIER REJECTION!
    Q3 = quantile(log(data), 0.75);
    IQR = iqr(log(data));
    cutoff = cfg.R2_rejection.Cutoff_IQR_factor*IQR + Q3;
    mask = log(data) > cutoff;
    fprintf("  . Rejecting %d outliers (Q3 = %0.8f, IQR = %0.8f => Cutoff = %0.8f; element 13547 has R² = %0.8f)\n", sum(mask), Q3, IQR, cutoff, log(data(13547)))
    data(mask) = 0;
end

end

