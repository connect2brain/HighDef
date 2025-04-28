
make_sigmoid4 = @(y_0, amp, r, x_0) (@(x) y_0 + (amp-y_0)./(1+exp(-r.*(x-x_0))));


sigmoid4 = make_sigmoid4(2, 5, 0.5, 2);
noise_sigmoid = make_sigmoid4(0.05, 2, 0.4, 4);

n = 300;

x = rand(n, 1) * 30 - 20;
y_raw = sigmoid4(x);


[fitresult_raw, gof] = regress(x,y_raw);
R2_raw = gof.rsquare;
nexttile()
plot(x,y_raw,'kx')
hold on
plot(fitresult_raw);




noise_homoskedastic = randn(n, 1) * 1;

y_homoskedastic = y_raw + noise_homoskedastic;

[fitresult_homoskedastic, gof] = regress(x,y_homoskedastic);
R2_homoskedastic = gof.rsquare;
nexttile()
plot(x,y_homoskedastic,'kx')
hold on
plot(fitresult_homoskedastic, 'b');
plot(fitresult_raw, 'r');
legend(["", "homoskedastic", "raw"])




noise_heteroskedastic = randn(n,1) .* noise_sigmoid(x);
noise_heteroskedastic = noise_heteroskedastic .* (std(noise_homoskedastic)/std(noise_heteroskedastic));

y_heteroskedastic = y_raw + noise_heteroskedastic;

[fitresult_heteroskedastic, gof] = regress(x,y_heteroskedastic);
R2_heteroskedastic = gof.rsquare;
nexttile()
plot(x,y_heteroskedastic,'kx')
hold on
plot(fitresult_heteroskedastic, 'g');
hold on
plot(fitresult_homoskedastic, 'b--');
hold on
plot(fitresult_raw, 'r');
legend(["", "heteroskedastic", "homoskedastic", "raw"])

fprintf("Homoskedastic:   R² = %.5f (like SIHI)\nHeteroskedastic: R² = %.5f (like MEP)\n", R2_homoskedastic, R2_heteroskedastic)



function [fitresult, gof] = regress(x,y)

[xData, yData] = prepareCurveData(x, y);
% Set up fittype and options.
ft = fittype( 'a/(1+exp(-b*(x-s))) + o', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [max(y) 1 mean(x) min(y)];
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
end