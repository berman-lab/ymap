function [G1_a, G1_b, G1_c, G2_a, G2_b, G2_c] = fit2GaussianModel2(dataX,dataY,func_type)
% Fit a 2 guassian mixed model to input data.
%    G1 and G2 are vertical lines...  broadened by noise.

time      = dataX;   % 1:length(dataX);
histogram = dataY;   % data.histogram{parameter};
G1_a      = nan;
G1_b      = nan;
G1_c      = nan;
G2_a      = nan;
G2_b      = nan;
G2_c      = nan;

% find max height in data.
[maxData,maxIndex] = max(dataY);

% a = height.
% b = location.
% c = width.

G1_ai   = dataY(1);
G1_bi   = 0;
G1_ci   = 5;
G2_ai   = maxData;
G2_bi   = maxIndex-1;
G2_ci   = 5;

initial = [G1_ai, G1_ci, G2_ai, G2_bi, G2_ci];
options = optimset('Display','off','FunValCheck','on','MaxFunEvals',500);
data    = histogram;

[Estimates,ignore,exitflag] = fminsearch(@fiterror, ...   % function to be fitted.
	initial, ...      % initial x-value.
	options, ...      % options for fitting algorithm.
	time, ...         % problem-specific parameter 1.
	data, ...         % problem-specific parameter 2.
	func_type ...     % problem-specific parameter 3.
	);
if (exitflag > 0)
	% > 0 : converged to a solution.
	G1_a = abs(Estimates(1));
	G1_b = 0;
	G1_c = abs(Estimates(2));
	G2_a = abs(Estimates(3));
	G2_b = abs(Estimates(4));
	G2_c = abs(Estimates(5));
else
	% = 0 : exceeded maximum iterations allowed.
	% < 0 : did not converge to a solution.
	% return last estimage in these cases...
	G1_a = abs(Estimates(1));
	G1_b = 0;
	G1_c = abs(Estimates(2));
	G2_a = abs(Estimates(3));
	G2_b = abs(Estimates(4));
	G2_c = abs(Estimates(5));
end;
end

function sse = fiterror(params,time,data,func_type)
G1_a = abs(params(1));       % G1_height.
G1_b = 0;                    % G1_location.
G1_c = abs(params(2));       % G1_width.
G2_a = abs(params(3));       % G2_height.
G2_b = abs(params(4));       % G2_location.
G2_c = abs(params(5));       % G2_width.

G1_fit = G1_a*exp(-0.5*((time-G1_b)./G1_c).^2);
G2_fit = G2_a*exp(-0.5*((time-G2_b)./G2_c).^2);
fitted = G1_fit + G2_fit;

switch(func_type)
    case 'cubic'
        Error_Vector = (fitted).^2 - (data).^2;
        sse          = sum(abs(Error_Vector));
    case 'linear'
        Error_Vector = (fitted) - (data);
        sse          = sum(Error_Vector.^2);
    case 'log'
        Error_Vector = log(fitted) - log(data);
        sse          = sum(abs(Error_Vector));
    otherwise
        error('Error: choice for fitting not implemented yet!');
end;
end
