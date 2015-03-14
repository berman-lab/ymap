function [p1_a,p1_b,p1_c, p2_a,p2_b,p2_c] = fit_Gaussian_model_monosomy_2(workingDir, saveName, data,locations,init_width,func_type)
	% attempt to fit a 2-gaussian model to data.

	show = false;
	p1_a = nan;   p1_b = nan;   p1_c = nan;
	p2_a = nan;   p2_b = nan;   p2_c = nan;

	if isnan(data)
		% fitting variables
		return
	end

	% find max height in data.
	datamax = max(data);
	%datamax(data ~= max(datamax)) = [];

	% if maxdata is final bin, then find next highest p
	if (find(data == datamax) == length(data))
		data(data == datamax) = 0;
		datamax = data;
		datamax(data ~= max(datamax)) = [];
	end;

	% a = height; b = location; c = width.
	p1_ai = datamax;   p1_bi = locations(1);   p1_ci = init_width;
	p2_ai = datamax;   p2_bi = locations(2);   p2_ci = init_width;

	initial = [p1_ai,p1_bi,p1_ci,p2_ai,p2_bi,p2_ci];
	options = optimset('Display','off','FunValCheck','on','MaxFunEvals',100000);
	time= 1:length(data);

	[Estimates,~,exitflag] = fminsearch(@fiterror, ...   % function to be fitted.
	                                    initial, ...     % initial values.
	                                    options, ...     % options for fitting algorithm.
	                                    time, ...        % problem-specific parameter 1.
	                                    data, ...        % problem-specific parameter 2.
	                                    func_type, ...   % problem-specific parameter 3.
	                                    locations, ...   % problem-specific parameter 4.
	                                    show ...         % problem-specific parameter 5.
	                            );
	if (exitflag > 0)
		% > 0 : converged to a solution.
	else
		% = 0 : exceeded maximum iterations allowed.
		% < 0 : did not converge to a solution.
		% return last best estimate anyhow.
	end;
	p1_a         = abs(Estimates(1));
	p1_b         = abs(Estimates(2));
	p1_c         = abs(Estimates(3));
	p2_a         = abs(Estimates(4));
	p2_b         = abs(Estimates(5));
	p2_c         = abs(Estimates(6));
end

function sse = fiterror(params,time,data,func_type,locations,show)
	p1_a         = abs(params(1));   % height.
	p1_b         = abs(params(2));   % location.
	p1_c         = abs(params(3));   % width.
	p2_a         = abs(params(4));   % height.
	p2_b         = abs(params(5));   % location.
	p2_c         = abs(params(6));   % width.
	if (p1_c == 0); p1_c = 0.001; end
	if (p2_c == 0); p2_c = 0.001; end
	if (p1_c < 2);   p1_c = 2;   end;
	if (p2_c < 2);   p2_c = 2;   end;
	p1_fit = p1_a*exp(-0.5*((time-p1_b)./p1_c).^2);
	p2_fit = p2_a*exp(-0.5*((time-p2_b)./p2_c).^2);
	fitted = p1_fit+p2_fit;

if (show ~= 0)
%----------------------------------------------------------------------
% show fitting in process.
figure(show);
% show data being fit.
plot(data,'x-','color',[0.75 0.75 1]);
hold on;
title('monosomy');
% show fit lines.
plot(p1_fit,'-','color',[0 0.75 0.75],'lineWidth',2);
plot(p2_fit,'-','color',[0 0.75 0.75],'lineWidth',2);
plot(fitted,'-','color',[0 0.50 0.50],'lineWidth',2);
hold off;
%----------------------------------------------------------------------
end;

	width = 0.5;
	switch(func_type)
		case 'cubic'
			Error_Vector = (fitted).^2 - (data).^2;
			sse  = sum(abs(Error_Vector));
		case 'linear'
			Error_Vector = (fitted) - (data);
			sse  = sum(Error_Vector.^2);
		case 'log'
			Error_Vector = log(fitted) - log(data);
			sse  = sum(abs(Error_Vector));
		case 'fcs'
			Error_Vector = (fitted) - (data);
			%Error_Vector(1:round(G1_b*(1-width))) = 0;
			%Error_Vector(round(G1_b*(1+width)):end) = 0;
			sse  = sum(Error_Vector.^2);
		otherwise
			error('Error: choice for fitting not implemented yet!');
			sse  = 1;
	end;
end
