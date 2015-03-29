function [p1_a,p1_b,p1_c, p2_a,p2_b,p2_c, p3_a,p3_b,p3_c, p4_a,p4_b,p4_c, p5_a,p5_b,p5_c, p6_a,p6_b,p6_c, p7_a,p7_b,p7_c, p8_a,p8_b,p8_c, p9_a,p9_b,p9_c] = fit_Gaussian_model_octasomy_2(workingDir, saveName, data,locations,init_width,func_type)
	% attempt to fit a 9-gaussian model to data.

	show = false;
	p1_a = nan;   p1_b = nan;   p1_c = nan;
	p2_a = nan;   p2_b = nan;   p2_c = nan;
	p3_a = nan;   p3_b = nan;   p3_c = nan;
	p4_a = nan;   p4_b = nan;   p4_c = nan;
	p5_a = nan;   p5_b = nan;   p5_c = nan;
	p6_a = nan;   p6_b = nan;   p6_c = nan;
	p7_a = nan;   p7_b = nan;   p7_c = nan;
	p8_a = nan;   p8_b = nan;   p8_c = nan;
	p9_a = nan;   p9_b = nan;   p9_c = nan;

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
	p3_ai = datamax;   p3_bi = locations(3);   p3_ci = init_width;
	p4_ai = datamax;   p4_bi = locations(4);   p4_ci = init_width;
	p5_ai = datamax;   p5_bi = locations(5);   p5_ci = init_width;
	p6_ai = datamax;   p6_bi = locations(6);   p6_ci = init_width;
	p7_ai = datamax;   p7_bi = locations(7);   p7_ci = init_width;
	p8_ai = datamax;   p8_bi = locations(8);   p8_ci = init_width;
	p9_ai = datamax;   p9_bi = locations(9);   p9_ci = init_width;

	initial = [p1_ai,p1_ci,p2_ai,p2_ci,p3_ai,p3_ci,p4_ai,p4_ci,p5_ai,p5_ci,p6_ai,p6_ci,p7_ai,p7_ci,p8_ai,p8_ci,p9_ai,p9_ci];
	options = optimset('Display','off','FunValCheck','on','MaxFunEvals',100000);
	time= 1:length(data);

	[Estimates,~,exitflag] = fminsearch(@fiterror, ...   % function to be fitted.
	                                    initial, ...     % initial values.
	                                    options, ...     % options for fitting algorithm.
	                                    time, ...        % problem-specific parameter 1.
	                                    data, ...        % problem-specific parameter 2.
	                                    func_type, ...   % problem-specific parameter 3.
	                                    locations, ...   % problem-specific parameter 4.
	                                    show ...     % problem-specific parameter 5.
	                            );
	if (exitflag > 0)
		% > 0 : converged to a solution.
	else
		% = 0 : exceeded maximum iterations allowed.
		% < 0 : did not converge to a solution.
		% return last best estimate anyhow.
	end;
	p1_a         = abs(Estimates(1));
	p1_b         = locations(1);
	p1_c         = abs(Estimates(2));
	p2_a         = abs(Estimates(3));
	p2_b         = locations(2);
	p2_c         = abs(Estimates(4));
	p3_a         = abs(Estimates(5));
	p3_b         = locations(3);
	p3_c         = abs(Estimates(6));
	p4_a         = abs(Estimates(7));
	p4_b         = locations(4);
	p4_c         = abs(Estimates(8));
	p5_a         = abs(Estimates(9));
	p5_b         = locations(5);
	p5_c         = abs(Estimates(10));
	p6_a         = abs(Estimates(11));
	p6_b         = locations(6);
	p6_c         = abs(Estimates(12));
	p7_a         = abs(Estimates(13));
	p7_b         = locations(7);
	p7_c         = abs(Estimates(14));
	p8_a         = abs(Estimates(15));
	p8_b         = locations(8);
	p8_c         = abs(Estimates(16));
	p9_a         = abs(Estimates(17));
	p9_b         = locations(9);
	p9_c         = abs(Estimates(18));
end

function sse = fiterror(params,time,data,func_type,locations,show)
	p1_a         = abs(params(1));   % height.
	p1_b         = locations(1);     % location.
	p1_c         = abs(params(2));   % width.
	p2_a         = abs(params(3));   % height.
	p2_b         = locations(2);     % location.
	p2_c         = abs(params(4));   % width.
	p3_a         = abs(params(5));   % height.
	p3_b         = locations(3);     % location.
	p3_c         = abs(params(6));   % width.
	p4_a         = abs(params(7));   % height.
	p4_b         = locations(4);     % location.
	p4_c         = abs(params(8));   % width.
	p5_a         = abs(params(9));   % height.
	p5_b         = locations(5);     % location.
	p5_c         = abs(params(10));  % width.
	p6_a         = abs(params(11));  % height.
	p6_b         = locations(6);     % location.
	p6_c         = abs(params(12));  % width.
	p7_a         = abs(params(13));  % height.
	p7_b         = locations(7);     % location.
	p7_c         = abs(params(14));  % width.
	p8_a         = abs(params(15));  % height.
	p8_b         = locations(8);     % location.
	p8_c         = abs(params(16));  % width.
	p9_a         = abs(params(17));  % height.
	p9_b         = locations(9);     % location.
	p9_c         = abs(params(18));  % width.
	if (p1_c == 0); p1_c = 0.001; end;
	if (p2_c == 0); p2_c = 0.001; end;
	if (p3_c == 0); p3_c = 0.001; end;
	if (p4_c == 0); p4_c = 0.001; end;
	if (p5_c == 0); p5_c = 0.001; end;
	if (p6_c == 0); p6_c = 0.001; end;
	if (p7_c == 0); p7_c = 0.001; end;
	if (p8_c == 0); p8_c = 0.001; end;
	if (p9_c == 0); p9_c = 0.001; end;
	if (p1_c < 2);   p1_c = 2;   end;
	if (p2_c < 2);   p2_c = 2;   end;
	if (p3_c < 2);   p3_c = 2;   end;
	if (p4_c < 2);   p4_c = 2;   end;
	if (p5_c < 2);   p5_c = 2;   end;
	if (p6_c < 2);   p6_c = 2;   end;
	if (p7_c < 2);   p7_c = 2;   end;
	if (p8_c < 2);   p8_c = 2;   end;
	if (p9_c < 2);   p9_c = 2;   end;
	p1_fit = p1_a*exp(-0.5*((time-p1_b)./p1_c).^2);
	p2_fit = p2_a*exp(-0.5*((time-p2_b)./p2_c).^2);
	p3_fit = p3_a*exp(-0.5*((time-p3_b)./p3_c).^2);
	p4_fit = p4_a*exp(-0.5*((time-p4_b)./p4_c).^2);
	p5_fit = p5_a*exp(-0.5*((time-p5_b)./p5_c).^2);
	p6_fit = p6_a*exp(-0.5*((time-p6_b)./p6_c).^2);
	p7_fit = p7_a*exp(-0.5*((time-p7_b)./p7_c).^2);
	p8_fit = p8_a*exp(-0.5*((time-p8_b)./p8_c).^2);
	p9_fit = p9_a*exp(-0.5*((time-p9_b)./p9_c).^2);
	fitted = p1_fit+p2_fit+p3_fit+p4_fit+p5_fit+p6_fit+p7_fit+p8_fit+p9_fit;

if (show ~= 0)
%----------------------------------------------------------------------
% show fitting in process.
figure(show);
% show data being fit.
plot(data,'x-','color',[0.75 0.75 1]);
hold on;
title('hexasomy');
% show fit lines.
plot(p1_fit,'-','color',[0 0.75 0.75],'lineWidth',2);
plot(p2_fit,'-','color',[0 0.75 0.75],'lineWidth',2);
plot(p3_fit,'-','color',[0 0.75 0.75],'lineWidth',2);
plot(p4_fit,'-','color',[0 0.75 0.75],'lineWidth',2);
plot(p5_fit,'-','color',[0 0.75 0.75],'lineWidth',2);
plot(p6_fit,'-','color',[0 0.75 0.75],'lineWidth',2);
plot(p7_fit,'-','color',[0 0.75 0.75],'lineWidth',2);
plot(p8_fit,'-','color',[0 0.75 0.75],'lineWidth',2);
plot(p9_fit,'-','color',[0 0.75 0.75],'lineWidth',2);
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
