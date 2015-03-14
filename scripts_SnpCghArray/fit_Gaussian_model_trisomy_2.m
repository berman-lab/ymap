function [p1_a,p1_b,p1_c, p2_a,p2_b,p2_c, p3_a,p3_b,p3_c, p4_a,p4_b,p4_c, skew_factor] = fit_Gaussian_model_trisomy_2(data,locations,init_width,fraction,skew_factor,func_type,show)
% attempt to fit a single-gaussian model to data.
%[G1_a, G1_b, G1_c, G2_a, G2_b, G2_c, S_a, S_c] = GaussianModel_G1SG2(tet_control,parameter,'fcs1','');
	p1_a = nan;   p1_b = nan;   p1_c = nan;
	p2_a = nan;   p2_b = nan;   p2_c = nan;
	p3_a = nan;   p3_b = nan;   p3_c = nan;
	p4_a = nan;   p4_b = nan;   p4_c = nan;
    
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
	p4_ai = datamax;   p4_bi = locations(3);   p4_ci = init_width;
   
	initial = [p1_ai,p1_ci,  p2_ai,  p3_ai,  p4_ai,  skew_factor];
	options = optimset('Display','off','FunValCheck','on','MaxFunEvals',100000);
	time    = 1:length(data);

	if (data == zeros(1,length(data)))  % curve fittings don't work with no data, curves instead are flat.
		p1_a = 1;
		p1_b = locations(1);
		p1_c = 5;
		p2_a = 1;
		p2_b = locations(2);
		p2_c = p2_a/p1_a*p1_c;             % peak width scales with peak height.
		p3_a = 1;
		p3_b = locations(3);
		p3_c = p3_a/p1_a*p1_c;             % peak width scales with peak height.
		p4_a = 1;
		p4_b = locations(4);
		p4_c = p4_a/p1_a*p1_c;             % peak width scales with peak height.
		skew_factor = 1;
	else
		[Estimates,~,exitflag] = fminsearch(@fiterror, ...   % function to be fitted.
		                                    initial, ...     % initial values.
		                                    options, ...     % options for fitting algorithm.
		                                    time, ...        % problem-specific parameter 1.
		                                    data, ...        % problem-specific parameter 2.
		                                    func_type, ...   % problem-specific parameter 3.
		                                    locations, ...   % problem-specific parameter 4.
		                                    show, ...        % problem-specific parameter 5.
		                                    fraction ...     % problem-specific parameter 6.
		                            );
		if (exitflag > 0)
			% > 0 : converged to a solution.
		else
			% = 0 : exceeded maximum iterations allowed.
			% < 0 : did not converge to a solution.
			% return last best estimate anyhow.
		end;
		p1_a = abs(Estimates(1));
		p1_b = locations(1);
		p1_c = abs(Estimates(2));
		if (p1_c < 2);   p1_c = 2;   end;
		p2_a = abs(Estimates(3));
		p2_b = locations(2);
		p2_c = p2_a/p1_a*p1_c;             % peak width scales with peak height.
		p3_a = abs(Estimates(4));
		p3_b = locations(3);
		p3_c = p3_a/p1_a*p1_c;             % peak width scales with peak height.
		p4_a = abs(Estimates(5));
		p4_b = locations(4);
		p4_c = p4_a/p1_a*p1_c;             % peak width scales with peak height.
		skew_factor = abs(Estimates(6));
	end;
	c1_  = p1_c/2 + p1_c*skew_factor/(100-abs(100-p1_b))/2;
	p1_c = p1_c*p1_c/c1_;
	c2_  = p2_c/2 + p2_c*skew_factor/(100-abs(100-p2_b))/2;
	p2_c = p2_c*p2_c/c2_;
	c3_  = p3_c/2 + p3_c*skew_factor/(100-abs(100-p3_b))/2;
	p3_c = p3_c*p3_c/c3_;
	c4_  = p4_c/2 + p4_c*skew_factor/(100-abs(100-p4_b))/2;
	p4_c = p4_c*p4_c/c4_;
end

function sse = fiterror(params,time,data,func_type,locations,show,fraction)
	p1_a = abs(params(1));
	p1_b = locations(1);
	p1_c = abs(params(2));
	if (p1_c < 2);   p1_c = 2;   end;
	p2_a = abs(params(3));
	p2_b = locations(2);
	p2_c = p2_a/p1_a*p1_c;             % peak width scales with peak height.
	p3_a = abs(params(4));
	p3_b = locations(3);
	p3_c = p3_a/p1_a*p1_c;             % peak width scales with peak height.
	p4_a = abs(params(5));
	p4_b = locations(4);
	p4_c = p4_a/p1_a*p1_c;             % peak width scales with peak height.
	skew_factor = abs(params(6));
	time1_1 = 1:floor(p1_b);
	time1_2 = ceil(p1_b):200;
	if (time1_1(end) == time1_2(1));    time1_1(end) = [];  end;
	time2_1 = 1:floor(p2_b);
	time2_2 = ceil(p2_b):200;
	if (time2_1(end) == time2_2(1));    time2_1(end) = [];  end;
	time3_1 = 1:floor(p3_b);
	time3_2 = ceil(p3_b):200;
	if (time3_1(end) == time3_2(1));    time3_2(1) = [];    end;
	time4_1 = 1:floor(p4_b);
	time4_2 = ceil(p4_b):200;
	if (time4_1(end) == time4_2(1));    time4_2(1) = [];    end;
	if (skew_factor/(100-abs(100-p1_b)) < 1)
		skew_factor = (100-abs(100-p1_b));
	end;
	c1_  = p1_c/2 + p1_c*skew_factor/(100-abs(100-p1_b))/2;
	p1_c = p1_c*p1_c/c1_;
	c2_  = p2_c/2 + p2_c*skew_factor/(100-abs(100-p2_b))/2;
	p2_c = p2_c*p2_c/c2_;
	c3_  = p3_c/2 + p3_c*skew_factor/(100-abs(100-p3_b))/2;
	p3_c = p3_c*p3_c/c3_;
	c4_  = p4_c/2 + p4_c*skew_factor/(100-abs(100-p4_b))/2;
	p4_c = p4_c*p4_c/c4_;
	p1_fit_L = p1_a*exp(-0.5*((time1_1-p1_b)./p1_c).^2);
	p1_fit_R = p1_a*exp(-0.5*((time1_2-p1_b)./p1_c/(skew_factor/(100-abs(100-p1_b))) ).^2);
	p2_fit_L = p2_a*exp(-0.5*((time2_1-p2_b)./p2_c).^2);
	p2_fit_R = p2_a*exp(-0.5*((time2_2-p2_b)./p2_c/(skew_factor/(100-abs(100-p2_b))) ).^2);
	p3_fit_L = p3_a*exp(-0.5*((time3_1-p3_b)./p3_c/(skew_factor/(100-abs(100-p3_b))) ).^2);
	p3_fit_R = p3_a*exp(-0.5*((time3_2-p3_b)./p3_c).^2);
	p4_fit_L = p4_a*exp(-0.5*((time4_1-p4_b)./p4_c/(skew_factor/(100-abs(100-p4_b))) ).^2);
	p4_fit_R = p4_a*exp(-0.5*((time4_2-p4_b)./p4_c).^2);
	p1_fit = [p1_fit_L p1_fit_R];
	p2_fit = [p2_fit_L p2_fit_R];
	p3_fit = [p3_fit_L p3_fit_R];
	p4_fit = [p4_fit_L p4_fit_R];
	fitted = p1_fit+p2_fit+p3_fit+p4_fit;

if (show ~= 0)
%----------------------------------------------------------------------
% show fitting in process.
figure(show);
% show data being fit.
plot(data,'x-','color',[0.75 0.75 1]);
hold on;
title('trisomy');
% show fit lines.
plot(p1_fit,'-','color',[0 0.75 0.75],'lineWidth',2);
plot(p2_fit,'-','color',[0 0.75 0.75],'lineWidth',2);
plot(p3_fit,'-','color',[0 0.75 0.75],'lineWidth',2);
plot(p4_fit,'-','color',[0 0.75 0.75],'lineWidth',2);
plot(fitted,'-','color',[0 0.50 0.50],'lineWidth',2);
hold off;
%----------------------------------------------------------------------
end;

	width = 0.5;
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
		case 'fcs'
			Error_Vector = (fitted) - (data);
			%Error_Vector(1:round(G1_b*(1-width))) = 0;
			%Error_Vector(round(G1_b*(1+width)):end) = 0;
			sse          = sum(Error_Vector.^2);
		otherwise
			error('Error: choice for fitting not implemented yet!');
			sse          = 1;            
	end;
end
