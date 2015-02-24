function [fitted] = fit_Gaussian_model_heterozygousPeak_initial(data,locations,func_type,show, workingDir)
	fitted = [];
	if isnan(data);
		return
	end

	% find max height in data.
	datamax = max(data);
	% if maxdata is final bin, then find next highest p
	if (find(data == datamax) == length(data))
		data(data == datamax) = 0;
		datamax = data;
		datamax(data ~= max(datamax)) = [];
	end;

	% a = height.
	% b = location.
	% c = width.
	pPlus_ai  = datamax;   pPlus_ci  = 15;

	initial = [pPlus_ai,pPlus_ci ...
%	           ,pPlus_ai,pPlus_ci ...
%	           ,pPlus_ai,pPlus_ci ...
%	           ,pPlus_ai,pPlus_ci ...
%	           ,pPlus_ai,pPlus_ci ...
%	           ,pPlus_ai,pPlus_ci ...
%	           ,pPlus_ai,pPlus_ci ...
%	           ,pPlus_ai,pPlus_ci ...
%	           ,pPlus_ai,pPlus_ci ...
	          ];
	options   = optimset('Display','off','FunValCheck','on','MaxFunEvals',1000000);
	time      = 1:length(data);
	func_type = 'test';

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
	p1_a  =  abs(Estimates(1));    p1_b  =  100.5;
	p1_c  =  abs(Estimates(2));
%	p2_a  =  abs(Estimates(3));    p2_b  =  100.5;
%	p2_c  =  abs(Estimates(4));
%	p3_a  =  abs(Estimates(5));    p3_b  =  100.5;
%	p3_c  =  abs(Estimates(6));
%	p4_a  =  abs(Estimates(7));    p4_b  =  100.5;
%	p4_c  =  abs(Estimates(8));
%	p5_a  =  abs(Estimates(9));    p5_b  =  100.5;
%	p5_c  =  abs(Estimates(10));
%	p6_a  =  abs(Estimates(11));   p6_b  =  100.5;
%	p6_c  =  abs(Estimates(12));
%	p7_a  =  abs(Estimates(13));   p7_b  =  100.5;
%	p7_c  =  abs(Estimates(14));
%	p8_a  =  abs(Estimates(15));   p8_b  =  100.5;
%	p8_c  =  abs(Estimates(16));
%	p9_a  =  abs(Estimates(17));   p9_b  =  100.5;
%	p9_c  =  abs(Estimates(18));

	p1_fit  = p1_a *exp(-0.5*((time-p1_b )./p1_c ).^2);
%	p2_fit  = p2_a *exp(-0.5*((time-p2_b )./p2_c ).^2);
%	p3_fit  = p3_a *exp(-0.5*((time-p3_b )./p3_c ).^2);
%	p4_fit  = p4_a *exp(-0.5*((time-p4_b )./p4_c ).^2);
%	p5_fit  = p5_a *exp(-0.5*((time-p5_b )./p5_c ).^2);
%	p6_fit  = p6_a *exp(-0.5*((time-p6_b )./p6_c ).^2);
%	p7_fit  = p7_a *exp(-0.5*((time-p7_b )./p7_c ).^2);
%	p8_fit  = p8_a *exp(-0.5*((time-p8_b )./p8_c ).^2);
%	p9_fit  = p9_a *exp(-0.5*((time-p9_b )./p9_c ).^2);

	fitted = p1_fit;% + p2_fit;% + p3_fit + p4_fit + p5_fit + p6_fit + p7_fit + p8_fit + p9_fit;

	%----------------------------------------------------------------------
	% Generate central peak Gaussian fitting figure.
	fig = figure(123);
	plot(data,'x-','color',[0.75 0.75 1]);
	hold on;
	title('disomy initial'); 
	plot(p1_fit ,'-','color',[0.2   0.8   0.2  ],'lineWidth',1);
%	plot(p2_fit ,'-','color',[0.2   0.8   0.2  ],'lineWidth',1);
%	plot(p3_fit ,'-','color',[0.2   0.8   0.2  ],'lineWidth',1);
%	plot(p4_fit ,'-','color',[0.2   0.8   0.2  ],'lineWidth',1);
%	plot(p5_fit ,'-','color',[0.2   0.8   0.2  ],'lineWidth',1);
%	plot(p6_fit ,'-','color',[0.2   0.8   0.2  ],'lineWidth',1);
%	plot(p7_fit ,'-','color',[0.2   0.8   0.2  ],'lineWidth',1);
%	plot(p8_fit ,'-','color',[0.2   0.8   0.2  ],'lineWidth',1);
%	plot(p9_fit ,'-','color',[0.2   0.8   0.2  ],'lineWidth',1);

	plot(fitted ,'-','color',[0.333 0.333 0.333],'lineWidth',2);
	hold off;
	% saveas(fig, [workingDir 'initGaussianFit.eps'], 'epsc');
	  saveas(fig, [workingDir 'initGaussianFit.png'], 'png');
	delete(fig);
	%----------------------------------------------------------------------
end

function sse = fiterror(params,time,data,func_type,locations,show)
	% a = height.
	% b = location.
	% c = width.
	p1_a  =  abs(params(1));    p1_b  =  100.5;
	p1_c  =  abs(params(2));
%	p2_a  =  abs(params(3));    p2_b  =  100.5;
%	p2_c  =  abs(params(4));
%	p3_a  =  abs(params(5));    p3_b  =  100.5;
%	p3_c  =  abs(params(6));
%	p4_a  =  abs(params(7));    p4_b  =  100.5;
%	p4_c  =  abs(params(8));
%	p5_a  =  abs(params(9));    p5_b  =  100.5;
%	p5_c  =  abs(params(10));
%	p6_a  =  abs(params(11));   p6_b  =  100.5;
%	p6_c  =  abs(params(12));
%	p7_a  =  abs(params(13));   p7_b  =  100.5;
%	p7_c  =  abs(params(14));
%	p8_a  =  abs(params(15));   p8_b  =  100.5;
%	p8_c  =  abs(params(16));
%	p9_a  =  abs(params(17));   p9_b  =  100.5;
%	p9_c  =  abs(params(18));

	p1_fit  = p1_a *exp(-0.5*((time-p1_b )./p1_c ).^2);
%	p2_fit  = p2_a *exp(-0.5*((time-p2_b )./p2_c ).^2);
%	p3_fit  = p3_a *exp(-0.5*((time-p3_b )./p3_c ).^2);
%	p4_fit  = p4_a *exp(-0.5*((time-p4_b )./p4_c ).^2);
%	p5_fit  = p5_a *exp(-0.5*((time-p5_b )./p5_c ).^2);
%	p6_fit  = p6_a *exp(-0.5*((time-p6_b )./p6_c ).^2);
%	p7_fit  = p7_a *exp(-0.5*((time-p7_b )./p7_c ).^2);
%	p8_fit  = p8_a *exp(-0.5*((time-p8_b )./p8_c ).^2);
%	p9_fit  = p9_a *exp(-0.5*((time-p9_b )./p9_c ).^2);

	fitted = p1_fit;% + p2_fit;% + p3_fit + p4_fit + p5_fit + p6_fit + p7_fit + p8_fit + p9_fit;

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
		case 'test'
			Error_Vector = zeros(200,1);
			for i = 1:200
				if (fitted(i) > data(i))  % fitting is above data...  bigger problem.
					Error_Vector(i) = abs(fitted(i) - data(i))^3;
				else                      % fitting is below data...  lesser problem.
					Error_Vector(i) = abs(fitted(i) - data(i));
				end;
			end;
			sse          = sum(abs(Error_Vector));
		otherwise
			error('Error: choice for fitting not implemented yet!');
			sse          = 1;            
	end;
end
