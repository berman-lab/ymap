function [x_peak,actual_cutoffs,mostLikelyGaussians] = FindGaussianCutoffs_3(workingDir,saveName, chromosome,segment,copyNum, smoothed_Histogram, MakeFigure);

monosomy_peaks  = [0, 1]*199+1;
disomy_peaks    = [0, 1/2, 1]*199+1;
trisomy_peaks   = [0, 1/3, 2/3, 1]*199+1;
tetrasomy_peaks = [0, 1/4, 2/4, 3/4, 1]*199+1;
pentasomy_peaks = [0, 1/5, 2/5, 3/5, 4/5, 1]*199+1;
hexasomy_peaks  = [0, 1/6, 2/6, 3/6, 4/6, 5/6, 1]*199+1;
heptasomy_peaks = [0, 1/7, 2/7, 3/7, 4/7, 5/7, 6/7, 1]*199+1;
octasomy_peaks  = [0, 1/8, 2/8, 3/8, 4/8, 5/8, 6/8, 7/8, 1]*199+1;
nonasomy_peaks  = [0, 1/9, 2/9, 3/9, 4/9, 5/9, 6/9, 7/9, 8/9, 1]*199+1;

%% Calculation of Gaussians against per chromosome data.
% Fits Gaussians to real data per chromomsome, then determines equal probability cutoffs between them.
sigma = 5;
%% FindGaussianCutoffs Finds cutoffs as intersections of Gaussians, fit to the data at each peak location.
ErrorType      = 'linear';
% Define range of fit curves.
range = 1:200;

% Generate figure.
if (MakeFigure == true)
	fig = figure(1);
	hold on;
end;
if (copyNum == 0)
	G                   = [];
	list                = [];
	x_peak              = [];
	actual_cutoffs      = [];
	mostLikelyGaussians = [];
elseif (copyNum == 1)
	G                   = [];
	[G{1}.a,G{1}.b,G{1}.c, G{2}.a,G{2}.b,G{2}.c] = ...
		fit_Gaussian_model_monosomy_2(workingDir,saveName, smoothed_Histogram,monosomy_peaks,sigma,ErrorType);
	[list]              = FindHighestGaussian_2(G);
	actual_cutoffs      = [];
	mostLikelyGaussians = [];
	for i = 1:199
		if (list(i) ~= list(i+1))   % we've found a boundary.
			actual_cutoffs = [actual_cutoffs (i+0.5)];
			% actual_cutoffs = [actual_cutoffs FindGaussianCrossover_2(G{list(i)},G{list(i+1)},i)];
			mostLikelyGaussians = [mostLikelyGaussians list(i)];
		end;
	end;
	mostLikelyGaussians = [mostLikelyGaussians list(200)];
	x_peak              = [];
	for i = 1:2; x_peak(i) = G{i}.b; end;

	% Construct curve
	fit_curve_1 = G{1}.a*exp(-0.5*((range-G{1}.b)./G{1}.c).^2);
	fit_curve_2 = G{2}.a*exp(-0.5*((range-G{2}.b)./G{2}.c).^2);
	plot(fit_curve_1,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_2,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	fit_curve_tot = fit_curve_1+fit_curve_2;
elseif (copyNum == 2)
	G                   = [];
	[G{1}.a,G{1}.b,G{1}.c, G{2}.a,G{2}.b,G{2}.c, G{3}.a,G{3}.b,G{3}.c] = ...
		fit_Gaussian_model_disomy_2(workingDir,saveName, smoothed_Histogram,disomy_peaks,sigma,ErrorType);
	[list]              = FindHighestGaussian_2(G);
	actual_cutoffs      = [];
	mostLikelyGaussians = [];
	for i = 1:199
		if (list(i) ~= list(i+1))   % we've found a boundary.
			actual_cutoffs = [actual_cutoffs (i+0.5)];
			% actual_cutoffs = [actual_cutoffs FindGaussianCrossover_2(G{list(i)},G{list(i+1)},i)];
			mostLikelyGaussians = [mostLikelyGaussians list(i)];
		end;
	end;
	mostLikelyGaussians = [mostLikelyGaussians list(200)];
	x_peak              = [];
	for i = 1:3; x_peak(i) = G{i}.b; end;

	% Construct curve
	fit_curve_1 = G{1}.a*exp(-0.5*((range-G{1}.b)./G{1}.c).^2);
	fit_curve_2 = G{2}.a*exp(-0.5*((range-G{2}.b)./G{2}.c).^2);
	fit_curve_3 = G{3}.a*exp(-0.5*((range-G{3}.b)./G{3}.c).^2);
	plot(fit_curve_1,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_2,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_3,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	fit_curve_tot = fit_curve_1+fit_curve_2+fit_curve_3;
elseif (copyNum == 3)
	G                   = [];
	[G{1}.a,G{1}.b,G{1}.c, G{2}.a,G{2}.b,G{2}.c, G{3}.a,G{3}.b,G{3}.c, G{4}.a,G{4}.b,G{4}.c] = ...
		fit_Gaussian_model_trisomy_2(workingDir,saveName, smoothed_Histogram,trisomy_peaks,sigma,ErrorType);
	list                = FindHighestGaussian_2(G);
	actual_cutoffs      = [];
	mostLikelyGaussians = [];
	for i = 1:199
		if (list(i) ~= list(i+1))   % we've found a boundary.
			actual_cutoffs = [actual_cutoffs (i+0.5)];
			% actual_cutoffs = [actual_cutoffs FindGaussianCrossover_2(G{list(i)},G{list(i+1)},i)];
			mostLikelyGaussians = [mostLikelyGaussians list(i)];
		end;
	end;
	mostLikelyGaussians = [mostLikelyGaussians list(200)];
	x_peak              = [];
	for i = 1:4; x_peak(i) = G{i}.b; end;

	% Construct curve
	fit_curve_1 = G{1}.a*exp(-0.5*((range-G{1}.b)./G{1}.c).^2);
	fit_curve_2 = G{2}.a*exp(-0.5*((range-G{2}.b)./G{2}.c).^2);
	fit_curve_3 = G{3}.a*exp(-0.5*((range-G{3}.b)./G{3}.c).^2);
	fit_curve_4 = G{4}.a*exp(-0.5*((range-G{4}.b)./G{4}.c).^2);
	plot(fit_curve_1,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_2,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_3,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_4,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	fit_curve_tot = fit_curve_1+fit_curve_2+fit_curve_3+fit_curve_4;
elseif (copyNum == 4)
	G                   = [];
	[G{1}.a,G{1}.b,G{1}.c, G{2}.a,G{2}.b,G{2}.c, G{3}.a,G{3}.b,G{3}.c, G{4}.a,G{4}.b,G{4}.c, G{5}.a,G{5}.b,G{5}.c] = ...
		fit_Gaussian_model_tetrasomy_2(workingDir,saveName, smoothed_Histogram,tetrasomy_peaks,sigma,ErrorType);
	[list]              = FindHighestGaussian_2(G);
	actual_cutoffs      = [];
	mostLikelyGaussians = [];
	for i = 1:199
		if (list(i) ~= list(i+1))   % we've found a boundary.
			actual_cutoffs = [actual_cutoffs (i+0.5)];
			% actual_cutoffs = [actual_cutoffs FindGaussianCrossover_2(G{list(i)},G{list(i+1)},i)];
			mostLikelyGaussians = [mostLikelyGaussians list(i)];
		end;
	end;
	mostLikelyGaussians = [mostLikelyGaussians list(200)];
	x_peak              = [];
	for i = 1:5; x_peak(i) = G{i}.b; end;

	% Construct curve
	fit_curve_1 = G{1}.a*exp(-0.5*((range-G{1}.b)./G{1}.c).^2);
	fit_curve_2 = G{2}.a*exp(-0.5*((range-G{2}.b)./G{2}.c).^2);
	fit_curve_3 = G{3}.a*exp(-0.5*((range-G{3}.b)./G{3}.c).^2);
	fit_curve_4 = G{4}.a*exp(-0.5*((range-G{4}.b)./G{4}.c).^2);
	fit_curve_5 = G{5}.a*exp(-0.5*((range-G{5}.b)./G{5}.c).^2);
	plot(fit_curve_1,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_2,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_3,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_4,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_5,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	fit_curve_tot = fit_curve_1+fit_curve_2+fit_curve_3+fit_curve_4+fit_curve_5;
elseif (copyNum == 5)
	G                   = [];
	[G{1}.a,G{1}.b,G{1}.c, G{2}.a,G{2}.b,G{2}.c, G{3}.a,G{3}.b,G{3}.c, G{4}.a,G{4}.b,G{4}.c, G{5}.a,G{5}.b,G{5}.c, G{6}.a,G{6}.b,G{6}.c] = ...
		fit_Gaussian_model_pentasomy_2(workingDir,saveName, smoothed_Histogram,pentasomy_peaks,sigma,ErrorType);
	[list]              = FindHighestGaussian_2(G);
	actual_cutoffs      = [];
	mostLikelyGaussians = [];
	for i = 1:199
		if (list(i) ~= list(i+1))   % we've found a boundary.
			actual_cutoffs = [actual_cutoffs (i+0.5)];
			% actual_cutoffs = [actual_cutoffs FindGaussianCrossover_2(G{list(i)},G{list(i+1)},i)];
			mostLikelyGaussians = [mostLikelyGaussians list(i)];
		end;
	end;
	mostLikelyGaussians = [mostLikelyGaussians list(200)];
	x_peak              = [];
	for i = 1:6; x_peak(i) = G{i}.b; end;

	% Construct curve
	fit_curve_1 = G{1}.a*exp(-0.5*((range-G{1}.b)./G{1}.c).^2);
	fit_curve_2 = G{2}.a*exp(-0.5*((range-G{2}.b)./G{2}.c).^2);
	fit_curve_3 = G{3}.a*exp(-0.5*((range-G{3}.b)./G{3}.c).^2);
	fit_curve_4 = G{4}.a*exp(-0.5*((range-G{4}.b)./G{4}.c).^2);
	fit_curve_5 = G{5}.a*exp(-0.5*((range-G{5}.b)./G{5}.c).^2);
	fit_curve_6 = G{6}.a*exp(-0.5*((range-G{6}.b)./G{6}.c).^2);
	plot(fit_curve_1,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_2,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_3,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_4,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_5,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_6,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	fit_curve_tot = fit_curve_1+fit_curve_2+fit_curve_3+fit_curve_4+fit_curve_5+fit_curve_6;
elseif (copyNum == 6)
	G                   = [];
	[G{1}.a,G{1}.b,G{1}.c, G{2}.a,G{2}.b,G{2}.c, G{3}.a,G{3}.b,G{3}.c, G{4}.a,G{4}.b,G{4}.c, G{5}.a,G{5}.b,G{5}.c, G{6}.a,G{6}.b,G{6}.c, G{7}.a,G{7}.b,G{7}.c] = ...
		fit_Gaussian_model_hexasomy_2(workingDir,saveName, smoothed_Histogram,hexasomy_peaks,sigma,ErrorType);
	[list]              = FindHighestGaussian_2(G);
	actual_cutoffs      = [];
	mostLikelyGaussians = [];
	for i = 1:199
		if (list(i) ~= list(i+1))   % we've found a boundary.
			actual_cutoffs = [actual_cutoffs (i+0.5)];
			% actual_cutoffs = [actual_cutoffs FindGaussianCrossover_2(G{list(i)},G{list(i+1)},i)];
			mostLikelyGaussians = [mostLikelyGaussians list(i)];
		end;
	end;
	mostLikelyGaussians = [mostLikelyGaussians list(200)];
	x_peak              = [];
	for i = 1:7; x_peak(i) = G{i}.b; end;

	% Construct curve
	fit_curve_1 = G{1}.a*exp(-0.5*((range-G{1}.b)./G{1}.c).^2);
	fit_curve_2 = G{2}.a*exp(-0.5*((range-G{2}.b)./G{2}.c).^2);
	fit_curve_3 = G{3}.a*exp(-0.5*((range-G{3}.b)./G{3}.c).^2);
	fit_curve_4 = G{4}.a*exp(-0.5*((range-G{4}.b)./G{4}.c).^2);
	fit_curve_5 = G{5}.a*exp(-0.5*((range-G{5}.b)./G{5}.c).^2);
	fit_curve_6 = G{6}.a*exp(-0.5*((range-G{6}.b)./G{6}.c).^2);
	fit_curve_7 = G{7}.a*exp(-0.5*((range-G{7}.b)./G{7}.c).^2);
	plot(fit_curve_1,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_2,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_3,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_4,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_5,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_6,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_7,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	fit_curve_tot = fit_curve_1+fit_curve_2+fit_curve_3+fit_curve_4+fit_curve_5+fit_curve_6+fit_curve_7;
elseif (copyNum == 7)
	G                   = [];
	[G{1}.a,G{1}.b,G{1}.c, G{2}.a,G{2}.b,G{2}.c, G{3}.a,G{3}.b,G{3}.c, G{4}.a,G{4}.b,G{4}.c, G{5}.a,G{5}.b,G{5}.c, G{6}.a,G{6}.b,G{6}.c, G{7}.a,G{7}.b,G{7}.c, G{8}.a,G{8}.b,G{8}.c] = ...
		fit_Gaussian_model_heptasomy_2(workingDir,saveName, smoothed_Histogram,heptasomy_peaks,sigma,ErrorType);	
	[list]              = FindHighestGaussian_2(G);
	%	fprintf(['||2 list            = ' num2str(list) '\n']);
	actual_cutoffs      = [];
	mostLikelyGaussians = [];
	for i = 1:199
		if (list(i) ~= list(i+1))   % we've found a boundary.
			actual_cutoffs = [actual_cutoffs (i+0.5)];
			% actual_cutoffs = [actual_cutoffs FindGaussianCrossover_2(G{list(i)},G{list(i+1)},i)];
			mostLikelyGaussians = [mostLikelyGaussians list(i)];
		end;
	end;
	mostLikelyGaussians = [mostLikelyGaussians list(200)];
	x_peak              = [];
	for i = 1:8; x_peak(i) = G{i}.b; end;

	% Construct curve
	fit_curve_1 = G{1}.a*exp(-0.5*((range-G{1}.b)./G{1}.c).^2);
	fit_curve_2 = G{2}.a*exp(-0.5*((range-G{2}.b)./G{2}.c).^2);
	fit_curve_3 = G{3}.a*exp(-0.5*((range-G{3}.b)./G{3}.c).^2);
	fit_curve_4 = G{4}.a*exp(-0.5*((range-G{4}.b)./G{4}.c).^2);
	fit_curve_5 = G{5}.a*exp(-0.5*((range-G{5}.b)./G{5}.c).^2);
	fit_curve_6 = G{6}.a*exp(-0.5*((range-G{6}.b)./G{6}.c).^2);
	fit_curve_7 = G{7}.a*exp(-0.5*((range-G{7}.b)./G{7}.c).^2);
	fit_curve_8 = G{8}.a*exp(-0.5*((range-G{8}.b)./G{8}.c).^2);
	plot(fit_curve_1,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_2,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_3,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_4,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_5,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_6,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_7,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_8,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	fit_curve_tot = fit_curve_1+fit_curve_2+fit_curve_3+fit_curve_4+fit_curve_5+fit_curve_6+fit_curve_7+fit_curve_8;
elseif (copyNum == 8)
	G                  = [];
	[G{1}.a,G{1}.b,G{1}.c, G{2}.a,G{2}.b,G{2}.c, G{3}.a,G{3}.b,G{3}.c, G{4}.a,G{4}.b,G{4}.c, G{5}.a,G{5}.b,G{5}.c, G{6}.a,G{6}.b,G{6}.c, G{7}.a,G{7}.b,G{7}.c, G{8}.a,G{8}.b,G{8}.c, G{9}.a,G{9}.b,G{9}.c] = ...
		fit_Gaussian_model_octasomy_2(workingDir,saveName, smoothed_Histogram,octasomy_peaks,sigma,ErrorType);
	[list]              = FindHighestGaussian_2(G);
	actual_cutoffs      = [];
	mostLikelyGaussians = [];
	for i = 1:199
		if (list(i) ~= list(i+1))   % we've found a boundary.
			actual_cutoffs = [actual_cutoffs (i+0.5)];
			% actual_cutoffs = [actual_cutoffs FindGaussianCrossover_2(G{list(i)},G{list(i+1)},i)];
			mostLikelyGaussians = [mostLikelyGaussians list(i)];
		end;
	end;
	mostLikelyGaussians = [mostLikelyGaussians list(200)];
	x_peak              = [];
	for i = 1:9; x_peak(i) = G{i}.b; end;

	% Construct curve
	Fit_curve_1 = G{1}.a*exp(-0.5*((range-G{1}.b)./G{1}.c).^2);
	fit_curve_2 = G{2}.a*exp(-0.5*((range-G{2}.b)./G{2}.c).^2);
	fit_curve_3 = G{3}.a*exp(-0.5*((range-G{3}.b)./G{3}.c).^2);
	fit_curve_4 = G{4}.a*exp(-0.5*((range-G{4}.b)./G{4}.c).^2);
	fit_curve_5 = G{5}.a*exp(-0.5*((range-G{5}.b)./G{5}.c).^2);
	fit_curve_6 = G{6}.a*exp(-0.5*((range-G{6}.b)./G{6}.c).^2);
	fit_curve_7 = G{7}.a*exp(-0.5*((range-G{7}.b)./G{7}.c).^2);
	fit_curve_8 = G{8}.a*exp(-0.5*((range-G{8}.b)./G{8}.c).^2);
	fit_curve_9 = G{9}.a*exp(-0.5*((range-G{9}.b)./G{9}.c).^2);
	plot(fit_curve_1,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_2,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_3,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_4,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_5,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_6,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_7,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_8,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_9,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	fit_curve_tot = fit_curve_1+fit_curve_2+fit_curve_3+fit_curve_4+fit_curve_5+fit_curve_6+fit_curve_7+fit_curve_8+fit_curve_9;
else % if (copyNum == 9+)
	G                   = [];
	[G{1}.a,G{1}.b,G{1}.c, G{2}.a,G{2}.b,G{2}.c, G{3}.a,G{3}.b,G{3}.c, G{4}.a,G{4}.b,G{4}.c, G{5}.a,G{5}.b,G{5}.c, G{6}.a,G{6}.b,G{6}.c, G{7}.a,G{7}.b,G{7}.c, G{8}.a,G{8}.b,G{8}.c, G{9}.a,G{9}.b,G{9}.c, G{10}.a,G{10}.b,G{10}.c] = ...
		fit_Gaussian_model_nonasomy_2(workingDir,saveName, smoothed_Histogram,nonasomy_peaks,sigma,ErrorType);
	[list]              = FindHighestGaussian_2(G);
	actual_cutoffs      = [];
	mostLikelyGaussians = [];
	for i = 1:199
		if (list(i) ~= list(i+1))   % we've found a boundary.
			actual_cutoffs = [actual_cutoffs (i+0.5)];
			% actual_cutoffs = [actual_cutoffs FindGaussianCrossover_2(G{list(i)},G{list(i+1)},i)];
			mostLikelyGaussians = [mostLikelyGaussians list(i)];
		end;
	end;
	mostLikelyGaussians = [mostLikelyGaussians list(200)];
	x_peak              = [];
	for i = 1:9; x_peak(i) = G{i}.b; end;

	% Construct curve
	fit_curve_1  = G{ 1}.a*exp(-0.5*((range-G{ 1}.b)./G{ 1}.c).^2);
	fit_curve_2  = G{ 2}.a*exp(-0.5*((range-G{ 2}.b)./G{ 2}.c).^2);
	fit_curve_3  = G{ 3}.a*exp(-0.5*((range-G{ 3}.b)./G{ 3}.c).^2);
	fit_curve_4  = G{ 4}.a*exp(-0.5*((range-G{ 4}.b)./G{ 4}.c).^2);
	fit_curve_5  = G{ 5}.a*exp(-0.5*((range-G{ 5}.b)./G{ 5}.c).^2);
	fit_curve_6  = G{ 6}.a*exp(-0.5*((range-G{ 6}.b)./G{ 6}.c).^2);
	fit_curve_7  = G{ 7}.a*exp(-0.5*((range-G{ 7}.b)./G{ 7}.c).^2);
	fit_curve_8  = G{ 8}.a*exp(-0.5*((range-G{ 8}.b)./G{ 8}.c).^2);
	fit_curve_9  = G{ 9}.a*exp(-0.5*((range-G{ 9}.b)./G{ 9}.c).^2);
	fit_curve_10 = G{10}.a*exp(-0.5*((range-G{10}.b)./G{10}.c).^2);
	plot(fit_curve_1,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_2,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_3,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_4,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_5,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_6,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_7,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_8,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_9,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	plot(fit_curve_10,'color',[1.00 0.50 0.50],'linestyle','-','linewidth',1);
	fit_curve_tot = fit_curve_1+fit_curve_2+fit_curve_3+fit_curve_4+fit_curve_5+fit_curve_6+fit_curve_7+fit_curve_8+fit_curve_9+fit_curve_10;
end;

%	fprintf('List = ');
%	if (length(list) > 0)
%		for i = 1:length(list)
%			fprintf(num2str(list(i)));
%		end;
%	end;
%	fprintf('\n');
%	fprintf('G    = ');
%	if (length(G) > 0)
%		for i = 1:length(G);
%			fprintf([num2str(G{i}.a) ':' num2str(G{i}.b) ':' num2str(G{i}.c) ' ']);
%		end;
%	end;
%	fprintf('\n');
%	fprintf(['^^^     Glist = ' num2str(list) '\n']);

if (MakeFigure == true)
	plot(smoothed_Histogram,'color',[0.50 0.50 1.00],'linestyle','-','linewidth',1);
	plot(fit_curve_tot,     'color',[1.00 0.00 0.00],'linestyle','-','linewidth',2);
	title(['allelicRatios.chr_' num2str(chromosome) '.segment_' num2str(segment)],'HorizontalAlign','center','VerticalAlign','middle');
	hold off;
	xlim([1,200]);
	% save then delete figures.
	saveas(fig, [workingDir 'hist_chr-' num2str(chromosome) '_seg-' num2str(segment) '.png'], 'png');
	delete(fig);
end;

end

