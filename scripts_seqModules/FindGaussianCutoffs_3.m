function [x_peak,actual_cutoffs,mostLikelyGaussians] = FindGaussianCutoffs_3(segment_smoothedHistogram, segment_copyNum);

monosomy_peaks  = [0, 1]*199+1;
disomy_peaks    = [0, 1/2, 1]*199+1;
trisomy_peaks   = [0, 1/3, 2/3, 1]*199+1;
tetrasomy_peaks = [0, 1/4, 1/2, 3/4, 1]*199+1;
pentasomy_peaks = [0, 1/5, 2/5, 3/5, 4/5, 1]*199+1;
hexasomy_peaks  = [0, 1/6, 1/3, 1/2, 2/3, 5/6, 1]*199+1;
heptasomy_peaks = [0, 1/7, 2/7, 3/7, 4/7, 5/7, 6/7, 1]*199+1;
octasomy_peaks  = [0, 1/8, 1/4, 3/8, 1/2, 5/8, 3/4, 7/8, 1]*199+1;

%% Calculation of Gaussians against per chromosome data.
% Fits Gaussians to real data per chromomsome, then determines equal probability cutoffs between them.
sigma = 5;

%% FindGaussianCutoffs Finds cutoffs as intersections of Gaussians, fit to the data at each peak location.
ErrorType      = 'linear';

if (segment_copyNum == 0)
	G                   = [];
	list                = [];
	x_peak              = [];
	actual_cutoffs      = [];
	mostLikelyGaussians = [];
elseif (segment_copyNum == 1)
	G                   = [];
	[G{1}.a,G{1}.b,G{1}.c, G{2}.a,G{2}.b,G{2}.c, skew_factor] = ...
	    fit_Gaussian_model_monosomy_2(segment_smoothedHistogram,monosomy_peaks,sigma,ErrorType);
	for i = 1:(segment_copyNum+1); G{i}.d = skew_factor; end;
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
elseif (segment_copyNum == 2)
	G                   = [];
	[G{1}.a,G{1}.b,G{1}.c, G{2}.a,G{2}.b,G{2}.c, G{3}.a,G{3}.b,G{3}.c, skew_factor] = ...
	    fit_Gaussian_model_disomy_2(segment_smoothedHistogram,disomy_peaks,sigma,ErrorType);
	for i = 1:(segment_copyNum+1); G{i}.d = skew_factor; end;
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
elseif (segment_copyNum == 3)
	G                   = [];
	[G{1}.a,G{1}.b,G{1}.c, G{2}.a,G{2}.b,G{2}.c, G{3}.a,G{3}.b,G{3}.c, G{4}.a,G{4}.b,G{4}.c, skew_factor] = ...
	    fit_Gaussian_model_trisomy_2(segment_smoothedHistogram,trisomy_peaks,sigma,ErrorType);
	for i = 1:(segment_copyNum+1); G{i}.d = skew_factor; end;
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
elseif (segment_copyNum == 4)
	G                   = [];
	[G{1}.a,G{1}.b,G{1}.c, G{2}.a,G{2}.b,G{2}.c, G{3}.a,G{3}.b,G{3}.c, G{4}.a,G{4}.b,G{4}.c, G{5}.a,G{5}.b,G{5}.c, skew_factor] = ...
	    fit_Gaussian_model_tetrasomy_2(segment_smoothedHistogram,tetrasomy_peaks,sigma,ErrorType);
	for i = 1:(segment_copyNum+1); G{i}.d = skew_factor; end;
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
elseif (segment_copyNum == 5)
	G                   = [];
	[G{1}.a,G{1}.b,G{1}.c, G{2}.a,G{2}.b,G{2}.c, G{3}.a,G{3}.b,G{3}.c, G{4}.a,G{4}.b,G{4}.c, G{5}.a,G{5}.b,G{5}.c, G{6}.a,G{6}.b,G{6}.c, skew_factor] = ...
	    fit_Gaussian_model_pentasomy_2(segment_smoothedHistogram,pentasomy_peaks,sigma,ErrorType);
	for i = 1:(segment_copyNum+1); G{i}.d = skew_factor; end;
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
elseif (segment_copyNum == 6)
	G                   = [];
	[G{1}.a,G{1}.b,G{1}.c, G{2}.a,G{2}.b,G{2}.c, G{3}.a,G{3}.b,G{3}.c, G{4}.a,G{4}.b,G{4}.c, G{5}.a,G{5}.b,G{5}.c, G{6}.a,G{6}.b,G{6}.c, G{7}.a,G{7}.b,G{7}.c, skew_factor] = ...
	    fit_Gaussian_model_hexasomy_2(segment_smoothedHistogram,hexasomy_peaks,sigma,ErrorType);
	for i = 1:(segment_copyNum+1); G{i}.d = skew_factor; end;
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
elseif (segment_copyNum == 7)
	G                   = [];
	[G{1}.a,G{1}.b,G{1}.c, G{2}.a,G{2}.b,G{2}.c, G{3}.a,G{3}.b,G{3}.c, G{4}.a,G{4}.b,G{4}.c, G{5}.a,G{5}.b,G{5}.c, G{6}.a,G{6}.b,G{6}.c, G{7}.a,G{7}.b,G{7}.c, G{8}.a,G{8}.b,G{8}.c, skew_factor] = ...
		fit_Gaussian_model_heptasomy_2(segment_smoothedHistogram,heptasomy_peaks,sigma,ErrorType);	
	for i = 1:(segment_copyNum+1); G{i}.d = skew_factor; end;
	[list]              = FindHighestGaussian_2(G);
	fprintf(['||2 list            = ' num2str(list) '\n']);
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
else % if (segment_copyNum == 8+)
	G                   = [];
	[G{1}.a,G{1}.b,G{1}.c, G{2}.a,G{2}.b,G{2}.c, G{3}.a,G{3}.b,G{3}.c, G{4}.a,G{4}.b,G{4}.c, G{5}.a,G{5}.b,G{5}.c, G{6}.a,G{6}.b,G{6}.c, G{7}.a,G{7}.b,G{7}.c, G{8}.a,G{8}.b,G{8}.c, G{9}.a,G{9}.b,G{9}.c, skew_factor] = ...
		fit_Gaussian_model_octasomy_2(segment_smoothedHistogram,octasomy_peaks,sigma,ErrorType);
	for i = 1:(length(G));   G{i}.d = skew_factor;   end;
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
end;

fprintf('List = ');
if (length(list) > 0)
	for i = 1:length(list)
		fprintf(num2str(list(i)));
	end;
end;
fprintf('\n');
fprintf('G    = ');
if (length(G) > 0)
	for i = 1:length(G);
		fprintf([num2str(G{i}.a) ':' num2str(G{i}.b) ':' num2str(G{i}.c) ':' num2str(G{i}.d) ' ']);
	end;
end;
fprintf('\n');

% fprintf(['^^^     Glist = ' num2str(list) '\n']);

end

