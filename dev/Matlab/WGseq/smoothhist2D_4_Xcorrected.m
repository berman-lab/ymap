function [imageX, imageY, imageC, imageD] = smoothhist2D_4_Xcorrected(x,y,lambda,nbins,rangeMax, XcorrectionVector, mean_val, scaler)
%
%  This version allows user control of the colormap in use.   Only standard matlab colormaps are valid.
%  This version also allows the user to apply a x-correction to the datapoints before smoothing.
%
colormap_ = 'hot';

% SMOOTHHIST2D Plot a smoothed histogram of bivariate data.
%	SMOOTHHIST2D(x,y,LAMBDA,NBINS) plots a smoothed histogram of the bivariate
%	data in the vectores x and y.  Rows correspond to observations.  LAMBDA
%	is a positive scalar smoothing parameter; higher values lead to more
%	smoothing, values close to zero lead to a plot that is essentially just
%	the raw data.  NBINS is a two-element vector that determines the number
%	of histogram bins in the horizontal and vertical directions.
%
%	SMOOTHHIST2D(X,LAMBDA,NBINS,CUTOFF) plots outliers in the data as points
%	overlaid on the smoothed histogram.  Outliers are defined as points in
%	regions where the smoothed density is less than (100*CUTOFF)% of the
%	maximum density.
%
%	SMOOTHHIST2D(X,LAMBDA,NBINS,[],'surf') plots a smoothed histogram as a
%	surface plot.  SMOOTHHIST2D ignores the CUTOFF input in this case, and
%	the surface plot does not include outliers.
%
%	SMOOTHHIST2D(X,LAMBDA,NBINS,CUTOFF,'image') plots the histogram as an
%	image plot, the default.
%
%	Example:
%		X = [mvnrnd([0 5], [3 0; 0 3], 2000);
%		     mvnrnd([0 8], [1 0; 0 5], 2000);
%		     mvnrnd([3 5], [5 0; 0 1], 2000)];
%		smoothhist2D(X,5,[100, 100],.05);
%		smoothhist2D(X,5,[100, 100],[],'surf');
%
%	Reference:
%		Eilers, P.H.C. and Goeman, J.J (2004) "Enhancing scaterplots with
%		smoothed densities", Bioinformatics 20(5):623-628.
%	Copyright 2009 The MathWorks, Inc.
%	Revision: 1.0  Date: 2006/12/12
%
%	Requires MATLABÆ R14.

outliercutoff = 0;
if (nargin < 6) || (isempty(outliercutoff))
	outliercutoff = 0;   %.02;
end;

% data limited to 500,000 points in original version.
data_limit = 1000000;
if (length(x) > data_limit)
	% randomize the order of data in the input vectors.
	randIndex = randperm(length(x));
	x_temp = x(randIndex);
	y_temp = y(randIndex);
	% Limit the randomized data to the data_limit.
	x2 = x_temp(1:data_limit);
	y2 = y_temp(1:data_limit);
else
	x2 = x;
	y2 = y;
end;


X = [x2;y2]';
%X = [x;y]';


minx = min(X,[],1);
maxx = max(X,[],1);
if (rangeMax(1) > 0)
	if (maxx(1) > rangeMax(1));
		maxx(1) = rangeMax(1);
	end;
end;
if (rangeMax(2) > 0)
	if (maxx(2) > rangeMax(2));
		maxx(2) = rangeMax(2);
	end;
end;
edges1 = linspace(minx(1), maxx(1), nbins(1)+1);
ctrs1  = edges1(1:end-1) + .5*diff(edges1);
edges1 = [-Inf edges1(2:end-1) Inf];
edges2 = linspace(minx(2), maxx(2), nbins(2)+1);
ctrs2  = edges2(1:end-1) + .5*diff(edges2);
edges2 = [-Inf edges2(2:end-1) Inf];

[n,p] = size(X);
bin = zeros(n,2);
% Reverse the columns of H to put the first column of X along the
% horizontal axis, the second along the vertical.
[dum,bin(:,2)] = histc(X(:,1),edges1);
[dum,bin(:,1)] = histc(X(:,2),edges2);

%% Accumulate points into 2D array.
H = accumarray(bin,1,nbins([2 1]))./n;


%% Apply the X correction factor: this intended for use to normalize SNP data by CNV estimate.
%    H                 : Accumulated data.
%    XcorrectionVector : Input correction per X bin. The correction factor being squared was determined empirically.
H_original = H;
for i = 1:size(H,2)
	H(:,i) = H(:,i)*(XcorrectionVector(i).^1);  %2?
end;


%% Crop 2D histogram to three times the median.
temp_H              = H;
temp_H(temp_H == 0) = [];
max_val             = mean_val*100;
H(H > max_val)      = max_val;


%% Trouble-shooting outputs.
%	testCorrection_length = length(XcorrectionVector)
%	testDataset_length    = length(H)
%	testCorrection        = XcorrectionVector
%	testDataset           = H
%	H_median              = median(H(:))
%	H_mode                = mode(H(:))
%	H_mean                = mean(H(:))


%% Perform smoothing of 2D histogram: Eiler's 1D smooth, twice
G = smooth1D(H,lambda)*10000;
F = smooth1D(G',lambda)';

%	relF = F./max(F(:));
relF = F/max(max(F));
if outliercutoff > 0
	outliers = (relF(nbins(2)*(bin(:,2)-1)+bin(:,1)) < outliercutoff);
end

%% Scale data by input variable "scaler", which is used to correct for per-chromosome differences in data levels.
relF = relF*scaler;

nc = 256;

imageX = ctrs1;
imageY = ctrs2;
imageC = floor(nc*2*log(relF+1));
imageD = H_original;

switch lower(colormap_)
	case 'jet'
		colormap(jet(nc));
	case 'hsv'
		colormap(hsv(nc));
	case 'hot'
		colormap(hot(nc));
	case 'cool'
		colormap(cool(nc));
	case 'spring'
		colormap(spring(nc));
	case 'summer'
		colormap(summer(nc));
	case 'autumn'   
		colormap(autumn(nc));
	case 'winter'
		colormap(winter(nc));
	case 'gray'
		colormap(gray(nc));
	case 'bone'
		colormap(bone(nc));
	case 'copper'
		colormap(copper(nc));
	case 'pink'
		colormap(pink(nc));
	case 'lines'
		colormap(lines(nc));
end;


%-----------------------------------------------------------------------------
function Z = smooth1D(Y,lambda)
[m,n] = size(Y);
E = eye(m);
D1 = diff(E,1);
D2 = diff(D1,1);
P = lambda.^2 .* D2'*D2 + 2.*lambda .* D1'*D1;
Z = (E + P) \ Y;
%-----------------------------------------------------------------------------
function Z = filter2D(Y,bw)
z = -1:(1/bw):1;
k = .75 * (1 - z.^2); % epanechnikov-like weights
k = k ./ sum(k);
Z = filter2(k'*k,Y);
