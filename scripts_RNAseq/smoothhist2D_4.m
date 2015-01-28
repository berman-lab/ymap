function [imageX, imageY, imageC] = smoothhist2D_4(x,y,lambda,nbins,rangeMax,outliercutoff)
%
%  This version allows user control of the colormap in use.   Only standard matlab colormaps are valid.
%
if (nargin<7) || isempty(colormap_)
	colormap_ = 'hot';
end;

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

if (nargin < 6) || (isempty(outliercutoff))
	outliercutoff = 0;%.02;
end;

if (length(x) > 499999)
	x2 = x(1:500000);
	y2 = y(1:500000);
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
H = accumarray(bin,1,nbins([2 1]))./n;

% Eiler's 1D smooth, twice
G = smooth1D(H,lambda)*10000;
F = smooth1D(G',lambda)';
% % An alternative, using filter2.  However, lambda means totally different things in this case: for smooth1D, it is a smoothness penalty parameter,  while for filter2D, it is a smoothing window halfwidth
% F = filter2D(H,lambda);

%	relF = F./max(F(:));
relF = F/max(max(F));
if outliercutoff > 0
	outliers = (relF(nbins(2)*(bin(:,2)-1)+bin(:,1)) < outliercutoff);
end

%	edge_gap = 5; % distance from edges to ignore max color value.
%	dd       = max(relF(edge_gap:(nbins-1-edge_gap),edge_gap:(nbins-1-edge_gap))');
%	maxCol   = max(dd);
trimmed_relF = relF;
maxCol   = max(max(trimmed_relF));

nc = 256;
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

%fighandle = image(ctrs1,ctrs2,floor(nc.*log(relF./maxCol)));
imageX = ctrs1;
imageY = ctrs2;
%	imageC = floor(nc*log(relF/maxCol+1));
imageC = floor(370*log(relF+1));

%	fighandle = image(imageX, imageY, imageC);
%	fighandle = image(ctrs1,ctrs2,floor(nc.*log(relF./maxCol+1.1)));
%	axis square;
%	set(gca,'YDir','normal');
%	hold on;
%	% plot the outliers
%	if (outliercutoff > 0)
%		plot(X(outliers,1),X(outliers,2),'.','MarkerEdgeColor',[1 1 1],'MarkerSize',1);
%	end
%	% % plot a subsample of the data
%	% Xsample = X(randsample(n,n/10),:);
%	% plot(Xsample(:,1),Xsample(:,2),'bo');
%	hold off;

%-----------------------------------------------------------------------------
function Z = smooth1D(Y,lambda)
[m,n] = size(Y);
E = eye(m);
D1 = diff(E,1);
D2 = diff(D1,1);
P = lambda.^2 .* D2'*D2 + 2.*lambda .* D1'*D1;
Z = (E + P) \ Y;
% This is a better solution, but takes a bit longer for n and m large
% opts.RECT = true;
% D1 = [diff(E,1); zeros(1,n)];
% D2 = [diff(D1,1); zeros(1,n)];
% Z = linsolve([E; 2.*sqrt(lambda).*D1; lambda.*D2],[Y; zeros(2*m,n)],opts);
%-----------------------------------------------------------------------------
function Z = filter2D(Y,bw)
z = -1:(1/bw):1;
k = .75 * (1 - z.^2); % epanechnikov-like weights
k = k ./ sum(k);
Z = filter2(k'*k,Y);
