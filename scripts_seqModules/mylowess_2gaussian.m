function [ys_1,ys_2,ys_3] = mylowess_2gaussian(xy,xs,span, yLimit)
%
%
% Not functional yet.
%
%

%MYLOWESS Lowess smoothing, preserving x values
%   YS=MYLOWESS(XY,XS) returns the smoothed version of the x/y data in the
%   two-column matrix XY, but evaluates the smooth at XS and returns the
%   smoothed values in YS.  Any values outside the range of XY are taken to
%   be equal to the closest values.
%
%   Performs a 2-part mixed Gaussian fit to the localized data.

if nargin<3 || isempty(span)
    span = 0.3;
end

% Sort and filter xy data to usable range.
xy             = sortrows(xy);
x              = xy(:,1);
y              = xy(:,2);
x(y >= yLimit) = [];
y(y >= yLimit) = [];


%% ===========================================================================================
% Perform dual-Gaussian fitting on locally-weighted data using variation of LOESS algorithm.
%---------------------------------------------------------------------------------------------
iter     = 5;
method   = 'loess';  % quadratic fit.
n        = length(y);
span     = ceil(span*n);
span     = floor(span);
span     = min(span,n);
if (span == 1)
	return;
end
c1       = y;
c2       = y;
c3       = y;
useLoess = false;
diffx    = diff(x); % Calculates differences between adjacent values of x.

% Turn off warnings when called from command line (already off if called from cftool).
ws                       = warning('off', 'all'); % save warning state
[lastwarnmsg,lastwarnid] = lastwarn;  % save last warning

ynan      = isnan(y);
anyNans   = any(ynan(:));
seps      = sqrt(eps);
theDiffs  = [1; diffx; 1];

% Compute the smooth for non-uniform x
for i = 1:n
	% if x(i) and x(i-1) are equal we just use the old value.
	if (theDiffs(i) == 0)
		c1(i) = c1(i-1);
		c2(i) = c2(i-1);
		c3(i) = c3(i-1);
		continue;
	end

	% Determine the data covered by the span of interest.
	left  = max(1,i-span+1);
	right = min(n,i+span-1);

	% now see if we have any adjacent_equal values that we need to take into account
	while (left > 0) && (theDiffs(left) == 0)
		left = left-1;
	end
	while (right <= n) && (theDiffs(right+1) == 0)
		right = right+1;
	end

	% Load current location value...  center around current point to improve conditioning.
	mx = x(i);

	% Look at the span interval around x(i).
	d           = abs(x(left:right)-mx);
	[dsort,idx] = sort(d);

	% Add back left value => generates sorted indexes for span of interest.
	idx         = idx + left - 1;

	% Remove indexes for any 'NaN' values.
	if anyNans
		idx = idx(dsort <= dsort(span) & ~ynan(idx));
	else
		idx = idx(dsort <= dsort(span));
	end

	% Checks if index list retains any values after filterin out 'NaN's.
	if isempty(idx)
		c(i) = NaN;
		continue
	end

	x1    = x(idx) - mx;        % X values of span, centered around coordinate.
	y1    = y(idx);             % Y values of span.
	dsort = d(idx - left + 1);  %
	dmax  = dsort(end);
	if (dmax == 0)
		dmax = 1;
	end

	% Define weights for smoothing of local values from scatterplot.
	weight = (1 - (dsort/dmax).^3).^1.5; % tri-cubic weight
	if all(weight < seps)
		weight(:) = 1;    % if all weights are 0, just skip weighting
	end

	% Scale localized Y values by weights.
	v  = [ones(size(x1)) x1];
	v  = weight(:,ones(1,size(v,2))) .* v;
%	y1 = weight .* y1;

	% Fit 2-Gaussian distribution to histogram of weighted Y data.
	[ G1_pos,G2_pos, G1G2_crossover ] = find_2_Gaussian_fits_from_raw_data( y1,weight );

	c1(i) = G1_pos;   % b(1);
	c2(i) = G2_pos;
	c3(i) = G1G2_crossover;
end

%% ===========================================================================================
% 
%---------------------------------------------------------------------------------------------

ys1     = c1;
ys2     = c2;
ys3     = c3;
% Remove repeats so we can interpolate
t1      = diff(x)==0;
x(t1)   = [];
ys1(t1) = [];
ys2(t1) = [];
ys3(t1) = [];

% Interpolate to evaluate this at the xs values
ys_1    = interp1(x,ys1,xs,'linear',NaN);
ys_2    = interp1(x,ys2,xs,'linear',NaN);
ys_3    = interp1(x,ys3,xs,'linear',NaN);

% Some of the original points may have x values outside the range of the
% resampled data.  Those are now NaN because we could not interpolate them.
% Replace NaN by the closest smoothed value.  This amounts to extending the
% smooth curve using a horizontal line.
if any(isnan(ys_1))
	ys_1(xs<x(1))   = ys1(1);
	ys_1(xs>x(end)) = ys1(end);
end
if any(isnan(ys_2))
	ys_2(xs<x(1))   = ys2(1);
	ys_2(xs>x(end)) = ys2(end);
end
if any(isnan(ys_3))
	ys_3(xs<x(1))   = ys3(1);
	ys_3(xs>x(end)) = ys3(end);
end
