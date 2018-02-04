function ys=mylowess(xy,xs,span)

%MYLOWESS Lowess smoothing, preserving x values
%   YS=MYLOWESS(XY,XS) returns the smoothed version of the x/y data in the
%   two-column matrix XY, but evaluates the smooth at XS and returns the
%   smoothed values in YS.  Any values outside the range of XY are taken to
%   be equal to the closest values.

if nargin<3 || isempty(span)
    span = .3;
end

% Sort and get smoothed version of xy data
arrayDim  = size(xy);
if (arrayDim(1) > arrayDim(2))
	xy    = sortrows(xy);
else
	xy    = sortrows(xy');
end;
x1        = xy(:,1);
y1        = xy(:,2);

ys1       = smooth(x1,y1,span,'rlowess');
old_ys1 = ys1;

% Remove repeats so we can interpolate
t         = diff(x1)==0;
x1(t)     = [];
ys1(t)    = [];

% Interpolate to evaluate this at the xs values
if length(ys1) > 1
    ys     = interp1(x1,ys1,xs,'linear',NaN);

    % Some of the original points may have x values outside the range of the
    % resampled data.  Those are now NaN because we could not interpolate them.
    % Replace NaN by the closest smoothed value.  This amounts to extending the
    % smooth curve using a horizontal line.
    if any(isnan(ys))
        ys(xs<x1(1))   = ys1(1);
        ys(xs>x1(end)) = ys1(end);
    end
else
    ys = old_ys1;
end;