function s = smooth_gaussian(data,sigma,size)
%    data  : input vector with raw data.
%    sigma : standard deviation of the gaussian distribution used in the smoothing.
%    size  : size of vector over which smoothing function is applied.   (2-3 sigmas is usually good.)

% Gaussian smoothing.
halfsize = round(size/2);
a        = 1/(sqrt(2*pi)*sigma);
b        = 1/(2*sigma^2);
w        = a*exp(-b*(-halfsize:1:halfsize).^2);

% filters data and shifts smoothing left to align with data.
% s1 = left end of smoothed data.
s1   = circshift(filter(w,1,data),[1 -halfsize]);

% filters right end of data in reverse orientation.
% s2  = right end of smoothed data.
if (length(data) > size)
	data2 = fliplr(data(length(data)-size:end));
else
	data2 = fliplr(data);
end;
s2    = fliplr(filter(w,1,data2));

% crops away the unreliable sections of the two smooth operations
s1(length(data)-halfsize:end) = [];
s2(halfsize+2:end)            = [];

% combines the two smooths to generate a consistent smooth for the user.
s = [s1 s2];
