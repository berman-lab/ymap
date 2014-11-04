function [ G1_pos,G2_pos, G1G2_crossover ] = find_2_Gaussian_fits_from_raw_data( data1,data2 )
% Calculates a Gaussian-smoothed weighted distribution
%   data1             : raw data.
%   data2             : weight of raw data.
resolution    = 100;
totalDistance = 0;
count         = 0;

% Calculate the average distance between values in 'data1'.
for i = 1:(length(data1)-1)
	for j = (i+1):length(data1)
		totalDistance = totalDistance + abs(data1(j)-data1(i));
		count         = count + 1;
	end;
end;
averageGap     = totalDistance/count;
GaussianWidth  = averageGap*4;

% Construct 'dataHistogram' vector as histogram of weighted data points.
dataHistogram  = zeros(1,round(max(data1)*resolution*1.2)+1);
for i = 1:length(data1)
	dataHistogram(round(data1(i)*resolution)+1) = dataHistogram(round(data1(i)*resolution)+1)+data2(i);
end;

% function s = smooth_gaussian(data,sigma,size)
% %    data  : input vector with raw data.
% %    sigma : standard deviation of the gaussian distribution used in the smoothing.
% %    size  : size of vector over which smoothing function is applied.   (2-3 sigmas is usually good.)

% smooth display vector 'dataHistogram'.
smoothedDataHistogram  = smooth_gaussian(dataHistogram, GaussianWidth, GaussianWidth*10);

newX = ((1:length(smoothedDataHistogram)) - 1)/resolution;
newY = smoothedDataHistogram;

G1 = [];
G2 = [];
[G1.a, G1.b, G1.c, G2.a, G2.b, G2.c] = fit2GaussianModels(newX,newY,'cubic');
% a = height.
% b = location.
% c = width.
G1_pos         = G1.b/resolution;
G2_pos         = G2.b/resolution;
G1G2_crossover = FindGaussianCrossover(G1,G2)/resolution;

end
