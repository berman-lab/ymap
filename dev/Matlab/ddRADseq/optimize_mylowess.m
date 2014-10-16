function [newX, newY] = optimize_mylowess(rawData_X,rawData_Y, numFits, maxX)

numDat        = length(rawData_X);
spans         = linspace(0.01,0.99, numFits);
sse           = zeros(size(spans));
cp            = cvpartition(numDat,'k',10);
X             = rawData_X';
Y             = rawData_Y';
final_length  = 250;     % length of resulting fitted lowess vector.
if (maxX == 0)
	maxX = max(X);
	minX = min(X);
else
	minX = 0;
end;

LOWESS_method = 3;

if (LOWESS_method == 1)
	% Attempts LOWESS fitting with [numFits] smoothing values evenly spaced from 0.01 to 0.99.
	sse  = zeros(size(spans));
	fprintf(['\t\tSquared-error minimization:\n']);
	for j = 1:numFits
		fprintf(['\t\t\tFit type1 #' num2str(j) '/' num2str(numFits) '.\t[']);
		f      = @(train,test) norm(test(:,2) - mylowess(train,test(:,1),spans(j)))^2;
		fprintf(':');
		sse(j) = sum(crossval(f,[X,Y],'partition',cp));
		fprintf(['](error = ' num2str(sse(j)) ')\n']);
	end;
elseif (LOWESS_method == 2)
	% Attempts LOWESS fitting with [numFits] smoothing values evenly spaced from 0.01 to 0.99.
	sse  = zeros(size(spans));
	fprintf(['\t\tSquared-error minimization:\n']);
	for j = 1:numFits
		fprintf(['\t\t\tFit type2 #' num2str(j) '/' num2str(numFits) '.\t[']);
		% perform LOWESS fitting with current span.
		arrayDim = size(X);
		if (arrayDim(1) > arrayDim(2))
			Ys  = mylowess([X, Y],X,spans(j));
		else
			Ys  = mylowess([X', Y'],X,spans(j));
		end;
		fprintf(':');
		% determine error between LOWESS fit and dataset.
		LSerrorSum  = 0;
		LSerrorSum2 = 0;
		for i = 1:length(Y)
			LSerrorSum  = LSerrorSum  + (Y(i)-Ys(i))^2;
		end;
		sse(j)  = LSerrorSum;
		fprintf(['](error = ' num2str(sse(j)) ')\n']);
	end;
elseif (LOWESS_method == 3)
	% Attempts LOWESS fitting with [numFits] smoothing values evenly spaced from 0.01 to 0.99, using 10-fold cross-validation of randomly partitioned data.
	fprintf(['\t\t10-fold cross validation, with squared-error minimization:\n']);

	fitCurves = cell(1,numFits);
	for j = 1:numFits
		fprintf(['\t\t\tFit type3 #' num2str(j) '/' num2str(numFits) '.\t[']);

		% Randomly sort the input data into 10 partitions.
		randIndex       = randperm(length(rawData_X));      % random order of length the length of the data.
		partitionData_X = cell(1,10);
		partitionData_Y = cell(1,10);
		partitionIndex  = 1;

		for i = 1:length(rawData_X);
			partitionData_X{partitionIndex} = [partitionData_X{partitionIndex} rawData_X(i)];
			partitionData_Y{partitionIndex} = [partitionData_Y{partitionIndex} rawData_Y(i)];
			partitionIndex                 = partitionIndex+1;
			if (partitionIndex == 11);
				partitionIndex = 1;
			end;
		end;
		fprintf(':');

		% Perform 10 fittings.
		LSerrorSum = 0;
		Fitting_Y  = cell(10,1);
		Fitting_X  = linspace(minX,maxX,final_length);
		for partitionIndex = 1:10
			% Define training dataset.
			Training_X = [];
			Training_Y = [];
			for i = 1:10
				if (i ~= partitionIndex)
					Training_X = [Training_X partitionData_X{partitionIndex}];
					Training_Y = [Training_Y partitionData_Y{partitionIndex}];
				end;
			end;

			% Define testing dataset.
			Testing_X   = partitionData_X{partitionIndex};
			Testing_Y   = partitionData_Y{partitionIndex};

			% Generate LOWESS fittings of training dataset.
			%    Training_Y_smoothed : fit curve to training dataset at Testing_X coordinates (used to compare fit curve to test dataset).
			%    Fitting_Y           : fit curve to training dataset at Fitting_X coordinates (a linear spacing of 250 points from min to max of raw data), 
			arrayDim = size(Training_X);
			if (arrayDim(1) > arrayDim(2))
				Training_Y_smoothed         = mylowess([Training_X,  Training_Y ],Testing_X,spans(j));
				Fitting_Y{partitionIndex}   = mylowess([Training_X,  Training_Y ],Fitting_X,spans(j));
			else
				Training_Y_smoothed         = mylowess([Training_X', Training_Y'],Testing_X,spans(j));
				Fitting_Y{partitionIndex}   = mylowess([Training_X', Training_Y'],Fitting_X,spans(j));
			end;

			% determine cumulative error between LOWESS fit and testing dataset.
			for i = 1:length(Testing_X)
				LSerrorSum  = LSerrorSum  + (Testing_Y(i)-Training_Y_smoothed(i))^2;
			end;
		end;

		% Perform fitting to all data at current span.
		arrayDim = size(rawData_X);
		if (arrayDim(1) > arrayDim(2))
			allData_Fitting_Y = mylowess([rawData_X, rawData_Y ],Fitting_X,spans(j));
		else
			allData_Fitting_Y = mylowess([rawData_X',rawData_Y'],Fitting_X,spans(j));
		end;
		fitCurves{j}      = allData_Fitting_Y;

%		Fitting_Y_average = zeros(1,final_length);
%		for partitionIndex = 1:10
%			Fitting_Y_average = Fitting_Y_average + Fitting_Y{partitionIndex};
%%			for i = 1:length(Fitting_Y{partitionIndex})
%%			end;
%%			fprintf(['\ntest_' num2str(partitionIndex-1) ' = [' num2str(Fitting_Y{partitionIndex}) ']\n']);
%		end;
%		fitCurves{j} = Fitting_Y_average/10;

		sse(j) = LSerrorSum;
		fprintf(['](error = ' num2str(sse(j)) ')\n']);
	end;
end;


% Find the smoothing value which produces the least error between the LOWESS fit and the raw data.
[minsse,minj] = min(sse);
span          = spans(minj);
%	if (span < 0.3)
%		span = spans(3);
%		minj = 3;
%	end;

X_range       = linspace(minX,maxX,final_length);

newX = X_range;
if ((LOWESS_method == 1) || (LOWESS_method == 2))
	% Generate final fit from LOWESS (on all data) with found best span.
	newY = mylowess([X,Y],X_range,span);
elseif (LOWESS_method == 3)
	% pull the curve fit produced earlier.
	newY = fitCurves{minj};
end;

end
