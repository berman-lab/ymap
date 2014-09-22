function [newX, newY] = optimize_mylowess(rawData_X,rawData_Y, numFits, maxX, minSpanKey)

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
		fprintf(['\t\t\tFit type0 #' num2str(j) '/' num2str(numFits) '.\t[']);
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
		fprintf(['\t\t\tFit type1 #' num2str(j) '/' num2str(numFits) '.\t[']);
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
		fprintf(['\t\t\tFit type2 #' num2str(j) '/' num2str(numFits) '.\t[']);
		% Randomly sort the input data into 10 partitions.
		randIndex      = randperm(length(rawData_X));      % random order of length the length of the data.
		partitionDataX = cell(1,10);
		partitionDataY = cell(1,10);
		partitionIndex = 1;
		for i = 1:length(rawData_X);
			partitionDataX{partitionIndex} = [partitionDataX{partitionIndex} rawData_X(i)];
			partitionDataY{partitionIndex} = [partitionDataY{partitionIndex} rawData_Y(i)];
			partitionIndex                 = partitionIndex+1;
			if (partitionIndex == 11);
				partitionIndex = 1;
			end;
		end;
		fprintf(':');
		LSerrorSum  = 0;
		fittingY = cell(10,1);
		fittingX = linspace(minX,maxX,final_length);
		for partitionIndex = 1:10
			OtherX = [];
			OtherY = [];
			for i = 1:10
				if (i ~= partitionIndex)
					OtherX = [OtherX partitionDataX{partitionIndex}];
					OtherY = [OtherY partitionDataY{partitionIndex}];
				end;
			end;
			ThisX   = partitionDataX{partitionIndex};
			ThisY   = partitionDataY{partitionIndex};
			arrayDim = size(OtherX);
			if (arrayDim(1) > arrayDim(2))
				ThisYs                   = mylowess([OtherX, OtherY],ThisX,spans(j));
				fittingY{partitionIndex} = mylowess([ThisX,  ThisY ],fittingX,spans(j));
			else
				ThisYs                   = mylowess([OtherX', OtherY'],ThisX,spans(j));
				fittingY{partitionIndex} = mylowess([ThisX',  ThisY' ],fittingX,spans(j));
			end;
			% determine cumulative error between LOWESS fit and dataset.
			for i = 1:length(ThisX)
				LSerrorSum  = LSerrorSum  + (ThisY(i)-ThisYs(i))^2;
			end;
		end;

		fittingY_average = zeros(1,final_length);
		for partitionIndex = 1:10
			fittingY_average = fittingY_average + fittingY{partitionIndex};
%			for i = 1:length(fittingY{partitionIndex})
%			end;
%			fprintf(['\ntest_' num2str(partitionIndex-1) ' = [' num2str(fittingY{partitionIndex}) ']\n']);
		end;
		fitCurves{j} = fittingY_average/10;

		sse(j) = LSerrorSum;
		fprintf(['](error = ' num2str(sse(j)) ')\n']);
	end;
end;


% Find the smoothing value which produces the least error between the LOWESS fit and the raw data.
[minsse,minj] = min(sse);
if (minj < minSpanKey)
	minj = minSpanKey;
end;
span          = spans(minj);
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
