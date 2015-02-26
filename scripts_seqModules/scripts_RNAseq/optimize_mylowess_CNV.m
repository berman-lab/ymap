function [newX, newY] = optimize_mylowess_CNV(rawData_X,rawData_Y)
addpath('../);

%numDat        = length(rawData_X);
%spans         = linspace(0.01,0.99, numFits);
%sse           = zeros(size(spans));
%cp            = cvpartition(numDat,'k',10);
X             = rawData_X';
Y             = rawData_Y';

%% Attempts LOWESS fitting with [numFits] smoothing values evenly spaced from 0.01 to 0.99.
%for j = 1:numFits
%	fprintf(['\t\tFit #' num2str(j) '/' num2str(numFits) '.  [']);
%	f      = @(train,test) norm(test(:,2) - mylowess(train,test(:,1),spans(j)))^2;
%	sse(j) = sum(crossval(f,[X,Y],'partition',cp));
%	fprintf(']\n');
%end;

% Find the smoothing value which produces the least error between the LOWESS fit and the raw data.
%[minsse,minj] = min(sse);
%span          = spans(minj);
span          = 0.7;
%	fprintf(['X      = ' num2str(X)      '\n']);
%	fprintf(['min(X) = ' num2str(min(X)) '\n']);
%	fprintf(['max(X) = ' num2str(max(X)) '\n']);
X_range       = linspace(min(X),max(X),400);

newX = X_range;
arrayDim = size(X);
if (arrayDim(1) > arrayDim(2))
	newY = mylowess([X, Y],X_range,span);
else
	newY = mylowess([X', Y'],X_range,span);
end;
%newY = mylowess([X,Y],X_range,span);

end
