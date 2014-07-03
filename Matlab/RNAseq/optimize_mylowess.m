function [newX, newY] = optimize_mylowess(rawData_X,rawData_Y, numFits)

numDat        = length(rawData_X);
spans         = linspace(0.01,0.99, numFits);
sse           = zeros(size(spans));
cp            = cvpartition(numDat,'k',10);
X             = rawData_X';
Y             = rawData_Y';

% Attempts LOWESS fitting with [numFits] smoothing values evenly spaced from 0.01 to 0.99.
for j = 1:numFits
    fprintf(['\t\tFit #' num2str(j) '/' num2str(numFits) '.  [']);
    f      = @(train,test) norm(test(:,2) - mylowess(train,test(:,1),spans(j)))^2;
    sse(j) = sum(crossval(f,[X,Y],'partition',cp));
    fprintf(']\n');
end;

% Find the smoothing value which produces the least error between the LOWESS fit and the raw data.
[minsse,minj] = min(sse);
span          = spans(minj);
X_range       = linspace(min(X),max(X),400);

newX = X_range;
newY = mylowess([X,Y],X_range,span);

end
