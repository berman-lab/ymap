function [newX, newY] = optimize_mylowess(rawData_X,rawData_Y, numFits)

numDat        = length(rawData_X);
spans         = linspace(0.01,0.99, numFits);
sse           = zeros(size(spans));
cp            = cvpartition(numDat,'k',10);
X             = rawData_X';
Y             = rawData_Y';

% % Attempts LOWESS fitting with [numFits] smoothing values evenly spaced from 0.01 to 0.99.
% %for j = 1:numFits
% %    fprintf(['\t\tFit #' num2str(j) '/' num2str(numFits) '.  [']);
% %    f      = @(train,test) norm(test(:,2) - mylowess(train,test(:,1),spans(j)))^2;
% %    sse(j) = sum(crossval(f,[X,Y],'partition',cp));
% %    fprintf(']\n');
% %end;


% Attempts LOWESS fitting with [numFits] smoothing values evenly spaced from 0.01 to 0.99.
for j = 1:numFits
    fprintf(['\t\tFit #' num2str(j) '/' num2str(numFits) '.  [']);
    % perform LOWESS fitting with current span.
    Ys = mylowess([X,Y],X,spans(j));
    fprintf(':');
    % determine error between LOWESS fit and dataset.
    LSerrorSum = 0;
    for i = 1:length(X)
        LSerrorSum = LSerrorSum + (Y(i)-Ys(i))^2;
    end;
    sse(j) = LSerrorSum;
    fprintf(']\n');
end;


% Attempts LOWESS fitting with [numFits] smoothing values evenly spaced from 0.01 to 0.99, using 10-fold cross-validation of randomly partitioned data.
%for j = 1:numFits
%	% Randomly sort the input data into 10 partitions.
%	randIndex      = randperm(length(rawData_X));
%	partitionDataX = cell(1,10);
%	partitionDataY = cell(1,10);
%	partitionIndex = 1;
%	for i = 1:length(rawData_X);
%		partitionDataX{partitionIndex} = [partitionDataX{partitionIndex} rawData_X(i)];
%		partitionDataY{partitionIndex} = [partitionDataY{partitionIndex} rawData_Y(i)];
%		partitionIndex                 = partitionIndex+1;
%		if (partitionIndex == 11);    partitionIndex = 1;   end;
%	end;
%
%	LSerrorSum = 0;
%	for partitionIndex = 1:10
%		OtherX = [];
%		OtherY = [];
%		for i = 1:10
%			if (i ~= partitionIndex)
%				OtherX = [OtherX partitionDataX{partitionIndex}];
%				OtherY = [OtherY partitionDataY{partitionIndex}];
%			end;
%		end;
%		OtherXY = [OtherX,OtherY]';
%		ThisX   = partitionDataX{partitionIndex};
%		ThisY   = partitionDataY{partitionIndex};
%		ThisXY  = [ThisX,ThisY]';
%		ThisYs  = mylowess([OtherX,OtherY],ThisX,spans(j));
%
%		% determine cumulative error between LOWESS fit and dataset.
%		for i = 1:length(ThisX)
%			LSerrorSum = LSerrorSum + (ThisY(i)-ThisYs(i))^2;
%		end;
%	end;
%	sse(j) = LSerrorSum/10;
%end;


% Find the smoothing value which produces the least error between the LOWESS fit and the raw data.
[minsse,minj] = min(sse);
span          = spans(minj);
X_range       = linspace(min(X),max(X),400);

newX = X_range;
newY = mylowess([X,Y],X_range,span);

end
