classdef cvpartition
    %CVPARTITION Create a cross-validation partition for data.
    %   An object of the CVPARTITION class defines a random partition on a
    %   set of data of a specified size.  This partition can be used to
    %   define test and training sets for validating a statistical model
    %   using cross-validation.
    %
    %   C = CVPARTITION(N,'Kfold',K) creates a CVPARTITION object C
    %   defining a random partition for K-fold cross-validation on N
    %   observations. The partition divides N observations into K disjoint
    %   subsamples (folds), chosen randomly but with roughly equal size.
    %   The default value of K is 10.
    %
    %   C = CVPARTITION(GROUP,'Kfold',K) creates a CVPARTITION object C
    %   defining a random partition for a stratified K-fold
    %   cross-validation. GROUP is a vector indicating the class
    %   information for each observation. GROUP can be a categorical
    %   variable, a numeric vector, a string array, or a cell array of
    %   strings. Each subsample has roughly equal size and roughly the same
    %   class proportions as in GROUP. CVPARTITION treats NaNs or empty
    %   strings in GROUP as missing values.
    %
    %   C = CVPARTITION(N,'Holdout',P) creates a CVPARTITION object C
    %   defining a random partition for holdout validation on N
    %   observations. This partition divides N observations into a training
    %   set and a test (holdout) set. P must be a scalar. When 0<P<1,
    %   CVPARTITION randomly selects approximately P*N observations for the
    %   test set. When P is an integer, it randomly selects P observations
    %   for the test set. The default value of P is 1/10.
    %
    %   C = CVPARTITION(GROUP,'Holdout',P) randomly partitions observations
    %   into a training set and a test set with stratification, using the
    %   class information in GROUP, i.e., both training and test sets have
    %   roughly the same class proportions as in GROUP.
    %
    %   C = CVPARTITION(N,'Leaveout') creates an object C defining a random
    %   partition for Leave-one-out cross-validation on N observations.
    %   Leave-one-out is a special case of K-fold in which the number of
    %   folds is equal to the number of observations N.
    %
    %   C = CVPARTITION(N,'Resubstitution') creates a CVPARTITION object C
    %   which does not partition the data. Both the training set and the
    %   test set contain all of the original N observations.
    %
    %   C has the following properties:
    %
    %      Type         The type of validation partition. It is 'kfold',
    %                   'holdout', 'leaveout' or 'resubstitution'.
    %      N            The number of observations (including observations
    %                   with missing GROUP values if GROUP is provided).
    %      NumTestSets  The number of test sets. Its value is the
    %                   number of folds in K-fold and Leave-one-out, 1 in
    %                   Holdout and Resubstitution.
    %      TrainSize    The size of each training set. It is a vector in
    %                   K-fold and Leave-one-out, a scalar in Holdout and
    %                   Resubstitution.
    %      TestSize     The size of each test set. It is a vector in K-fold
    %                   and Leave-one-out or a scalar in Holdout and
    %                   Resubstitution.
    %
    %   Example: Use 10-fold stratified cross-validation to compute the
    %   mis-classification error for CLASSIFY on iris data.
    %
    %     load('fisheriris');
    %     CVO = cvpartition(species,'k',10);
    %     err = zeros(CVO.NumTestSets,1);
    %     for i = 1:CVO.NumTestSets
    %         trIdx = CVO.training(i);
    %         teIdx = CVO.test(i);
    %         ytest = classify(meas(teIdx,:),meas(trIdx,:),species(trIdx,:));
    %         err(i) = sum(~strcmp(ytest,species(teIdx)));
    %     end
    %     cvErr = sum(err)/sum(CVO.TestSize);
    %
    %   See also CVPARTITION/TEST, CVPARTITION/TRAINING, CVPARTITION/REPARTITION, CROSSVAL.
    %   Copyright 2007-2009 The MathWorks, Inc.
    %   $Revision: 1.1.6.8 $  $Date: 2009/07/18 15:55:08 $

    properties(GetAccess = 'public', SetAccess = 'protected')
        Type = '';
        N = [];
        NumTestSets = [];
        TrainSize = [];
        TestSize = [];
    end

    properties(GetAccess = 'protected', SetAccess = 'protected')
        indices = [];
        Group = [];
        holdoutT = [];
    end

    methods(Access = 'public')
        function cv = cvpartition(N,method,T,varargin)
            
            if isempty(varargin)
                s = RandStream.getGlobalStream;
                stdargin = nargin;
            else
                if length(varargin)>1
                    error('stats:cvpartition:cvpartition', ...
                        'CVPARTITION can have at most one optional argument');
                else
                    stdargin = nargin-1;
                    s = varargin{1};
                    if ~isa(s,'RandStream')
                        error('stats:cvpartition:cvpartition', ...
                            'CVPARTITION optional argument must be a RandStream object');
                    end
                end
            end
            
            if stdargin < 2
                error('stats:cvpartition:TooFewInputs',...
                    'At least two arguments are needed.');
            end

            if ischar(method) && size(method,1) == 1
                methodNames = {'kfold','holdout','leaveout','resubstitution'};
                j = strmatch(lower(method),methodNames);
                if length(j) > 1
                    error('stats:cvpartition:AmbiguousMethod', ...
                        'The second input has an ambiguous value:  %s.', method);
                elseif isempty(j)
                    error('stats:cvpartition:UnknownMethod', ...
                        'The second input must be either ''Kfold'', ''Holdout'',''Leaveout'' or ''Resubstitution''.');
                end
            else
                error('stats:cvpartition:InvalidType', ...
                    'The second input must be a string.');
            end

            cv.Type = methodNames{j};

            if isscalar(N)
                if ~isnumeric(N) || N <= 1 || N ~= round(N) || ~isfinite(N)
                    error('stats:cvpartition:BadN',...
                        'The number of observations must be a positive integer greater than one.');
                end
                cv.N = N;
            else
                cv.Group = grp2idx(N);
                cv.N = length(cv.Group); % the number of observations including NaNs
                [ignore,wasnan,cv.Group] = statremovenan(cv.Group);
                hadNaNs = any(wasnan);
                if hadNaNs
                    warning('stats:cvpartition:MissingGroupsRemoved',...
                        'Ignoring rows in GROUP with missing values.');
                    if length (cv.Group) <= 1
                        error('stats:cvpartition:BadN',...
                            'The number of rows in GROUP with non-missing values must be greater than one.');
                    end
                end
            end

            dftK = 10; % the default number of subsamples(folds) for Kfold
            P  = 1/10; % the default holdout ratio

            switch cv.Type
                case 'kfold'
                    if stdargin == 2 || isempty(T)
                        T = dftK;
                    elseif ~isscalar(T) || ~isnumeric(T) || T <= 1 || ...
                            T ~= round(T) || ~isfinite(T)
                        error('stats:cvpartition:BadK',...
                            'For K-fold, the number of folds must be an integer greater than one.');
                    end

                    if  isempty(cv.Group) && T > cv.N
                        warning('stats:cvpartition:KfoldTooLarge',...
                            ['The number of folds K is greater than ',...
                            'the number of observations N. ',...
                            'K will be set to the value of N.']);
                        T = cv.N;
                    elseif ~isempty(cv.Group) && T > length(cv.Group)
                        warning('stats:cvpartition:KfoldTooLarge',...
                            ['The number of folds K is greater than the number of observations with non-missing'...
                            'group values. K will be set to the number of '...
                            'observations with non-missing group values.']);
                        T = length(cv.Group);
                    end

                    cv.NumTestSets = T; %set the number of fold
                    cv = cv.rerandom(s);

                case 'leaveout'
                    if stdargin > 2 && T ~= 1
                        error('stats:cvpartition:UnsupportedLeaveout',...
                            'Currently the only supported value of leave-out is 1.');
                    end

                    if isempty(cv.Group)
                        cv.NumTestSets = cv.N;
                    else
                        cv.NumTestSets = length(cv.Group);
                    end

                    [~,cv.indices] = sort(rand(s,cv.NumTestSets,1));

                    cv.TrainSize = (cv.NumTestSets-1) * ones(1,cv.NumTestSets);
                    cv.TestSize = ones(1,cv.NumTestSets);

                case 'resubstitution'
                    if stdargin > 2 
                        error('stats:cvpartition:UnsupportedCV',...
                            'Only two inputs are accepted when the second input is ''resubstitution''.');
                    end

                    if isempty(cv.Group)
                        numObs = N;
                    else
                        numObs = length(cv.Group);
                    end

                    cv.indices = (1: numObs)';
                    cv.NumTestSets = 1;
                    cv.TrainSize =  numObs;
                    cv.TestSize =  numObs;

                case 'holdout'
                    if stdargin == 2 || isempty(T)
                        T = P;
                    elseif ~isscalar(T) || ~ isnumeric(T) || T <= 0 || ~isfinite(T)
                        error('stats:cvpartition:BadP',...
                            'For holdout, the third input must be a non-negative numeric scalar.');
                    end

                    if T >= 1 %hold-T observations out
                        if T ~=round(T)
                            error('stats:cvpartition:BadP',...
                                ['For holdout, the third input must be either a scalar in the range of ',...
                                '(0,1) or a positive integer.']);
                        end
                        if isempty(cv.Group)
                            if T >= cv.N
                                error('stats:cvpartition:BadP',...
                                    ['The number of observations for testing (holdout) should be ',...
                                    'smaller than the number of observations.']);
                            end
                        else
                            if T>= length(cv.Group)
                                error('stats:cvpartition:BadP',...
                                    ['The number of observations for testing (holdout) should be ',...
                                    'smaller than the number of observations with non-missing group values.']);
                            end
                        end
                    else
                        if (isempty(cv.Group) && floor(cv.N *T) == 0) ||...
                                (~isempty(cv.Group) && floor(length(cv.Group) * T) == 0)
                            error('stats:cvpartition:BadP',...
                                'P is too small to have a non-empty test set.');

                        end
                    end

                    cv.holdoutT = T;
                    cv = cv.rerandom(s);
                    cv.NumTestSets = 1;
            end

            %add NaNs back
            if ~isempty(cv.Group) && hadNaNs
                [cv.indices, cv.Group] =...
                    statinsertnan(wasnan, cv.indices, cv.Group);
            end
        end % cvpartition constructor

        function  cv = repartition(cv,varargin)
        %REPARTITION Rerandomize a cross-validation partition. 
        %   D = REPARTITION(C) creates a new random cross-validation partition D
        %   of the same type and size as C.  Use REPARTITION to perform multiple
        %   Monte-Carlo repetitions of cross-validation.
        %   D = REPARTITION(C,S) uses the RandStream object S as its
        %   random number generator.
        %   
        %   See also CVPARTITION.

            if isempty(varargin)
                s = RandStream.getGlobalStream;
            else
                if length(varargin)>1
                    error('stats:cvpartition:Repartition', ...
                        'REPARTITION can only at most one optional argument');
                else
                    s = varargin{1};
                    if ~isa(s,'RandStream')
                        error('stats:cvpartition:Repartition', ...
                            'REPARTITION optional argument must be a RandStream object');
                    end
                end
            end
        
            if strcmp(cv.Type,'resubstitution')
                warning('stats:cvpartition:RepartNone',...
                    'C keeps the same test and training sets as before if the ''Type'' is ''resubstitution''.');
                return;
            end
            %remove NaNs from cv.Group
            if ~isempty(cv.Group)
                [ignore,wasnan,cv.Group] = statremovenan(cv.Group);
                hadNaNs = any(wasnan);
            end

            %  regenerate the data partition
            cv = cv.rerandom(s);

            %add NaNs back into cv.indices and cv.Group
            if ~isempty(cv.Group) && hadNaNs
                [cv.indices, cv.Group] =...
                    statinsertnan(wasnan, cv.indices, cv.Group);
            end
        end % repartition
       
        function trainIndices = training(cv,i)
        %TRAINING Training set for a cross-validation partition.
        %   TRIDX = TRAINING(C) returns a logical vector TRIDX that selects
        %   the observations in the training set for the hold-out
        %   cross-validation partition C.  C may also be a partition for
        %   resubstitution, in which case TRIDX is a logical vector that
        %   selects all observations.
        %
        %   TRIDX = TRAINING(C,I) returns a logical vector TRIDX that selects
        %   the observations in the I-th training set for a K-fold or
        %   leave-one-out cross-validation partition C.  In K-fold
        %   cross-validation, C divides a data set into K disjoint folds with
        %   roughly equal size.  The I-th training set consists of all
        %   observations not contained in the I-th fold.  In leave-one-out
        %   cross-validation, the I-th training set consists of the entire
        %   data set except the I-th observation.
        %
        %   See also CVPARTITION, CVPARTITION/TEST.
            switch cv.Type
                case {'kfold', 'leaveout'}
                    if nargin ~= 2
                        error('stats:cvpartition:WrongNumInputs', ...
                            'Requires two input arguments.');
                    end
                    checkindex(i,cv.NumTestSets);
                    trainIndices = (cv.indices ~= i & ~isnan(cv.indices));

                case 'holdout'
                    if nargin == 2 && i~=1
                        error('stats:cvpartition:InvalidIndex',...
                            'I must be 1 for Holdout validation.');
                    end
                    trainIndices = (cv.indices == 1);
                case 'resubstitution'
                    if nargin == 2 && i~= 1
                        error('stats:cvpartition:InvalidIndex',...
                            'I must be 1 if cv.Type is ''resubstitution''.');
                    end
                    trainIndices = ~isnan(cv.indices);
            end
        end

        function testIndices = test(cv,i)
        %TEST Test set for a cross-validation partition.
        %   TEIDX = TEST(C) returns a logical vector TEIDX that selects the
        %   observations in the test set for the hold-out cross-validation
        %   partition C.  C may also be a partition for resubstitution, in
        %   which case TEIDX is a logical vector that selects all
        %   observations.
        %
        %   TEIDX = TEST(C,I) returns a logical vector TEIDX that selects the
        %   observations in the I-th test set for a K-fold or leave-one-out
        %   cross-validation partition C.  In K-fold cross-validation, C
        %   divides a data set into K disjoint folds with roughly equal size.
        %   The I-th test set consists of the I-th fold.  In leave-one-out
        %   cross-validation, the I-th test set consists of the I-th
        %   observation.
        %
        %   See also CVPARTITION, CVPARTITION/TRAINING.
            switch cv.Type
                case {'kfold','leaveout'}
                    if nargin ~= 2
                        error('stats:cvpartition:WrongNumInputs', ...
                            'Requires two input arguments.');
                    end
                    checkindex(i,cv.NumTestSets);

                    testIndices = (cv.indices == i);

                case 'holdout'
                    if nargin == 2 && i~= 1
                        error('stats:cvpartition:InvalidIndex',...
                            'I must be 1 for Holdout validation.');
                    end
                    testIndices = (cv.indices == 2);
                case 'resubstitution'
                    if nargin == 2 && i ~= 1
                        error('stats:cvpartition:InvalidIndex',...
                            'I must be 1 if cv.Type is ''resubstitution''.');
                    end
                    testIndices = ~isnan(cv.indices);
            end
        end

        % Display methods
        function display(cv)
            isLoose = strcmp(get(0,'FormatSpacing'),'loose');

            objectname = inputname(1);
            if isempty(objectname)
                objectname = 'ans';
            end

            if (isLoose)
                fprintf('\n');
            end
            fprintf('%s = \n', objectname);
            disp(cv);
        end
        function disp(cv)
            isLoose = strcmp(get(0,'FormatSpacing'),'loose');

            if (isLoose)
                fprintf('\n');
            end
            switch cv.Type
                case 'kfold'
                    disp('K-fold cross validation partition');
                case 'holdout'
                    disp('Hold-out cross validation partition');
                case 'leaveout'
                    disp('Leave-one-out cross validation partition');
                case 'resubstitution'
                    disp ('Resubstitution (no partition of data)');
            end
            disp(['             N: ' num2str(cv.N)]);
            disp(['   NumTestSets: ' num2str(cv.NumTestSets)]);
            Ndisp = 10;
            if cv.NumTestSets <= Ndisp
                disp(['     TrainSize: ' num2str(cv.TrainSize )]);
                disp(['      TestSize: ' num2str(cv.TestSize )]);
            else
                disp(['     TrainSize: ' num2str(cv.TrainSize(1:Ndisp)), ' ...']);
                disp(['      TestSize: ' num2str(cv.TestSize(1:Ndisp)), ' ...']);
            end
            %             end
        end
    end % public methods block


    methods(Access = 'private')
        %re-generate the data partition using the RandStream object s
        function cv = rerandom(cv,s)

            switch cv.Type
                case 'kfold'
                    if isempty(cv.Group)
                        if cv.NumTestSets == cv.N
                            %special case of K-fold -- loocv
                            [~,cv.indices] = sort(rand(s,cv.NumTestSets,1)); % randperm 
                            cv.TestSize = ones(1,cv.NumTestSets);
                        else
                            cv.indices = kfoldcv(cv.N,cv.NumTestSets,s);
                            cv.TestSize = accumarray(cv.indices,1)';
                        end
                    else
                        if cv.NumTestSets == length(cv.Group)
                            %special case of K-fold -- loocv
                            [~,cv.indices] = sort(rand(s,cv.NumTestSets,1)); % randperm
                            cv.TestSize = ones(1,cv.NumTestSets);
                        else
                            cv.indices = stra_kfoldcv(cv.Group,cv.NumTestSets,s);
                            cv.TestSize = accumarray(cv.indices,1)';
                        end
                    end

                    cv.TrainSize = size(cv.indices,1) - cv.TestSize;

                case 'holdout'
                    if cv.holdoutT >= 1
                        if isempty(cv.Group)
                            cv.indices = holdoutcv(cv.N,cv.holdoutT,s);
                            cv.TestSize = cv.holdoutT;
                            cv.TrainSize = cv.N-cv.TestSize;
                        else
                            cv.indices = stra_holdoutcv(cv.Group, cv.holdoutT/length(cv.Group), s);
                            cv.TestSize = sum(cv.indices == 2);
                            cv.TrainSize = sum(cv.indices == 1);
                        end
                    else %hold cv.holdoutT*N out
                        if isempty(cv.Group)
                            cv.TestSize = floor(cv.N * cv.holdoutT);
                            cv.TrainSize = cv.N-cv.TestSize;
                            cv.indices = holdoutcv(cv.N,cv.TestSize,s);
                        else
                            cv.indices = stra_holdoutcv(cv.Group,cv.holdoutT,s);
                            cv.TestSize = sum(cv.indices == 2);
                            cv.TrainSize = sum(cv.indices == 1);
                        end
                    end

                case 'leaveout'
                    [~,cv.indices] = sort(rand(s,cv.NumTestSets,1));
            end
        end
    end % private methods block
    
    methods(Hidden = true)
        function b = fieldnames(a)
            b = properties(a);
        end
        
        % Methods that we inherit, but do not want
        function a = fields(varargin),     throwUndefinedError(); end
        function a = ctranspose(varargin), throwUndefinedError(); end
        function a = transpose(varargin),  throwUndefinedError(); end
        function a = permute(varargin),    throwUndefinedError(); end
        function a = reshape(varargin),    throwUndefinedError(); end
        function a = cat(varargin),        throwNoCatError(); end
        function a = horzcat(varargin),    throwNoCatError(); end
        function a = vertcat(varargin),    throwNoCatError(); end
    end
    methods(Hidden = true, Static = true)
        function a = empty(varargin)
            error(['stats:' mfilename ':NoEmptyAllowed'], ...
                  'Creation of empty %s objects is not allowed.',upper(mfilename));
        end
    end
   
end % classdef

function throwNoCatError()
error(['stats:' mfilename ':NoCatAllowed'], ...
      'Concatenation of %s objects is not allowed.  Use a cell array to contain multiple objects.',upper(mfilename));
end

function throwUndefinedError()
st = dbstack;
name = regexp(st(2).name,'\.','split');
error(['stats:' mfilename ':UndefinedFunction'], ...
      'Undefined function or method ''%s'' for input arguments of type ''%s''.',name{2},mfilename);
end

%----------------------------------------------------
%stratified k-fold cross-validation
function cvid=stra_kfoldcv(group,nfold,s)
size_groups = accumarray(group(:),1);
if any(size_groups < nfold & size_groups > 0)
    warning('stats:cvpartition:TestZero',...
        'One or more folds do not contain points from all the groups.');
end
N = size(group,1);
cvid = 1 + mod((1:N)',nfold);
idrand = group + rand(s,N,1);
[ignore,idx] = sort(idrand);
cvid = cvid(idx);
end

%----------------------------------------------------
%kfold cross-validation without stratification
function cvid = kfoldcv(N,nfold,s)
cvid = 1 + mod((1:N)',nfold);
[~,indices] = sort(rand(s,1,N)); % randperm
cvid = cvid(indices);
end

%-----------------------------------------------------
%holdout without stratification
function  idx= holdoutcv(N,num_test,s)
idx = 2*ones(N,1);
idx(1:N-num_test) = 1;
[~,indices] = sort(rand(s,1,N)); % randperm
idx = idx(indices);
end

%-----------------------------------------------------
%stratified holdout
function idx = stra_holdoutcv(group,test_ratio,s)
N = length(group);
size_groups = accumarray(group(:),1);
num_test = floor(size_groups * test_ratio);

test_diff = floor(N * test_ratio) - sum(num_test);
%add 1 for groups which are not in the test set
if any(num_test == 0)
    v=(num_test == 0);
    v(cumsum(v) > test_diff) = false;
    num_test(v) = num_test(v) + 1;
    test_diff = test_diff - sum(v);
end


if test_diff > 0
    ng= numel(size_groups);
    wasfull  =(num_test == size_groups);
    full_len = sum(wasfull);
    add = [ones(test_diff,1);zeros(ng - test_diff - full_len,1)];
    [~,indices] =  sort(rand(s,1,(ng-full_len))); % randperm
    add = add(indices);
    x = zeros(size(wasfull));
    x(~wasfull,:) = add;
    num_test = num_test + x;

end

if any(num_test == 0)
    warning('stats:cvpartition:testZero',...
        'The test set does not contain points from all groups.');
end

if any(num_test == size_groups)
    warning('stats:cvpartition:testZero',...
        'The training set does not contain points from all groups.');
end

idx = 2*ones(N,1);
for i = 1:numel(size_groups)
    g_idx = find(group == i);
    idx(g_idx(1:size_groups(i)-num_test(i))) = 1;
    [~,indices] = sort(rand(s,1,(size_groups(i)))); % randperm
    idx(g_idx) = idx(g_idx( indices ));
end

end

%-----------------------------
function checkindex(i,imax)
if ~(isnumeric(i) && isscalar(i) && i == round(i) && 1 <= i && i <= imax)
    error('stats:cvpartition:InvalidIndex', ...
        'Index must be a positive integer less than or equal to %d.',imax);
end
end

