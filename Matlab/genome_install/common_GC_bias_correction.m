function [output] = common_GC_bias_correction(input)
%% Corrects GC bias per standard display bin in Common_CNV files for analysis pipeline.

num_chrs = length(input);

%% figure out how much data was provided.
dataLength = 0;
for chr = 1:num_chrs
	dataLength = dataLength + length(input{chr});
end;

%% Load GC_ratio data...
X_data   = zeros(1,dataLength);

%% gather CNV data from input...
Y_data   = zeros(1,dataLength);
current  = 0;
for chr = 1:num_chrs;
	for bin = 1:length(input{chr})
		current = current + 1;
		Y_data(current) = input{chr,2}(bin);
	end;
end;

fprintf('\tPreparing for LOWESS fitting : GC_ratio vs CNV data.\n');

%% Perform LOWESS fitting.
fprintf('\tLOWESS fitting to reference data.\n');
[newX2, newY2] = optimize_mylowess(X_data,Y_data,10, 0);
fprintf('\tLOWESS fitting to reference data complete.\n');

Y_target          = 1;
for chr = 1:num_chrs
	Y_fitCurve    = interp1(newX2,newY2,input{chr,2},'spline');
	output{chr,1} = 1;
	output{chr,2} = input{chr,2}./Y_fitCurve*Y_target;
end;
