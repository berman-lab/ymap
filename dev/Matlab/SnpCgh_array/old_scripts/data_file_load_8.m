%% ========================================================================
% Reads in Tab-delimited text format output from microarray image analysis
%    programs.
% The default is to analyze data files produced by BlueFuse.
% Alternate data file formats can be used if control variables are
%    adjusted as needed.
%==========================================================================
function data_file_load_8(file_dir,file_name,design,experiment_name)
% control variables describing microarray analysis package output files.
header           = 46; % lines of header to be skipped before data lines.
Name_column      = 1;  % column containing probe ID/name.
Ch1_column       = 4;  % column containing probe channel 1 data.
Ch2_column       = 5;  % column containing probe channel 1 data.
Ratio_column     = 6;  % column containing probe channel ratio data.
Log2Ratio_column = 7;  % column containing probe channel log2ratio data.

% Control variables for SNP and CGH probe lists.
input_file1   = ['designs/' design '/SNP_probe_list.txt'];
input_file2   = ['designs/' design '/CGH_probe_list.txt'];
details_file  = ['designs/' design '/details.txt'];
fid   = fopen(details_file);
line1 = fgetl(fid);
line2 = fgetl(fid);
fclose(fid);
num_lines1 = str2num(sscanf(line1, '%s',1));
num_lines2 = str2num(sscanf(line2, '%s',1));
%%=========================================================================
final      = [];
raw        = [];

% initializing variables used.
probeset1       = [];
probeset2       = [];
IDset1          = {};
IDset2          = {};

%% load SNP probe list
fid  = fopen(input_file1, 'r');
% analyze SNP text file, line by line.
fprintf(['\nAnalyzing text file: \"' input_file1 '\"\n'])
i = 0;
lines_analyzed = 0;
while not (feof(fid))
    i              = i+1;
    line           = fgetl(fid);
    lines_analyzed = lines_analyzed+1;
    
    % take of interest data fields from each line.
    ProbeID              = sscanf(line, '%s',1);
    % add data fields into structure.
    probeset1(i).probe_ID = ProbeID;
    IDset1{i}             = char(ProbeID);
    
    % interpret probeID to determine probe chromosome number and location.
    probeset1(i).probe_chromosome           = str2double(ProbeID(9));
    probeset1(i).probe_location             = str2double(ProbeID(11:length(ProbeID)-2));
end;
fclose(fid);

%% load CGH probe list
fid2 = fopen(input_file2, 'r');
% analyze CGH text file, line by line.
fprintf(['\nAnalyzing text file: \"' input_file2 '\"\n'])
i = 0;
lines_analyzed = 0;
while not (feof(fid2))
    i              = i+1;
    line           = fgetl(fid);
    lines_analyzed = lines_analyzed+1;
    
    % take of interest data fields from each line.
    ProbeID               = sscanf(line, '%s',1);
    % add data fields into structure.
    probeset2(i).probe_ID = ProbeID;
    IDset2{i}             = char(ProbeID);
    
    % interpret probeID to determine probe chromosome number and location.
    probeset2(i).probe_chromosome           = str2double(ProbeID(9));
    probeset2(i).probe_location             = str2double(ProbeID(11:length(ProbeID)));
end;
fclose(fid2);

SNP_num_probes = length(IDset1);
CGH_num_probes = length(IDset2);
%cleanup.
clear input_file1 input_file2 lines_to_skip* fid* ProbeID* Sequence i j line ans;
clear counter lines_analyzed num_lines*;

% =========================================================================
lines_to_skip  = header; %number of header lines.

%% Process the input files for CGH and SNP information.
%% arbitrary data file format.
% Could be improved by presorting the input file then using
% algorithm for Bluefuse fused output file.
% Load CGH and SNP information in one step... very slow.

% load SNP & CGH information.
fprintf(['[data_file_load_8.m]\n']);
fprintf(['    ' file_dir '/' file_name '\n']);

fid = fopen([file_dir '/' file_name], 'r');
%skip header lines defined by 'lines_to_skip'.
if (lines_to_skip > 0)
    for j = 1:lines_to_skip
        discard = fgetl(fid);
    end;
end;
clear discard;

% analyze text file, line by line.
fprintf(['\nAnalyzing \"', file_name, '\" for SNP info.']);
lines_analyzed = 0;
j = 0;
counter = 0;
rawData = [];
while not (feof(fid))   % build raw data structure.
    counter = counter+1;
    j = j+1;
    line           = fgetl(fid);
    % take of interest data fields from each line.
    if (length(line) > 3)
	lines_analyzed = lines_analyzed+1;
	ProbeID		= sscanf(line, '%s',Name_column);
	for k = 1:length(sscanf(line,'%s',Name_column-1));
	    ProbeID(1) = [];
	end;
	ch1		= sscanf(line, '%s',Ch1_column);
	for k = 1:length(sscanf(line,'%s',Ch1_column-1));
	    ch1(1) = [];
	end;
	ch2		= sscanf(line, '%s',Ch2_column);
	for k = 1:length(sscanf(line,'%s',Ch2_column-1));
	    ch2(1) = [];
	end;
	Ratio		= sscanf(line, '%s',Ratio_column);
	for k = 1:length(sscanf(line,'%s',Ratio_column-1));
	    Ratio(1) = [];
	end;
	Log2Ratio	= sscanf(line, '%s',Log2Ratio_column);
	for k = 1:length(sscanf(line,'%s',Log2Ratio_column-1));
	    Log2Ratio(1) = [];
	end;

	rawData(lines_analyzed).ProbeID	  = ProbeID;
	rawData(lines_analyzed).ch1	  = str2num(ch1);
	rawData(lines_analyzed).ch2	  = str2num(ch2);
	rawData(lines_analyzed).Ratio	  = str2num(Ratio);
	rawData(lines_analyzed).Log2Ratio = str2num(Log2Ratio);
    end;
end;

% Converts loaded data structure into a cell array, which is sorted and then reformated back to the original structure form.
dataFields	= fieldnames(rawData);
dataCell	= struct2cell(rawData);
dataSize	= size(dataCell);
dataMatrix	= reshape(dataCell,dataSize(1),[]);
dataMatrix	= dataMatrix';
dataMatrix	= sortrows(dataMatrix,1);
dataMatrix	= reshape(dataMatrix',dataSize);
sortedData	= cell2struct(dataMatrix,dataFields,1);

% Used to help speed up matching loaded data entries into probeset.
%    Each time a CGH or SNP probe is looked for in the design file, we start from the last found location.
%    This requires the list of experimental probe data to be sorted by probeID, as is done in the block above.
CGH_search_start = 1;
SNP_search_start = 1;
countCGH1 = 0;
countCGH2 = 0;
countSNP1 = 0;
countSNP2 = 0;
countOther = 0;
for i = 1:length(sortedData)
    % Add data fields into structure: very slow step.
    %    CGH goes into probeset2.
    %    SNP goes into probeset1.
    ProbeID	= sortedData(i).ProbeID;
    ch1		= sortedData(i).ch1;
    ch2		= sortedData(i).ch2;
    Ratio	= sortedData(i).Ratio;
    Log2Ratio	= sortedData(i).Log2Ratio;

    if (strcmp(ProbeID(1:3),'SNP') == 1)
	countSNP2=countSNP2+1;
	for k = SNP_search_start:SNP_num_probes
	    if (strcmp(probeset1(k).probe_ID,ProbeID) == 1)
		probeset1(k).probe_ch1       = ch1;
		probeset1(k).probe_ch2       = ch2;
		probeset1(k).probe_Ratio     = Ratio;
		probeset1(k).probe_Log2Ratio = Log2Ratio;
		SNP_search_start = k+1;
		countSNP1=countSNP1+1;
		break;
	    end;
	end;
    elseif (strcmp(ProbeID(1:3),'CGH') == 1)
	countCGH2=countCGH2+1;
	for k = CGH_search_start:CGH_num_probes
	    if (strcmp(probeset2(k).probe_ID,ProbeID) == 1)
		probeset2(k).probe_ch1       = ch1;
		probeset2(k).probe_ch2       = ch2;
		probeset2(k).probe_Ratio     = Ratio;
		probeset2(k).probe_Log2Ratio = Log2Ratio;
		CGH_search_start = k+1;
		countCGH1=countCGH1+1;
		break;
	    end;
	end;
    end;
end;
fclose(fid);
fprintf(['\n Total   = ' num2str(countOther+countCGH1+countCGH2+countSNP1+countSNP2) ]);
fprintf(['\n   CGH   = ' num2str(countCGH1) '/' num2str(countCGH2) '(' num2str(CGH_num_probes) ')' ]);
fprintf(['\n   SNP   = ' num2str(countSNP1) '/' num2str(countSNP2) '(' num2str(SNP_num_probes) ')' ]);
fprintf(['\n   other = ' num2str(countOther) ]);
fprintf('\n');

% save *.MAT file of probeset.
save([file_dir '/' strrep(experiment_name,' ','_') '.' design '.SNP_data.mat'],'probeset1');
save([file_dir '/' strrep(experiment_name,' ','_') '.' design '.CGH_data.mat'],'probeset2');

end
