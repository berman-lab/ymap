%% ========================================================================
% Reads in Tab-delimited text format output from microarray image analysis
%    programs.
% The default is to analyze data files produced by BlueFuse.
% Alternate data file formats can be used if control variables are
%    adjusted as needed.
%==========================================================================
function data_file_load_online(output_file_dir,file_dir,file_name,design,experiment_name)
if (exist([output_file_dir '/' experiment_name '.' design '.CGH_data.mat'], 'file') == 0) && ...
   (exist([output_file_dir '/' experiment_name '.' design '.SNP_data.mat'], 'file') == 0)
	fprintf('\n## Processing input data file into CGH and SNP data structures.\n');

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
	probeset1       = [];   % stores SNP data.
	probeset2       = [];   % stores CGH data.
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
	    lineData       = fgetl(fid);
	    lines_analyzed = lines_analyzed+1;

	    % take of interest data fields from each line.
	    ProbeID              = sscanf(lineData, '%s',1);
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
	    lineData       = fgetl(fid);
	    lines_analyzed = lines_analyzed+1;

	    % take of interest data fields from each line.
	    ProbeID               = sscanf(lineData, '%s',1);
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
	clear input_file1 input_file2 lines_to_skip* fid* ProbeID* Sequence i j lineData ans;
	clear lines_analyzed num_lines*;

	% =========================================================================
	%% Process the input files for CGH and SNP information.
	%% arbitrary data file format.
	% Could be improved by presorting the input file then using
	% algorithm for Bluefuse fused output file.
	% Load CGH and SNP information in one step... very slow.

	% load SNP & CGH information.
	fprintf('\n');
	fprintf(['[data_file_load_online.m] : input file.\n']);
	fprintf(['    ' file_dir '/' file_name '\n']);

	currentDir = pwd;
	fprintf('\n');
	fprintf(['[data_file_load_online.m] : current folder.\n']);
	fprintf(['    ' currentDir '\n']);
	fprintf('\n');

	% analyze text file, line by line.
	fprintf(['## Analyzing \"', file_name, '\" for CGH & SNP info.\n']);

	%% Open pre-processed files for CGH data.
	CGHfid = fopen([file_dir '/CGH_rows.xls'], 'r');
	fprintf(['##\t' file_dir '/CGH_rows.xls\n##\tfid = ' num2str(CGHfid) '\n']);

	fprintf('## loading CGH subfile into [MATLAB] data structure.\n');
	lines_analyzed = 0;
	rawCghData     = [];
	while not (feof(CGHfid))   % build data structure for raw CGH
		lineData       = fgetl(CGHfid);
		% take of interest data fields from each line.
	    lines_analyzed = lines_analyzed+1;
	    ProbeID    = sscanf(lineData, '%s',Name_column);         for k = 1:length(sscanf(lineData,'%s',Name_column-1));         ProbeID(1) = [];      end;
	    ch1        = sscanf(lineData, '%s',Ch1_column);          for k = 1:length(sscanf(lineData,'%s',Ch1_column-1));          ch1(1) = [];          end;
	    ch2        = sscanf(lineData, '%s',Ch2_column);          for k = 1:length(sscanf(lineData,'%s',Ch2_column-1));          ch2(1) = [];          end;
	    Ratio      = sscanf(lineData, '%s',Ratio_column);        for k = 1:length(sscanf(lineData,'%s',Ratio_column-1));        Ratio(1) = [];        end;
	    Log2Ratio  = sscanf(lineData, '%s',Log2Ratio_column);    for k = 1:length(sscanf(lineData,'%s',Log2Ratio_column-1));    Log2Ratio(1) = [];    end;
	    rawCghData(lines_analyzed).ProbeID         = ProbeID;
	    rawCghData(lines_analyzed).ch1             = str2num(ch1);
	    rawCghData(lines_analyzed).ch2             = str2num(ch2);
	    rawCghData(lines_analyzed).probe_Ratio     = str2num(Ratio);
	    rawCghData(lines_analyzed).probe_Log2Ratio = str2num(Log2Ratio);
	end;
	fclose(CGHfid);

	%% Open pre-processed files for SNP data.
	SNPfid = fopen([file_dir '/SNP_rows.xls'], 'r');
	fprintf(['##\t' file_dir '/SNP_rows.xls\n##\tfid = ' num2str(SNPfid) '\n']);

	fprintf('## loading SNP subfile into [MATLAB] data structure.\n');
	lines_analyzed = 0;
	rawSnpData     = [];
	while not (feof(SNPfid))   % build data structure for raw SNP data.
	    lineData       = fgetl(SNPfid);
		lines_analyzed = lines_analyzed+1;
		ProbeID    = sscanf(lineData, '%s',Name_column);         for k = 1:length(sscanf(lineData,'%s',Name_column-1));         ProbeID(1) = [];      end;
		ch1        = sscanf(lineData, '%s',Ch1_column);          for k = 1:length(sscanf(lineData,'%s',Ch1_column-1));          ch1(1) = [];          end;
		ch2        = sscanf(lineData, '%s',Ch2_column);          for k = 1:length(sscanf(lineData,'%s',Ch2_column-1));          ch2(1) = [];          end;
		Ratio      = sscanf(lineData, '%s',Ratio_column);        for k = 1:length(sscanf(lineData,'%s',Ratio_column-1));        Ratio(1) = [];        end;
		Log2Ratio  = sscanf(lineData, '%s',Log2Ratio_column);    for k = 1:length(sscanf(lineData,'%s',Log2Ratio_column-1));    Log2Ratio(1) = [];    end;
		rawSnpData(lines_analyzed).ProbeID         = ProbeID;
		rawSnpData(lines_analyzed).ch1             = str2num(ch1);
		rawSnpData(lines_analyzed).ch2             = str2num(ch2);
		rawSnpData(lines_analyzed).probe_Ratio     = str2num(Ratio);
		rawSnpData(lines_analyzed).probe_Log2Ratio = str2num(Log2Ratio);
	end;
	fclose(SNPfid);

	% Used to help speed up matching loaded data entries into probeset.
	%    Each time a CGH or SNP probe is looked for in the design file, we start from the last found location.
	%    This requires the list of experimental probe data to be sorted by probeID, as is done in the block above.
	CGH_search_start = 1;
	countCGH1 = 0;
	countCGH2 = 0;
	fprintf(['## matching CGH data entries with probe designs. (' num2str(length(rawCghData)) '): ']);
	for i = 1:length(rawCghData)
		% Add data fields into structure (very slow step) : CGH goes into probeset2.
		ProbeID     = rawCghData(i).ProbeID;
	    ch1         = rawCghData(i).ch1;
	    ch2         = rawCghData(i).ch2;
	    Ratio       = rawCghData(i).probe_Ratio;
	    Log2Ratio   = rawCghData(i).probe_Log2Ratio;
		countCGH2=countCGH2+1;
		for k = CGH_search_start:CGH_num_probes
			% fprintf(['\n' strrep(ProbeID,'v1','') ' : ' probeset2(k).probe_ID]);
	        if (strcmp(probeset2(k).probe_ID,strrep(ProbeID,'v1','')) == 1)
	        	probeset2(k).probe_ch1       = ch1;
	        	probeset2(k).probe_ch2       = ch2;
	        	probeset2(k).probe_Ratio     = Ratio;
	        	probeset2(k).probe_Log2Ratio = Log2Ratio;
	        	CGH_search_start = k+1;
	        	countCGH1=countCGH1+1;
	        	break;
			end;
		end;
		if (mod(i,1000) == 0)
			fprintf(num2str(i));
		elseif (mod(i,100) == 0)
			fprintf('.')
	    end;
	end;

	SNP_search_start = 1;
	countSNP1 = 0;
	countSNP2 = 0;
	fprintf(['\n## matching SNP data entries with probe designs. (' num2str(length(rawSnpData)) '): ']);
	for i = 1:length(rawSnpData)
	    % Add data fields into structure (very slow step) :  SNP goes into probeset1.
	    ProbeID		= rawSnpData(i).ProbeID;
	    ch1			= rawSnpData(i).ch1;
	    ch2			= rawSnpData(i).ch2;
	    Ratio		= rawSnpData(i).probe_Ratio;
	    Log2Ratio	= rawSnpData(i).probe_Log2Ratio;

		countSNP2=countSNP2+1;
		for k = SNP_search_start:SNP_num_probes
			% fprintf(['\n' strrep(ProbeID,'v1','') ' : ' probeset1(k).probe_ID]);
		    if (strcmp(probeset1(k).probe_ID,strrep(ProbeID,'v1','')) == 1)
				probeset1(k).probe_ch1       = ch1;
				probeset1(k).probe_ch2       = ch2;
				probeset1(k).probe_Ratio     = Ratio;
				probeset1(k).probe_Log2Ratio = Log2Ratio;
				SNP_search_start = k+1;
				countSNP1=countSNP1+1;
				break;
		    end;
		end;
		if (mod(i,1000) == 0)
			fprintf(num2str(i));
		elseif (mod(i,100) == 0)
			fprintf('.');
		end;
	end;
	fprintf('\n');

	fprintf(['\n Total   = ' num2str(countCGH1+countCGH2+countSNP1+countSNP2) ]);
	fprintf( '\n   key   = matches/data (designed)');
	fprintf(['\n   CGH   = ' num2str(countCGH1) '/' num2str(countCGH2) '(' num2str(CGH_num_probes) ')' ]);
	fprintf(['\n   SNP   = ' num2str(countSNP1) '/' num2str(countSNP2) '(' num2str(SNP_num_probes) ')' ]);
	fprintf('\n');

	% output current and file output locations to log file.
	fprintf(['\nactive directory = ' pwd]);
	fprintf(['\noutput_file_dir  = ' output_file_dir]);
	fprintf('\n');

	% save *.MAT files of probeset.
	save([output_file_dir '/' strrep(experiment_name,' ','_') '.' design '.SNP_data.mat'],'probeset1');
	save([output_file_dir '/' strrep(experiment_name,' ','_') '.' design '.CGH_data.mat'],'probeset2');
else
	fprintf('\n## Input data file already processed into CGH and SNP data structures.\n');
	fprintf(['##\t' output_file_dir strrep(experiment_name,' ','_') '.' design '.SNP_data.mat\n']);
	fprintf(['##\t' output_file_dir strrep(experiment_name,' ','_') '.' design '.CGH_data.mat\n\n']);
end;
end
