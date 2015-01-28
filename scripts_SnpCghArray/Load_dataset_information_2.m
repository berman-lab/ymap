function [segmental_aneuploidy] = Load_dataset_information_1(projectName,workingDir)
% Load centromere definition file.
%    This is text file containing one header line and two columns.
%    The two columns hold the start and end bp for the centromeres, with
%       respect to each chromosome.

%%=========================================================================
% Load Common_ChARM files for project.
%--------------------------------------------------------------------------
% CNV data   : Common_ChARM.mat           : 'segmental_aneuploidy'
% SNP ratios : Common_ChARM_SNPratios.mat : 'segmental_aneuploidy'
%--------------------------------------------------------------------------
if (exist([workingDir 'Common_ChARM.mat'],'file') ~= 0) && (exist([workingDir 'Common_ChARM_SNPratios.mat'],'file') ~= 0)
	%%%% Load ChARM:CNV breakpoints.
	dataFile1 = [workingDir 'Common_ChARM.mat'];
	fprintf(['\nLoading Common_ChARM.mat file for "' projectName '" : ' dataFile1 '\n']);
	load(dataFile1);
	segmental_aneuploidy_1 = segmental_aneuploidy;

	for i = 1:length(segmental_aneuploidy)
		var1 = segmental_aneuploidy(i).chr;
		var2 = segmental_aneuploidy(i).break;
		fprintf(['a1.breakpoint = chr' num2str(var1) ':' num2str(var2) '\n']);
	end;

	%%%% Load ChARM:SNP-ratio breakpoints.
	dataFile2 = [workingDir 'Common_ChARM_SNPratios.mat'];
	fprintf(['\nLoading Common_ChARM_SNPratios.mat file for "' projectName '" : ' dataFile2 '\n']);
	load(dataFile2);
	segmental_aneuploidy_2 = segmental_aneuploidy;

	for i = 1:length(segmental_aneuploidy)
		var1 = segmental_aneuploidy(i).chr;
		var2 = segmental_aneuploidy(i).break;
		fprintf(['a2.breakpoint = chr' num2str(var1) ':' num2str(var2) '\n']);
	end;

	%%%% integrate breakpoints.
	segmental_aneuploidy = segmental_aneuploidy_1;
	num_seg_aneu_1       = length(segmental_aneuploidy_1);
	for i = 1:length(segmental_aneuploidy_2)
		segmental_aneuploidy(i+num_seg_aneu_1).chr   = segmental_aneuploidy_2(i).chr;     % chromosome being examined.
		segmental_aneuploidy(i+num_seg_aneu_1).break = segmental_aneuploidy_2(i).break;   % percent along chromosome of edge.
	end;

elseif (exist([workingDir 'Common_ChARM.mat'],'file') ~= 0)
	%%%% Load ChARM:CNV breakpoints.
	dataFile = [workingDir 'Common_ChARM.mat'];
	fprintf(['\nLoading Common_ChARM.mat file for "' projectName '" : ' dataFile '\n']);
	load(dataFile);
elseif (exist([workingDir 'Common_ChARM_SNPratios.mat'],'file') ~= 0)
	%%%% Load ChARM:SNP-ratio breakpoints.
    dataFile = [workingDir 'Common_ChARM_SNPratios.mat'];
    fprintf(['\nLoading Common_ChARM_SNPratios.mat file for "' projectName '" : ' dataFile '\n']);
    load(dataFile);
else
    fprintf(['\nThe ChARM files for "' projectName '" were not found.\n']);
    fprintf(['Analyze your dataset with "analyze_ChARM.sh" first.\n']);
    segmental_aneuploidy = [];
end;


% %% ====================================================================
% % Get segmental aneuploidy data for experiment from "segmental_aneuploidy"
% % This allows for manual annotation of interesting segments of the genome which
% % differ in apparent copy number...  this is not currently in use and is planned
% % to be replaced with an automated analysis function based on the ChARM algorithm.
% %----------------------------------------------------------------------
% if (exist(['links_dir/pileup_dir/' projectName '_segments.txt'],'file') ~= 0)
%     fid = fopen(['links_dir/pileup_dir/' projectName '_segments.txt'],'r');
%     % skip header line.
%     discard = fgetl(fid);
%     clear discard;
% 
%     % initialize some variables.
%     i = 0;
%     lines_analyzed = 0;
%     counter = 0;
%     segmental_aneuploidy = [];
% 
%     % procss each line of the file.
%     while not (feof(fid))
%         i              = i+1;
%         line           = fgetl(fid);
%         lines_analyzed = lines_analyzed+1;
%         
%         % take of interest data fields from each line.
%         segAneu_chr        = sscanf(line, '%s',1);
%         segAneu_break      = sscanf(line, '%s',2);
%         for k = 1:length(sscanf(line,'%s',1));
%             segAneu_break(1) = [];
%         end;
%         
%         % interpret probeID to determine probe chromosome number and location.
%         segmental_aneuploidy(i).chr     = str2double(segAneu_chr);
%         segmental_aneuploidy(i).break   = str2double(segAneu_break);
%     end;
%     fclose(fid);
% else
%     segmental_aneuploidy = [];
% end;

end
