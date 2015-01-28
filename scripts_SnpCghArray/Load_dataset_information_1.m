function [segmental_aneuploidy] = Load_dataset_information_1(projectName,workingDir)
% Load centromere definition file.
%    This is text file containing one header line and two columns.
%    The two columns hold the start and end bp for the centromeres, with
%       respect to each chromosome.

%%=========================================================================
% Load Common_ChARM file for project : 'segmental_aneuploidy'.
%--------------------------------------------------------------------------
if (exist([workingDir 'Common_ChARM.mat'],'file') ~= 0)
    dataFile = [workingDir 'Common_ChARM.mat'];
    fprintf(['\nLoading Common_ChARM file for "' projectName '" : ' dataFile '\n']);
    load(dataFile);
else
    fprintf(['\nThe Common_ChARM file for "' projectName '" was not found.\n']);
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
