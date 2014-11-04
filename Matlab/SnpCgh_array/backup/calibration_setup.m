function calibration_setup();
%
%    This function contains the code used to generate a hapmap file from datasets.   It was run offline and has not yet been adapted for online use.
% There is limited motivation to adapt it for online use, since the SnpCghArray-generated-hapmap has already been calculated and generally it is not
% expected that it will need to be redone.
%
%    This offline-script expects interaction with the user about which datasets to derive calibration data from.
%
%    The following issues would need to be resolved before this could be run.
%    1) Define the following variable.
%        dataset_files : cell list containing file names for datasets for the currently being processed data directory.
%    2) Calibration files in each data set directory refer to the file names instead of the numerical order of the datasets indicated in the header of 'calibration.txt' files.

num_datasets = length(dataset_files);
num_chrs     = 8;

% ask the user which calibration file to use... or to use datasets.
if (calibration_complete == 0)
    v = 0;
    while (v == 0)
        d   = dir;
        str = {d.name};
        j = 1;
        results = [];
        results{1} = '[ Experimental Datasets ]';
        for i = length(str):-1:1
            if (isnan(findstr(str{i}, 'calibration')) == 0)
                j = j+1;
                results{j} = str{i};
            end;
        end;
        [s,v] = listdlg('PromptString','Choose calibration method:',...
            'SelectionMode','single',...
            'ListString',results);
    end;
    
    if (s == 1)
        % load calibration datasets in next code block.
        file_dir = '';
    else
        files_found = 1;
        % get experimental data directory.
        file_dirs = uipickfiles('Prompt','Select data directories to process.','Output','cell');
        for i=1:length(file_dirs)
            if ispc  % Windows
                file_dirs{i} = [file_dirs{i} '\'];
            else     % MacOS
                file_dirs{i} = [file_dirs{i} '/'];
            end;
        end;
        calibration_file = results{s};
        load(calibration_file);
        calibration_complete = 1;
    end;
end;

% load calibration data from experimental datasets.
while (calibration_complete == 0)
    calibration_loaded_from_data = 1;
    % change to directory containing calibration data.
    if ispc  % Windows
        cal_dir = [uigetdir('..','Select SNP/CGH calibration directory.') '\'];
    else     % MacOS
        cal_dir = [uigetdir('..','Select SNP/CGH calibration directory.') '/'];
    end;
    if (exist([cal_dir 'SNP_data.mat'],'file') == 0) || (exist([cal_dir 'CGH_data.mat'],'file') == 0)
        files_found = data_file_load_9(cal_dir);   % saves output as 'CGH_data.mat' & 'SNP_data.mat'.
    end;
    % Load CGH and CNP data to help determine cutoffs
    load([cal_dir 'CGH_data.mat']);
    load([cal_dir 'SNP_data.mat']);
    % define global variables.
    SNP_probeset_length = length(probeset1);
    CGH_probeset_length = length(probeset2);
    bases_per_bin = SNP_bases_per_bin;
    
    % Get flow ploidy values for calibration from "ploidy.txt"
    if (exist([cal_dir 'ploidy.txt'],'file') == 0)
        for dataset = 1:num_datasets
            flow_ploidy{dataset} = 0;
        end;
    else
        fid = fopen([cal_dir 'ploidy.txt'],'r');
        % skip header lines defined as starting with a '#' character.
        firstChar = '#';
        while (firstChar == '#')
            line = fgetl(fid);
            firstChar = line(1);
        end;
        lines_analyzed = 0;
        
        while not (feof(fid))
            if (lines_analyzed > 0)
                line           = fgetl(fid);
            end;
            lines_analyzed = lines_analyzed+1;
            % take of interest data fields from each line.
            exp_number         = sscanf(line, '%s',1);
            exp_ploidy         = sscanf(line, '%s',2);
            for k = 1:length(sscanf(line,'%s',1));
                exp_ploidy(1)  = [];
            end;
            exp_type           = sscanf(line, '%s',3);
            for k = 1:length(sscanf(line,'%s',2));
                exp_type(1)    = [];
            end;
            % place experiment names into used array.
            flow_ploidy{str2num(exp_number)}      = str2double(exp_ploidy);
            flow_ploidy_type{str2num(exp_number)} = exp_type;
        end;
        fclose(fid);
    end;
    
    % Load calibration data from file: "calibration.txt"
    if (exist([cal_dir 'calibration.txt'],'file') ~= 0)
        fprintf(['\nLoading calibration file from: ' strrep(cal_dir,'\','\\') '\n']);
        lines_to_skip = 18;  % number of header lines to be skipped before data analysis.
        fid = fopen([cal_dir 'calibration.txt'],'r');
        
        % skip header lines defined as starting with a '#' character.
        firstChar = '#';
        while (firstChar == '#')
            line = fgetl(fid);
            firstChar = line(1);
        end;
        i = 0;
        lines_analyzed = 0;
        h = waitbar(0.0,'Please wait...','CreateCancelBtn','delete(h); clear h');
        set(h,'Name','Parsing calibration file');
        counter = 0;
        while not (feof(fid))
            counter = counter+1;
            if (counter == 1000)
                waitbar(i/100,h,[num2str(i) '/' num2str(num_lines)]);
                counter = 0;
            end;
            i              = i+1;
            if (lines_analyzed > 0)
                line           = fgetl(fid);
            end;
            lines_analyzed = lines_analyzed+1;
            
            % take of interest data fields from each line.
            cal_chr        = sscanf(line, '%s',1);
            cal_dataset    = sscanf(line, '%s',2);
            for k = 1:length(sscanf(line,'%s',1));
                cal_dataset(1) = [];
            end;
            cal_homologs   = sscanf(line, '%s',3);
            for k = 1:length(sscanf(line,'%s',2));
                cal_homologs(1) = [];
            end;
            
            % interpret probeID to determine probe chromosome number and location.
            calibration(i).chr      = str2double(cal_chr);
            calibration(i).dataset  = cal_dataset;
            calibration(i).homologs = cal_homologs;
        end;
        if exist('h','var'); delete(h); clear h; end;
        fclose(fid);
    end;
    
    %% ====================================================================
    % Get segmental aneuploidy data for experiment from "segmental_aneuploidy.txt"
    %----------------------------------------------------------------------
    if (exist([cal_dir 'segmental_aneuploidy.txt'],'file') ~= 0)
        fprintf(['\nLoading segmental aneuploidy file from: ' strrep(cal_dir,'\','\\')]);
        lines_to_skip = 18;  % number of header lines to be skipped before data analysis.
        fid = fopen([cal_dir 'segmental_aneuploidy.txt'],'r');
        
        % skip header lines defined as starting with a '#' character.
        firstChar = '#';
        while (firstChar == '#')
            line = fgetl(fid);
            firstChar = line(1);
        end;
        i = 0;
        lines_analyzed = 0;
        h = waitbar(0.0,'Please wait...','CreateCancelBtn','delete(h); clear h');
        set(h,'Name','Parsing segmental aneuploidy file');
        counter = 0;
        segmental_aneuploidy = [];
        while not (feof(fid))
            counter = counter+1;
            if (counter == 1000)
                waitbar(i/100,h,[num2str(i) '/' num2str(num_lines)]);
                counter = 0;
            end;
            i              = i+1;
            if (lines_analyzed > 0)
                line           = fgetl(fid);
            end;
            lines_analyzed = lines_analyzed+1;
            
            % take of interest data fields from each line.
            segAneu_dataset    = sscanf(line, '%s',1);
            segAneu_chr        = sscanf(line, '%s',2);
            for k = 1:length(sscanf(line,'%s',1));
                segAneu_chr(1) = [];
            end;
            segAneu_break      = sscanf(line, '%s',3);
            for k = 1:length(sscanf(line,'%s',2));
                segAneu_break(1) = [];
            end;
            
            % interpret probeID to determine probe chromosome number and location.
            segmental_aneuploidy(i).chr     = str2double(segAneu_chr);
            segmental_aneuploidy(i).dataset = segAneu_dataset;
            segmental_aneuploidy(i).break   = str2double(segAneu_break);
        end;
        if exist('h','var'); delete(h); clear h; end;
        fclose(fid);
    end;
    
    for dataset = 1:num_datasets
        h = waitbar(0.0,'Please wait...','CreateCancelBtn','delete(h); clear h');
        set(h,'Name','Plotting CGH data.');
        % Initializes vectors used to hold number of SNPs in each interpretation catagory.
        for chr = 1:num_chrs   % eight chromosomes.
            for j = 1:4   % four SNP interpretation catagories tracked.
                chr_SNPdata{chr,j} = zeros(1,ceil(chr_size(chr)/bases_per_bin));
            end;
            for j = 1:2   % two CGH data catagories tracked.
                chr_CGHdata{chr,j} = zeros(1,ceil(chr_size(chr)/bases_per_bin));
            end;
        end;
        counter = 0;
        for i = 1:CGH_probeset_length
            counter = counter+1;
            if (counter == 500)
                waitbar(i/CGH_probeset_length,h,[num2str(i) '/' num2str(CGH_probeset_length)]);
                counter = 0;
            end;
            % val is the genomic location of the probePair.
            val  = ceil((probeset2(i).probe_location)/bases_per_bin);
            if (strcmp(scale_type,'Ratio') == 1)
                val2 = probeset2(i).probe_Ratio(dataset);
            else
                val2 = probeset2(i).probe_Log2Ratio(dataset);
            end;
            % Determines distribution of CGH data.
            if (isnan(val2) ~= 1)
                % count of data points at CGH locus.
                chr_CGHdata{probeset2(i).probe_chromosome,1}(val) = chr_CGHdata{probeset2(i).probe_chromosome,1}(val)+1;
                % total data at locus.
                chr_CGHdata{probeset2(i).probe_chromosome,2}(val) = chr_CGHdata{probeset2(i).probe_chromosome,2}(val)+val2;
            end;
        end;
        if exist('h','var'); delete(h); clear h; end;
        % divide total Ratio values by number of points per bin.
        for chr = 1:num_chrs
            for j = 1:length(chr_CGHdata{chr,2})
                if (chr_CGHdata{chr,1}(j) == 0)
                    chr_CGHdata{chr,2}(j) = 1;
                else
                    chr_CGHdata{chr,2}(j) = chr_CGHdata{chr,2}(j)/chr_CGHdata{chr,1}(j);
                end;
            end;
        end;
        % precalculation of chromosome copy numbers.
        [chr_breaks{dataset}, chrCopyNum{dataset}] = ...
            FindChrSizes(segmental_aneuploidy,CGH_probeset_length,probeset2,dataset,chr_size,flow_ploidy);
        % precalculation of SNP peaks and cutoffs.
        [realHomozygous_peak{dataset}, disomy_fit{dataset}] = ...
            FindRealHomozygousPeaks_2(chrCopyNum{dataset},SNP_probeset_length,probeset1,dataset,chr_breaks{dataset},chr_size,show_unnassigned,DataTypeToUse,show_fitting);
        [monosomy_peak{dataset},disomy_peak{dataset},trisomy_peak{dataset},tetrasomy_peak{dataset},pentasomy_peak{dataset},hexasomy_peak{dataset} ] = ...
            FindPeaks(realHomozygous_peak{dataset});
        [monosomy_cutoff{dataset},disomy_cutoff{dataset},trisomy_cutoff{dataset},tetrasomy_cutoff{dataset},pentasomy_cutoff{dataset},hexasomy_cutoff{dataset} ] = ...
            FindCutoffs(monosomy_peak{dataset},disomy_peak{dataset},trisomy_peak{dataset},tetrasomy_peak{dataset},pentasomy_peak{dataset},hexasomy_peak{dataset});
    end;
    
    % Reset probe assignments only on the first dataset.
    if (calibration_set == 1)
        SNP_probeset = probeset1;
        h = waitbar(0.0,'Please wait...','CreateCancelBtn','delete(h); clear h');
        set(h,'Name','Resetting SNP probe polarity');
        counter = 0;
        for i = 1:SNP_probeset_length
            counter = counter+1;
            if (counter == 500)
                waitbar(i/SNP_probeset_length,h,[num2str(i) '/' num2str(SNP_probeset_length)]);
                counter = 0;
            end;
            %unassign probes.
            SNP_probeset(i).probe_polarity   = 5;
            SNP_probeset(i).assigned_correct = 0;  % correct_assignment
            SNP_probeset(i).assigned_wrong   = 0;  % incorrect_assignment.
        end;
        if exist('h','var'); delete(h); clear h; end;
    end;
    
    %--------------------------------------------------------------------------
    % Define SNP probe pair polarity.
    %--------------------------------------------------------------------------
    fprintf('Assigning SNPS to homologs based on calibration data.\n');
    h = waitbar(0.0,'Please wait...','CreateCancelBtn','delete(h); clear h');
    set(h,'Name','Assigning SNP probe pairs');
    counter = 0;
    for i = 1:2:SNP_probeset_length
        counter = counter+1;
        if (counter == 500)
            waitbar(i/SNP_probeset_length,h,[num2str(i) '/' num2str(SNP_probeset_length)]);
            counter = 0;
        end;
        % determines if probe pair is useful; both probes have data.
        %if (isnan(probeset1(i).probe_Ratio(dataset)) == 0) && (isnan(probeset1(i+1).probe_Ratio(dataset)) == 0)
        % Assigns probe pairs to specific homologs, using calibration strain data.
        %    Calibration strains from (Forche 2008 [PMID: 18162019])
        for calibrator = 1:length(calibration)
            if (probeset1(i).probe_chromosome == calibration(calibrator).chr)
		dataset_file = calibration(calibrator).dataset;
		for ii = 1:length(dataset_files)
		    if (strcmp(dataset_file,dataset_files{ii}) == 1)
			dataset = ii;
		    end;
		end;
                if (isnan(probeset1(i).probe_Ratio(dataset)) == 0) && (isnan(probeset1(i+1).probe_Ratio(dataset)) == 0)
                    % both probes are usable in this calibration
                    % dataset...  we can use this to inform SNP
                    % phasing.
                    if (DataTypeToUse == 1)   % (1)AllelicFraction; (2)Angle.
                        dataPoint = probeset1(i).probe_Ratio(dataset)/(probeset1(i).probe_Ratio(dataset)+probeset1(i+1).probe_Ratio(dataset));
                    elseif (DataTypeToUse == 2)
                        dataPoint = atan2(probeset1(i).probe_Ratio(dataset),probeset1(i+1).probe_Ratio(dataset));
                    end;
                    if (strcmp(calibration(calibrator).homologs,'b/b') == 1)
                        if (dataPoint < disomy_cutoff{dataset}(1))
                            SNP_probeset(i  ).assigned_correct    = SNP_probeset(i  ).assigned_correct    + 1; % correct_assignment += 1;
                        elseif (dataPoint > disomy_cutoff{dataset}(2))
                            SNP_probeset(i  ).assigned_wrong = SNP_probeset(i  ).assigned_wrong + 1; % incorrect_assignment += 1;
                        end;
                    elseif (strcmp(calibration(calibrator).homologs,'a/b/b') == 1)
                        if (dataPoint < trisomy_cutoff{dataset}(2))
                            SNP_probeset(i  ).assigned_correct    = SNP_probeset(i  ).assigned_correct    + 1; % correct_assignment += 1;
                        elseif (dataPoint > trisomy_cutoff{dataset}(2))
                            SNP_probeset(i  ).assigned_wrong = SNP_probeset(i  ).assigned_wrong + 1; % incorrect_assignment += 1;
                        end;
                    elseif (strcmp(calibration(calibrator).homologs,'a/a') == 1)
                        if (dataPoint < disomy_cutoff{dataset}(1))
                            SNP_probeset(i  ).assigned_wrong = SNP_probeset(i  ).assigned_wrong + 1; % incorrect_assignment += 1;
                        elseif (dataPoint > disomy_cutoff{dataset}(2))
                            SNP_probeset(i  ).assigned_correct    = SNP_probeset(i  ).assigned_correct    + 1; % correct_assignment += 1;
                        end;
                    elseif (strcmp(calibration(calibrator).homologs,'a/a/b') == 1)
                        if (dataPoint < trisomy_cutoff{dataset}(2))
                            SNP_probeset(i  ).assigned_wrong = SNP_probeset(i  ).assigned_wrong + 1; % incorrect_assignment += 1;
                        elseif (dataPoint > trisomy_cutoff{dataset}(2))
                            SNP_probeset(i  ).assigned_correct    = SNP_probeset(i  ).assigned_correct    + 1; % correct_assignment += 1;
                        end;
                    end;
                end;
            end;
        end;
    end;
    if exist('h','var'); delete(h); clear h; end;
    
    choice1 = questdlg('Would you like to use additional calibration data?','Calibrate more', 'Yes','No','No');
    % Handle response
    switch choice1
        case 'Yes'
            calibration_complete = 0;
            calibration_set      = calibration_set + 1;
        case 'No'
            calibration_complete = 1;
    end;
    clear line;
end;

%% Save calibration data.
if (calibration_complete == 1) && (calibration_loaded_from_data == 1)
    choice1 = questdlg('Would you like to save calibration data?','Save', 'Yes','No','Yes');
    switch choice1
        case 'Yes'
            prompt = {'enter name.   [*.mat]'};
            dlg_title = 'Save Calibration Data';
            num_lines = 1;
            def = {'calibration','hsv'};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            save([char(answer) '.mat'], 'SNP_probeset');
    end;
end;
