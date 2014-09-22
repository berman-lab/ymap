%% ========================================================================
% Generate graphs involved in analysis of SNP probe data from 1' Agilent
% slide.
%==========================================================================
clear;

% Attempting to improve cutoff calculations:
% 1) Cutoffs at 0.4 and 0.6 allelic fraction.
% 2) Cutoffs at 50% between adjacent peaks.
% 3) Cutoffs at equal probability marks between adjacent peaks, using
% Gaussian fittings.
% 4) Implement a skewed Gaussian fitting equation.
% 5) Update Gaussian cutoffs to account for skew Gaussian.
%
% Displays SNP and CGH data in a common figure format.
%    CGH_Genomic_display        : Draws CGH data in figure.
%    colorBars                  : Draws SNP information as background colors on chromosomes.
%    blendColorBars             : Blends adjacent colorBars to remove high contrast edges.
%    infillColorBars            : Fills in "no-data" areas with adjacent colored areas.
%    AnglePlot                  : Produces angle histogram of scatterplots, to left of chromosome.
%       FillColors              : Fills angle histogram with homolog identity colors.
%    HistPlot                   : Produces histogram of CGH data, to right of chromosomes.
%    ChrNum                     : Draws a large numeral of copy# estimage to right of chromosome.
%    Yscale_nearest_even_ploidy : automatically alters y-axis from [0..4] to [0..8] if ploidy is ~tetraploid.
%    Show_Genomic_LOH_fraction  : Shows an automatically determined %LOH, from how many SNP probes are homozygous.
%    show_unnassigned           : Shows unnassigned SNP pair data in histogram/chromosome plots.
%    Centromere_format          : Controls how centromeres are depicted.   [0..2]   '2' is pinched cartoon default.
%    SNP_bases_per_bin          : Controls bin sizes for SNP/CGH fractions of plot.
%    scale_type                 : 'Ratio' or 'Log2Ratio' y-axis scaling of copy number.
%                                 'Log2Ratio' does not properly scale CGH data by ploidy.
%    Chr_max_width              : max width of chromosomes as fraction of figure width.
SNP_Genomic_display              = true;
    CGH_Genomic_display          = true;
    colorBars                    = true;
    blendColorBars               = false;
    infillColorBars              = false;
    AnglePlot                    = true;
        FillColors               = true;
    HistPlot                     = true;
    ChrNum                       = true;
    Yscale_nearest_even_ploidy   = true;
    Show_Genomic_LOH_fraction    = true;
    show_unnassigned             = false;
    show_Xlabels                 = false;
    show_MRS                     = true;
    Centromere_format            = 2;
    SNP_bases_per_bin            = 5000;
    scale_type                   = 'Ratio';
    Chr_max_width                = 0.8;

    %not in place yet... uninformative currently looks like unnassigned to script.
    show_uninformative           = false;

% Makes figures showing Gaussian fits for each chromosome segment.
Gaussian_fit_display             = false;
    DataTypeToUse                = 1;   % (1)AllelicFraction; (2)Angle.
    show_fitting                 = 0;   % (0)false; (~0)figure number to use.
    
% Analyze incidence of SNP interpretation runs.
SNP_Runs_analysis                = false;

% Displays ratio scatterplot: sub-options collect all data per set or per chr.
Scatter_display                  = false;    % make scatter plots.
    main_scatter                 = true;    % actually make main scatter plot.
    per_chr_scatter              = true;    % for each chromosome instead of for entire dataset.
    distanceCutoff               = false;   % when making angleplots, ignores pairs with magnitude less than half of peak.
    % Alternate scatter plots used in literature.   (All are equivalent.)
    AllelicFraction_display      = false;   % Forche-2005.
    LogFraction_display          = false;   % O'Meara-2002, Lovmar-2003, Fan-2000.
    Intensity_v_ratio_display    = false;   % Pastinen-2000.

% Scatter plot with probe length indicated by spot radius.
Scatter_probeLength_display      = false;

% Displays distribution of probes across genome.
Display_distribution             = false;
logScale                         = false;   % 'false' means a linear scale is used.

%%=========================================================================
% Control variables for unfinished figure methods.
%--------------------------------------------------------------------------

% Attempt to simulate SC5314 used as reference and experimental strain.
% Results are poor due to lack of proper normalization using this method.
SC5314_display                   = false;

% Attempt to allow user to zoom into a chromosome for close-up examination
% of SNP/CGH data.
Zoom_view                        = false;       % should a close-up look be done?
Zoom_chr                         = 1;           % which chromosome to look closely at.
Zoom_range                       = [0.0 100.0]; % percent start and stop.

%%=========================================================================
% Control variables for Candida albicans SC5314.
%--------------------------------------------------------------------------
% Determines for which chromosomes to construct figures from.
chromosomes_to_analyze = 1:8;

% Defines chromosome sizes in bp. (diploid total=28,567,7888)
% Defines centromere locations in bp.
% Defines MRS locations in bp.
[centromeres, chr_sizes,MRSs] = Load_information();
for i = 1:8
    chr_size(i)  = chr_sizes(i).size;
    cen_start(i) = centromeres(i).start;
    cen_end(i)   = centromeres(i).end;
end;
for i = 1:length(MRSs)
    MRS_chr(i)   = MRSs(i).chr;
    MRS_location(i) = (MRSs(i).start+MRSs(i).end)/2;
end;
clear centromeres chr_sizes MRSs;

%%=========================================================================
%%=========================================================================
%%= No further control variables below. ===================================
%%=========================================================================
%%=========================================================================

%% ========================================================================
% Define SNP probe pair polarity.
%==========================================================================
choice1 = questdlg('Would you like to use calibration data?','Calibrate', 'Yes','No','Yes');

% Handle response
switch choice1
    case 'Yes'
        calibration_complete = 0;
        calibration_set      = 1;
        no_calibration       = 0;
        show_uncalibrated    = false;
        file_dirs            = [];
    case 'No'
        show_uncalibrated    = true;
        calibration_complete = 1;
        calibration_set      = 0;
        no_calibration       = 1;
        files_found          = 1;
        % get experimental data directory.
        file_dirs = uipickfiles('Prompt','Select data directories to process.','Output','cell');
        for i=1:length(file_dirs)
            if ispc  % Windows
                file_dirs{i} = [file_dirs{i} '\'];
            else     % MacOS
                file_dirs{i} = [file_dirs{i} '/'];
            end;
        end;
end;
calibration_loaded_from_data = 0;

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
        files_found         = data_file_load_6(cal_dir);
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
        for i = 1:8
            flow_ploidy{i} = 0;
        end;
    else
        fid = fopen([cal_dir 'ploidy.txt'],'r');
        while not (feof(fid))
            line           = fgetl(fid);
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
        %skip header lines defined by 'lines_to_skip'.
        if (lines_to_skip > 0)
            for j = 1:lines_to_skip
                discard = fgetl(fid);
            end;
            clear discard;
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
            line           = fgetl(fid);
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
            calibration(i).dataset  = str2double(cal_dataset);
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
        %skip header lines defined by 'lines_to_skip'.
        if (lines_to_skip > 0)
            for j = 1:lines_to_skip
                discard = fgetl(fid);
            end;
            clear discard;
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
            line           = fgetl(fid);
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
            segmental_aneuploidy(i).dataset = str2double(segAneu_dataset);
            segmental_aneuploidy(i).break   = str2double(segAneu_break);
        end;
        if exist('h','var'); delete(h); clear h; end;
        fclose(fid);
    end;
    
    % sequentially process each dataset, adding to SNP_probeset as we go.
    for dataset = 1:8
        h = waitbar(0.0,'Please wait...','CreateCancelBtn','delete(h); clear h');
        set(h,'Name','Plotting CGH data.');
        % Initializes vectors used to hold number of SNPs in each interpretation catagory.
        for i = 1:8   % eight chromosomes.
            for j = 1:4   % four SNP interpretation catagories tracked.
                chr_SNPdata{i,j} = zeros(1,ceil(chr_size(i)/bases_per_bin));
            end;
            for j = 1:2   % two CGH data catagories tracked.
                chr_CGHdata{i,j} = zeros(1,ceil(chr_size(i)/bases_per_bin));
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
        for i = 1:8
            for j = 1:length(chr_CGHdata{i,2})
                if (chr_CGHdata{i,1}(j) == 0)
                    chr_CGHdata{i,2}(j) = 1;
                else
                    chr_CGHdata{i,2}(j) = chr_CGHdata{i,2}(j)/chr_CGHdata{i,1}(j);
                end;
            end;
        end;
        % precalculation of chromosome copy numbers. dddd
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
            SNP_probeset(i).probe_polarity = 0;
            SNP_probeset(i).assign_cyan    = 0;  % correct_assignment
            SNP_probeset(i).assign_magenta = 0;  % incorrect_assignment.
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
                for dataset = calibration(calibrator).dataset
                    if (isnan(probeset1(i).probe_Ratio(dataset)) == 0) && (isnan(probeset1(i+1).probe_Ratio(dataset)) == 0)
                        if (DataTypeToUse == 1)   % (1)AllelicFraction; (2)Angle.
                            dataPoint = probeset1(i).probe_Ratio(dataset)/(probeset1(i).probe_Ratio(dataset)+probeset1(i+1).probe_Ratio(dataset));
                        elseif (DataTypeToUse == 2)
                            dataPoint = atan2(probeset1(i).probe_Ratio(dataset),probeset1(i+1).probe_Ratio(dataset));
                        end;
                        if (strcmp(calibration(calibrator).homologs,'b/b') == 1)
                            if (dataPoint < disomy_cutoff{dataset}(1))
                                SNP_probeset(i  ).assign_cyan    = SNP_probeset(i  ).assign_cyan    + 1; % correct_assignment += 1;
                            elseif (dataPoint > disomy_cutoff{dataset}(2))
                                SNP_probeset(i  ).assign_magenta = SNP_probeset(i  ).assign_magenta + 1; % incorrect_assignment += 1;
                            end;
                        elseif (strcmp(calibration(calibrator).homologs,'a/b/b') == 1)
                            if (dataPoint < trisomy_cutoff{dataset}(2))
                                SNP_probeset(i  ).assign_cyan    = SNP_probeset(i  ).assign_cyan    + 1; % correct_assignment += 1;
                            elseif (dataPoint > trisomy_cutoff{dataset}(2))
                                SNP_probeset(i  ).assign_magenta = SNP_probeset(i  ).assign_magenta + 1; % incorrect_assignment += 1;
                            end;
                        elseif (strcmp(calibration(calibrator).homologs,'a/a') == 1)
                            if (dataPoint < disomy_cutoff{dataset}(1))
                                SNP_probeset(i  ).assign_magenta = SNP_probeset(i  ).assign_magenta + 1; % incorrect_assignment += 1;
                            elseif (dataPoint > disomy_cutoff{dataset}(2))
                                SNP_probeset(i  ).assign_cyan    = SNP_probeset(i  ).assign_cyan    + 1; % correct_assignment += 1;
                            end;
                        elseif (strcmp(calibration(calibrator).homologs,'a/a/b') == 1)
                            if (dataPoint < trisomy_cutoff{dataset}(2))
                                SNP_probeset(i  ).assign_magenta = SNP_probeset(i  ).assign_magenta + 1; % incorrect_assignment += 1;
                            elseif (dataPoint > trisomy_cutoff{dataset}(2))
                                SNP_probeset(i  ).assign_cyan    = SNP_probeset(i  ).assign_cyan    + 1; % correct_assignment += 1;
                            end;
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

if isempty(file_dirs)
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
end;

%% ========================================================================
% process each data directory, sequentially.
%--------------------------------------------------------------------------
for file_dir_num = 1:length(file_dirs)
    file_dir = file_dirs{file_dir_num};

    %% ====================================================================
    % Setup for experimental data analysis using calibration data.
    % Swap experimental SNP data per hapmap/calibration file.
    %----------------------------------------------------------------------
    %info from design considerations.
    load SNP_probeset_2.mat;
    SNP_probeset_length = length(probeset_2);
    
    % Assign SNP probes to homologs, based on calibration data collected earlier.
    if (no_calibration == 0)
        for i = 1:2:SNP_probeset_length
            if (strcmp(probeset_2(i).probe_sequence,probeset_2(i+1).probe_sequence) == 1)
                % Gresham-designed SNP probes are identical.
                % There was a bug in Gresham's script which did not force the
                % SNP to be in the center of the probe pair.
                % 1197 probe pairs are garbage because of this.
                SNP_probeset(i  ).probe_polarity = 4;
                SNP_probeset(i+1).probe_polarity = 4;
            else
                if (SNP_probeset(i).assign_cyan > SNP_probeset(i).assign_magenta) % if (correct_assignment > incorrect_assignment)
                    SNP_probeset(i  ).probe_polarity = 1;
                    SNP_probeset(i+1).probe_polarity = 1;
                elseif (SNP_probeset(i).assign_cyan < SNP_probeset(i).assign_magenta) % if (correct_assignment < incorrect_assignment)
                    SNP_probeset(i  ).probe_polarity = 2;
                    SNP_probeset(i+1).probe_polarity = 2;
                elseif (SNP_probeset(i).assign_cyan == SNP_probeset(i).assign_magenta) % if (correct_assignment = incorrect_assignment)
                    SNP_probeset(i  ).probe_polarity = 0;
                    SNP_probeset(i+1).probe_polarity = 0;
                end;
            end;
        end;
    end;
    
    % load experimental values for all the probes.
    % Analyze results files if this has not been done.
    if (exist([file_dir 'SNP_data.mat'],'file') == 0) || (exist([file_dir 'CGH_data.mat'],'file') == 0)
        files_found = data_file_load_6(file_dir);
    else
        files_found = 1;
    end;
    if (files_found == 0)
        fprintf('\nNo results files found, perhaps they are in a different location.\n');
    end;
    
    load([file_dir 'SNP_data.mat']);
    load([file_dir 'CGH_data.mat']);
    
    SNP_probeset_length = length(probeset1);
    CGH_probeset_length = length(probeset2);
    
    %% ====================================================================
    % Get experiment names from "names.txt"
    %----------------------------------------------------------------------
    if (exist([file_dir 'names.txt'],'file') == 0)
        for i = 1:8
            names{i} = ['Experiment ' num2str(i)];   %default expected size;
        end;
    else
        fprintf(['\nLoading experiment names file from: ' strrep(file_dir,'\','\\')]);
        fid = fopen([file_dir 'names.txt'],'r');
        while not (feof(fid))
            line           = fgetl(fid);
            % take of interest data fields from each line.
            exp_number         = sscanf(line, '%s',1);
            exp_name           = line;
            for k = 1:length(sscanf(line,'%s',1));
                exp_name(1) = [];
            end;
            exp_name(1) = [];
            
            % place experiment names into used array.
            names{str2num(exp_number)} = exp_name;
        end;
        fclose(fid);
    end;
    
    %% ====================================================================
    % Get flow ploidy values for experiment from "ploidy.txt"
    %----------------------------------------------------------------------
    if (exist([file_dir 'ploidy.txt'],'file') == 0)
        for i = 1:8
            flow_ploidy{i} = 0;
        end;
    else
        fprintf(['\nLoading ploidy estimates file from: ' strrep(file_dir,'\','\\')]);
        fid = fopen([file_dir 'ploidy.txt'],'r');
        while not (feof(fid))
            line           = fgetl(fid);
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
    
    %% ====================================================================
    % Get segmental aneuploidy data for experiment from "segmental_aneuploidy.txt"
    %----------------------------------------------------------------------
    if (exist([file_dir 'segmental_aneuploidy.txt'],'file') ~= 0)
        fprintf(['\nLoading segmental aneuploidy file from: ' strrep(file_dir,'\','\\')]);
        lines_to_skip = 18;  % number of header lines to be skipped before data analysis.
        fid = fopen([file_dir 'segmental_aneuploidy.txt'],'r');
        %skip header lines defined by 'lines_to_skip'.
        if (lines_to_skip > 0)
            for j = 1:lines_to_skip
                discard = fgetl(fid);
            end;
            clear discard;
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
            line           = fgetl(fid);
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
            segmental_aneuploidy(i).dataset = str2double(segAneu_dataset);
            segmental_aneuploidy(i).break   = str2double(segAneu_break);
        end;
        if exist('h','var'); delete(h); clear h; end;
        fclose(fid);
    end;
    
    %% ====================================================================
    % Apply probe polarity assignment from calibration data to experimental data.
    %----------------------------------------------------------------------
    fprintf(['\nPhasing experimental SNP data from: ' strrep(file_dir,'\','\\')]);
    h = waitbar(0.0,'Please wait...','CreateCancelBtn','delete(h); clear h');
    set(h,'Name','Assigning SNP probe pairs');
    counter = 0;
    for i = 1:SNP_probeset_length
        counter = counter+1;
        if (counter == 500)
            waitbar(i/SNP_probeset_length,h,[num2str(i) '/' num2str(SNP_probeset_length)]);
            counter = 0;
        end;
        if (no_calibration == 1)
            probeset1(i).probe_polarity = 0;
        else
            probeset1(i).probe_polarity = SNP_probeset(i).probe_polarity;
        end;
    end;
    if exist('h','var'); delete(h); clear h; end;

    %% ====================================================================
    % Determine chromosome copy numbers.
    %----------------------------------------------------------------------
    datasetDetails = [];
    fprintf('\nDetermining chromsome copy numbers for microarray: ');
    for dataset = 1:length(names)
        fprintf([num2str(dataset) ' ']);
        [chr_breaks, chrCopyNum] = ...
            FindChrSizes(segmental_aneuploidy,CGH_probeset_length,probeset2,dataset,chr_size,flow_ploidy);
        datasetDetails{dataset}.chr_breaks = chr_breaks;
        datasetDetails{dataset}.chrCopyNum = chrCopyNum;
    end;
    
    %% ====================================================================
    % Determine cutoffs for experimental dataset chromosome segments
    %----------------------------------------------------------------------    
    fprintf('\nDetermining SNP interpretation cutoffs for microarray: ');
    for dataset = 1:length(names)
        fprintf([num2str(dataset) ' ']);
        % Finds initial homozygous peak locations.
        [realHomozygous_peak, disomy_fit,skew_factor] = ...
            FindRealHomozygousPeaks_2(chrCopyNum,SNP_probeset_length,probeset1,dataset,chr_breaks,chr_size,show_unnassigned,DataTypeToUse,show_fitting);
        % Finds real peak locations.
        [monosomy_peak,disomy_peak,trisomy_peak,tetrasomy_peak,pentasomy_peak,hexasomy_peak ] = ...
            FindPeaks(realHomozygous_peak);
        % Determine cutoffs between peaks for each datasets:chromosome:segment.
        for chromosome = chromosomes_to_analyze;
            for segment = 1:length(chrCopyNum{chromosome})
                chrCopyNum = datasetDetails{dataset}.chrCopyNum;
                chr_breaks = datasetDetails{dataset}.chr_breaks;
                % Determins cutoffs for a single chromosome segment, using intersections of Gaussian fit curves.
                MakeFigure = Gaussian_fit_display;
                [raw,smoothed,peaks,actual_cutoffs,mostLikelyGaussians,chrCopyNum] = ...
                    FindGaussianCutoffs_2(probeset1,chrCopyNum,chr_breaks,chr_size,dataset,chromosome,...
                    segment,monosomy_peak,disomy_peak,trisomy_peak,tetrasomy_peak,pentasomy_peak,...
                    hexasomy_peak,skew_factor,names{dataset},file_dir,Gaussian_fit_display,show_fitting,...
                    DataTypeToUse);

                datasetDetails{dataset}.histogram_raw{chromosome,segment}      = raw;
                datasetDetails{dataset}.histogram_smooth{chromosome,segment}   = smoothed;
                datasetDetails{dataset}.peaks{chromosome,segment}              = peaks/200;
                datasetDetails{dataset}.cutoffs{chromosome,segment}            = actual_cutoffs/200;
                datasetDetails{dataset}.importantGaussians{chromosome,segment} = mostLikelyGaussians;
                datasetDetails{dataset}.chrCopyNum                             = chrCopyNum;
            end;
        end;
    end;
    
    %% ====================================================================
    % Save datasetDetails.
    %----------------------------------------------------------------------    
    save([file_dir 'datasetDetails.mat'], 'datasetDetails');
    
    %% ========================================================================
    % Plotting probes across genome by interpretation catagory after accounting
    % for polarity of SNP pairs.
    %    CGH info plotted as histogram.
    %    SNP info plotted as colorbar.
    %==========================================================================
    if (SNP_Genomic_display == 1)
        fprintf(['\nGenerating SNP/CGH figures from: ' strrep(file_dir,'\','\\')]);
        % define number of and labels for chromosomes.
        chr_labels = {'Chr1','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7','ChrR'};
        num_chr = length(chr_labels);
        
        bases_per_bin = SNP_bases_per_bin;
        chr_length_scale_multiplier = 1/bases_per_bin;
        maxY = 10;
        for dataset = 1:length(names);
            % Initializes vectors used to hold number of SNPs in each
            %    interpretation catagory.
            for chromosome = 1:8   % eight chromosomes.
                for j = 1:14   % 14 SNP interpretation catagories tracked.
                    chr_SNPdata{chromosome,j} = zeros(1,ceil(chr_size(chromosome)/bases_per_bin));
                end;
                for j = 1:2   % two CGH data catagories tracked.
                    chr_CGHdata{chromosome,j} = zeros(1,ceil(chr_size(chromosome)/bases_per_bin));
                end;
            end;

            % tally CGH data into bins for main figure histogram.
            h = waitbar(0.0,'Please wait...','CreateCancelBtn','delete(h); clear h');
            set(h,'Name','Plotting CGH data.');
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
            for i = 1:8
                for j = 1:length(chr_CGHdata{i,2})
                    if (chr_CGHdata{i,1}(j) == 0)
                        chr_CGHdata{i,2}(j) = 1;
                    else
                        chr_CGHdata{i,2}(j) = chr_CGHdata{i,2}(j)/chr_CGHdata{i,1}(j);
                    end;
                end;
            end;

            h = waitbar(0.0,'Please wait...','CreateCancelBtn','delete(h); clear h');
            set(h,'Name','Plotting SNP data.');
            counter = 0;
            SNPs_hom = 0;
            SNPs_total = 0;
            for i = 1:2:SNP_probeset_length
                % probeset1 consists of sequential pairs, representing each allele for each SNP.
                %    The odd numbered probes are the first of each pair.
                counter = counter+1;
                if (counter == 500)
                    waitbar(i/SNP_probeset_length,h,[num2str(i) '/' num2str(SNP_probeset_length)]);
                    counter = 0;
                end;
                
                % determines if probe pair is useful; both probes have data.
                if (isnan(probeset1(i).probe_Ratio(dataset)) == 0) && (isnan(probeset1(i+1).probe_Ratio(dataset)) == 0)
                    % Determines distribution of SNP pairs in four catagories.
                    %  0:'ab'
                    %  1:'a'/'aa'/'aaa'/'aaaa'/'aaaaa'/'aaaaaa'
                    %  2:'b'/'bb'/'bbb'/'bbbb'/'bbbbb'/'bbbbbb'
                    %  3:'aab'
                    %  4:'abb'
                    %  5:'aaab'
                    %  6:'abbb'
                    %  7:'aaaab'
                    %  8:'aaabb'
                    %  9:'aabbb'
                    % 10:'abbbb'
                    % 11:'aaaaab'
                    % 12:'abbbbb'
                    % 13: hom unnassigned.
                    
                    % collect SNP data into interpretation arrays.
                    for iii = 1 %ddddd
                        % Find chromosome of probe.
                        chromosome = probeset1(i).probe_chromosome;
                        
                        % Find which segment the probe is found on.
                        if (length(datasetDetails{dataset}.chr_breaks{chromosome}) > 2)
                            location = probeset1(i).probe_location/chr_size(chromosome);
                            for last_break = length(datasetDetails{dataset}.chr_breaks{chromosome}):-1:2
                                if (location < datasetDetails{dataset}.chr_breaks{chromosome}(last_break))
                                    segment = last_break-1;
                                end;
                            end;
                        else
                            segment = 1;
                        end;
                        
                        % Calculate value of SNP probe pair.
                        [UsedData] = calculateValue(probeset1,i,DataTypeToUse,dataset);

                        % Determine which section of the allelic fraction histogram the probe pair falls into.
                        done = false;
                        section = datasetDetails{dataset}.importantGaussians{chromosome,segment}(length(datasetDetails{dataset}.importantGaussians{chromosome,segment}));
                        for cutoff = 1:length(datasetDetails{dataset}.cutoffs{chromosome,segment})
                            if (done == false) && (UsedData <= datasetDetails{dataset}.cutoffs{chromosome,segment}(cutoff))
                                done = true;
                                section = datasetDetails{dataset}.importantGaussians{chromosome,segment}(cutoff);
                            end;
                        end;

                        % Determines probe assignment of SNP pair.
                        ChrCopyNumber = datasetDetails{dataset}.chrCopyNum{probeset1(i).probe_chromosome}(segment);
                        if (probeset1(i).probe_polarity == 0)
                            probeset1(i).probe_assignment   = 13;
                            probeset1(i+1).probe_assignment = 13;
                            if (show_unnassigned == true) && (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
                        elseif (probeset1(i).probe_polarity == 4)
                            % null action for when (probe_polarity == 4) due to probe design error; probes are identical.
                        elseif (ChrCopyNumber <= 1)     %monosomy
                            % Assigns probe pairs to specific interpretations, collects data for each interpretation.
                            if (section == 1)   % monosomy:'a'
                                probeset1(i).probe_assignment   = 1;
                                probeset1(i+1).probe_assignment = 1;
                                if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
                            else % (section == 2)   % monosomy:'b'
                                probeset1(i).probe_assignment   = 2;
                                probeset1(i+1).probe_assignment = 2;
                                if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
                            end;
                        elseif (ChrCopyNumber <= 2) %disomy
                            % Assigns probe pairs to specific interpretations, collects data for each interpretation.
                            if (section == 1)   % disomy:'aa'
                                probeset1(i).probe_assignment   = 1;
                                probeset1(i+1).probe_assignment = 1;
                                if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
                            elseif (section == 3)   % disomy:'bb'
                                probeset1(i).probe_assignment   = 2;
                                probeset1(i+1).probe_assignment = 2;
                                if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
                            else % (section == 2)   % disomy:'ab'
                                probeset1(i).probe_assignment   = 0;
                                probeset1(i+1).probe_assignment = 0;
                                if (Show_Genomic_LOH_fraction == true);   SNPs_total = SNPs_total+1;   end;
                            end;
                        elseif (ChrCopyNumber <= 3) %trisomy
                            % Assigns probe pairs to specific interpretations, collects data for each interpretation.
                            if (section == 1)   % trisomy:'aaa'
                                probeset1(i).probe_assignment   = 1;
                                probeset1(i+1).probe_assignment = 1;
                                if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
                            elseif (section == 4)   % trisomy:3 'bbb'
                                probeset1(i).probe_assignment   = 2;
                                probeset1(i+1).probe_assignment = 2;
                                if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
                            elseif (section == 2)   % trisomy:'aab'
                                probeset1(i).probe_assignment   = 3;
                                probeset1(i+1).probe_assignment = 3;
                                if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
                            else % (section == 3)   % trisomy:'abb'
                                probeset1(i).probe_assignment   = 4;
                                probeset1(i+1).probe_assignment = 4;
                                if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
                            end;
                        elseif (ChrCopyNumber <= 4) %tetrasomy
                            % Assigns probe pairs to specific interpretations, collects data for each interpretation.
                            if (section == 1)   % tetrasomy:'aaaa'
                                probeset1(i).probe_assignment   = 1;
                                probeset1(i+1).probe_assignment = 1;
                                if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
                            elseif (section == 5)   % tetrasomy:'bbbb'
                                probeset1(i).probe_assignment   = 2;
                                probeset1(i+1).probe_assignment = 2;
                                if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
                            elseif (section == 2)   % tetrasomy:'aaab'
                                probeset1(i).probe_assignment   = 5;
                                probeset1(i+1).probe_assignment = 5;
                                if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
                            elseif (section == 4)   % tetrasomy:'abbb'
                                probeset1(i).probe_assignment   = 6;
                                probeset1(i+1).probe_assignment = 6;
                                if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
                            else % (section == 3)   % tetrasomy:'aabb'
                                probeset1(i).probe_assignment   = 0;
                                probeset1(i+1).probe_assignment = 0;
                                if (Show_Genomic_LOH_fraction == true);   SNPs_total = SNPs_total+1;   end;
                            end;
                        elseif (ChrCopyNumber <= 5) %pentasomy
                            % Assigns probe pairs to specific interpretations, collects data for each interpretation.
                            if (section == 1)   % pentasomy:'aaaaa'
                                probeset1(i).probe_assignment   = 1;
                                probeset1(i+1).probe_assignment = 1;
                                if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
                            elseif (section == 6)   % pentasomy:'bbbbb'
                                probeset1(i).probe_assignment   = 2;
                                probeset1(i+1).probe_assignment = 2;
                                if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
                            elseif (section == 2)   % pentasomy:'aaaab'
                                probeset1(i).probe_assignment   = 7;
                                probeset1(i+1).probe_assignment = 7;
                                if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
                            elseif (section == 5)   % pentasomy:'abbbb'
                                probeset1(i).probe_assignment   = 10;
                                probeset1(i+1).probe_assignment = 10;
                                if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
                            elseif (section == 3)   % pentasomy:'aaabb'
                                probeset1(i).probe_assignment   = 8;
                                probeset1(i+1).probe_assignment = 8;
                                if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
                            else % (section == 4)   % pentasomy:'aabbb'
                                probeset1(i).probe_assignment   = 9;
                                probeset1(i+1).probe_assignment = 9;
                                if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
                            end;
                        else %treat as hexasomy: if (ChrCopyNumber <= 6) %hexasomy
                            % Assigns probe pairs to specific interpretations, collects data for each interpretation.
                            if (section == 1)   % hexasomy:'aaaaaa'
                                probeset1(i).probe_assignment   = 1;
                                probeset1(i+1).probe_assignment = 1;
                                if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
                            elseif (section == 7)   % hexasomy:'bbbbbb'
                                probeset1(i).probe_assignment   = 2;
                                probeset1(i+1).probe_assignment = 2;
                                if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
                            elseif (section == 2)   % hexasomy:'aaaaab'
                                probeset1(i).probe_assignment   = 11;
                                probeset1(i+1).probe_assignment = 11;
                                if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
                            elseif (section == 6)   % hexasomy:'abbbbb'
                                probeset1(i).probe_assignment   = 12;
                                probeset1(i+1).probe_assignment = 12;
                                if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
                            elseif (section == 3)   % hexasomy:'aaaabb'
                                probeset1(i).probe_assignment   = 3;
                                probeset1(i+1).probe_assignment = 3;
                                if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
                            elseif (section == 5)   % hexasomy:'aabbbb'
                                probeset1(i).probe_assignment   = 4;
                                probeset1(i+1).probe_assignment = 4;
                                if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
                            else % (section == 4)   % hexasomy:'aaabbb'
                                probeset1(i).probe_assignment   = 0;
                                probeset1(i+1).probe_assignment = 0;
                                if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
                            end;
                        end;
                    end;
                    
                    % val is the genomic location of the probePair.
                    val = ceil(probeset1(i).probe_location/bases_per_bin);

                    % increment the appropriate interpretation category.
                    if (length(probeset1(i).probe_assignment) > 0)
                        chr_SNPdata{probeset1(i).probe_chromosome,probeset1(i).probe_assignment+1}(val) = ...
                            chr_SNPdata{probeset1(i).probe_chromosome,probeset1(i).probe_assignment+1}(val)+1;
                    end;
                end;
            end;
            if exist('h','var'); delete(h); clear h; end;
            
            % basic plot parameters.
            left            = 0.15;
            height          = 0.5/num_chr;
            base            = 0.1;
            vertical_margin = 0.3/num_chr;
            TickSize        = -0.005;  %negative for outside, percentage of longest chr figure.
            
            %define colors for colorBars plot
            colorNoData = [1.0    1.0    1.0   ]; %used when no data is available for the bin.
            colorInit   = [0.5    0.5    0.5   ]; %external; used in blending at ends of chr.
            %define colors for angleplot
            [colorA,colorB, colorAB, colorAAB,colorABB, colorAAAB,colorABBB, colorAAAAB,colorAAABB,...
                colorAABBB,colorABBBB, colorAAAAAB,colorABBBBB, colorPeak,colorCutoff ] = DefineColors();
            if (show_unnassigned == true) || (show_uncalibrated == true)
                colorUn_Hom  = [0.5   0.5   1.0  ]; % homozygous unassigned.
                colorUn_DHet = [0.667 0.667 0.667]; % disomy homozygous unassigned.
                colorUn_THet = [0.583 0.583 0.833]; % trisomy homozygous unassigned.
            else
                colorUn_Hom  = [1.0   1.0   1.0  ]; % homozygous unassigned.
            end;
            
            % main figure calculation and generation.
            fig = figure(dataset);
            set(gcf, 'Position', [0 70 1024 600]); %ddddd
            for chromosome = [8 1:7]
                usedPlot0   = chr_SNPdata{chromosome,1};  % 'ab'
                usedPlot1   = chr_SNPdata{chromosome,2};  % 'a'/'aa'/'aaa'/'aaaa'/'aaaaa'/'aaaaaa'
                usedPlot2   = chr_SNPdata{chromosome,3};  % 'b'/'bb'/'bbb'/'bbbb'/'bbbbb'/'bbbbbb'
                usedPlot3   = chr_SNPdata{chromosome,4};  % 'aab'
                usedPlot4   = chr_SNPdata{chromosome,5};  % 'abb'
                usedPlot5   = chr_SNPdata{chromosome,6};  % 'aaab'
                usedPlot6   = chr_SNPdata{chromosome,7};  % 'abbb'
                usedPlot7   = chr_SNPdata{chromosome,8};  % 'aaaab'
                usedPlot8   = chr_SNPdata{chromosome,9};  % 'aaabb'
                usedPlot9   = chr_SNPdata{chromosome,10}; % 'aabbb'
                usedPlot10  = chr_SNPdata{chromosome,11}; % 'abbbb'
                usedPlot11  = chr_SNPdata{chromosome,12}; % 'aaaaab'
                usedPlot12  = chr_SNPdata{chromosome,13}; % 'abbbbb'
                usedPlot13  = chr_SNPdata{chromosome,14}; % hom unassigned.
                usedPlotCGH = chr_CGHdata{chromosome,2}*maxY/2;
                % make CGH histograms to the right of the main chr cartoons.
                if (HistPlot == true)
                    width = 0.020;
                    if (chromosome == 8)
                        bottom = base + (chromosome)*(height+vertical_margin);
                    else
                        bottom = base + (8-chromosome)*(height+vertical_margin);
                    end;
                    for segment = 1:length(datasetDetails{dataset}.chrCopyNum{chromosome})
                        subplot('Position',[(0.15 + Chr_max_width*chr_size(chromosome)/max(chr_size) + 0.005)+width*(segment-1) bottom width height]);
                        histAll = [];
                        histAll2 = [];
                        smoothed = [];
                        smoothed2 = [];
                        for i = round(1+length(usedPlotCGH)*datasetDetails{dataset}.chr_breaks{chromosome}(segment)):round(length(usedPlotCGH)*datasetDetails{dataset}.chr_breaks{chromosome}(segment+1))
                            if (strcmp(scale_type,'Ratio') == 1)
                                if (flow_ploidy{dataset} == 0) % no value; no scaling.
                                    y_ = usedPlotCGH(i);
                                else
                                    if (Yscale_nearest_even_ploidy == true)
                                        if (abs(2-flow_ploidy{dataset}) < abs(4-flow_ploidy{dataset})) % nearer to diploid
                                            y_ = usedPlotCGH(i)*flow_ploidy{dataset}/2;
                                        else
                                            y_ = usedPlotCGH(i)*flow_ploidy{dataset}/4;
                                        end;
                                    else
                                        y_ = usedPlotCGH(i)*flow_ploidy{dataset}/2;
                                    end;
                                end;
                            else % scale_type = Log2Ratio.
                                if (flow_ploidy{dataset} == 0) % no value; no scaling.
                                    y_ = usedPlotCGH(i)+maxY/2;
                                else
                                    y_ = log2(pow2(usedPlotCGH(i))*flow_ploidy{dataset}/2)+maxY/2;
                                end;
                            end;
                            histAll{segment}(i) = y_;
                        end;
                        % make a histogram of CGH data, then smooth it for display.
                        histAll{segment}(histAll{segment}==0) = [];
                        histAll{segment}(length(histAll{segment})+1) = 0;   % endpoints added to ensure histogram bounds.
                        histAll{segment}(length(histAll{segment})+1) = 15;
                        histAll{segment}(histAll{segment}<0) = [];
                        histAll{segment}(histAll{segment}>15) = 15;
                        smoothed{segment} = smooth_gaussian(hist(histAll{segment},300),5,20);
                        % make a smoothed version of just the endpoints used to ensure histogram bounds.
                        histAll2{segment}(1) = 0;
                        histAll2{segment}(2) = 15;
                        smoothed2{segment} = smooth_gaussian(hist(histAll2{segment},300),5,20)*4;
                        % subtract the smoothed endpoints from the histogram to remove the influence of the added endpoints.
                        smoothed{segment} = (smoothed{segment}-smoothed2{segment});
                        smoothed{segment} = smoothed{segment}/max(smoothed{segment});
                        
                        plot([0; 1],[50; 50],'color',[0.75 0.75 0.75]);
                        hold on;
                        plot([0; 1],[100; 100],'color',[0.50 0.50 0.50]);
                        plot([0; 1],[150; 150],'color',[0.75 0.75 0.75]);
                        area(smoothed{segment},1:length(smoothed{segment}),'FaceColor',[0 0 0]);
                        hold off;
                        set(gca,'YTick',[]);    set(gca,'XTick',[]);
                        xlim([0,1]);            ylim([0,200]);
                    end;
                end;

                % places chr copy number to the right of the main chr cartoons.
                if (ChrNum == true)
                    % subplot to show chromosome copy number value.
                    width = 0.020;
                    if (chromosome == 8)
                        bottom = base + (chromosome)*(height+vertical_margin);
                    else
                        bottom = base + (8-chromosome)*(height+vertical_margin);
                    end;
                    
                    if (HistPlot == true)
                        subplot('Position',[(0.15 + Chr_max_width*chr_size(chromosome)/max(chr_size) + 0.005 + width*(length(datasetDetails{dataset}.chrCopyNum{chromosome})-1) + width+0.001) bottom width height]);
                    else
                        subplot('Position',[(0.15 + Chr_max_width*chr_size(chromosome)/max(chr_size) + 0.005) bottom width height]);
                    end;
                    axis off square;
                    set(gca,'YTick',[]);
                    set(gca,'XTick',[]);
                    if (length(datasetDetails{dataset}.chrCopyNum{chromosome}) == 1)
                        chr_string = num2str(datasetDetails{dataset}.chrCopyNum{chromosome}(1));
                    else
                        chr_string = num2str(datasetDetails{dataset}.chrCopyNum{chromosome}(1));
                        for i = 2:length(datasetDetails{dataset}.chrCopyNum{chromosome})
                            chr_string = [chr_string ',' num2str(datasetDetails{dataset}.chrCopyNum{chromosome}(i))];
                        end;
                    end;
                    text(0.1,0.5, chr_string,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',24 );
                end;

                % make allelic fraction histograms to the left of the main chr cartoons.
                if (AnglePlot == true)
                    width = 0.075;
                    if (chromosome == 8);  bottom = base + (chromosome)*(height+vertical_margin);
                    else                   bottom = base + (8-chromosome)*(height+vertical_margin);
                    end;
                    height = 0.5/num_chr;
                    
                    for segment = 1:length(datasetDetails{dataset}.chrCopyNum{chromosome})
                        if (segment == 1)
                            subplot('Position',[0.03 bottom width (height/length(datasetDetails{dataset}.chrCopyNum{chromosome}))]);
                        else
                            subplot('Position',[0.03 (bottom+height/length(datasetDetails{dataset}.chrCopyNum{chromosome})*(segment-1)) width (height/length(datasetDetails{dataset}.chrCopyNum{chromosome}))]);
                        end;
                        histAll = [];
                        histAll2 = [];
                        %for i = round(1+length(usedPlotCGH)*chr_breaks{chromosome}(segment))
                        %        round(  length(usedPlotCGH)*chr_breaks{chromosome}(segment+1))
                        for i = 1:2:SNP_probeset_length
                            if (probeset1(i).probe_chromosome == chromosome) && ...
                                    (probeset1(i).probe_location >= round(1+chr_size(chromosome)*datasetDetails{dataset}.chr_breaks{chromosome}(segment))) && ...
                                    (probeset1(i).probe_location <= round(chr_size(chromosome)*datasetDetails{dataset}.chr_breaks{chromosome}(segment+1)))
                                if (isnan(probeset1(i).probe_Ratio(dataset)) == 0) && (isnan(probeset1(i+1).probe_Ratio(dataset)) == 0)
                                    % Calculate value of SNP probe pair.
                                    [UsedData] = calculateValue(probeset1,i,DataTypeToUse,dataset);

                                    if (probeset1(i).probe_polarity == 0)
                                        if (show_unnassigned == true)
                                            histAll(i) = UsedData;
                                        else
                                            histAll(i) = 0;
                                        end;
                                    elseif (probeset1(i).probe_polarity == 4)
                                        % null action for when (probe_polarity == 4)
                                        % due to probe design error; probes are identical.
                                    else
                                        histAll(i) = UsedData;
                                    end;
                                end;
                            end;
                        end;
                        % make a histogram of SNP allelic fractions in segment, then smooth for display.
                        histAll(histAll==0) = [];
                        histAll(length(histAll)+1) = 0;
                        histAll(length(histAll)+1) = 1;
                        smoothed = smooth_gaussian(hist(histAll,200),3,20);
                        % make a smoothed version of just the endpoints used to ensure histogram bounds.
                        histAll2(1) = 0;
                        histAll2(2) = 1;
                        smoothed2 = smooth_gaussian(hist(histAll2,200),3,20);
                        % subtract the smoothed endpoints from the histogram to remove the influence of the added endpoints.
                        smoothed = (smoothed-smoothed2);
                        smoothed = smoothed/max(smoothed);
                        
                        hold on;
                        copynum    = datasetDetails{dataset}.chrCopyNum{chromosome}(segment);
                        region_    = 0;
                        for region = datasetDetails{dataset}.importantGaussians{chromosome,segment}
                            region_ = region_+1;
                            if (FillColors == true)
                                if (show_uncalibrated == true)
                                    color = colorUn_Hom;
                                else
                                    if (copynum <= 1) %monosomy
                                        if (region == 1); color = colorA;
                                        else              color = colorB;
                                        end;
                                        if (segment == 1)
                                            set(gca,'XTick',[0 200]);
                                            set(gca,'XTickLabel',{'a','b'});
                                        end;
                                    elseif (copynum <= 2) %disomy
                                        if (region == 1);     color = colorA;
                                        elseif (region == 2); color = colorAB;
                                        else                  color = colorB;
                                        end;
                                        if (segment == 1)
                                            set(gca,'XTick',0:100:200);
                                            set(gca,'XTickLabel',{'a','ab','b'});
                                        end;
                                    elseif (copynum <= 3) %trisomy
                                        if (region == 1);     color = colorA;
                                        elseif (region == 2); color = colorAAB;
                                        elseif (region == 3); color = colorABB;
                                        else                  color = colorB;
                                        end;
                                        if (segment == 1)
                                            set(gca,'XTick',[0 66.667 133.333 200]);
                                            set(gca,'XTickLabel',{'a','aab','abb','b'});
                                        end;
                                    elseif (copynum <= 4) %tetrasomy
                                        if (region == 1);     color = colorA;
                                        elseif (region == 2); color = colorAAAB;
                                        elseif (region == 3); color = colorAB;
                                        elseif (region == 4); color = colorABBB;
                                        else                  color = colorB;
                                        end;
                                        if (segment == 1)
                                            set(gca,'XTick',0:50:200);
                                            set(gca,'XTickLabel',{'a', 'aaab', 'ab', 'abbb' 'b'});
                                        end;
                                    elseif (copynum <= 5) %pentasomy
                                        if (region == 1);     color = colorA;
                                        elseif (region == 2); color = colorAAAAB;
                                        elseif (region == 3); color = colorAAABB;
                                        elseif (region == 4); color = colorAABBB;
                                        elseif (region == 5); color = colorABBBB;
                                        else                  color = colorB;
                                        end;
                                        if (segment == 1)
                                            set(gca,'XTick',0:40:200);
                                            set(gca,'XTickLabel',{'a', 'aaaab', 'aaabb', 'aabbb', 'abbbb' 'b'});
                                        end;
                                    else % if (copynum <= 6) %hexasomy
                                        if (region == 1);     color = colorA;
                                        elseif (region == 2); color = colorAAAAAB;
                                        elseif (region == 3); color = colorAAB;
                                        elseif (region == 4); color = colorAB;
                                        elseif (region == 5); color = colorABB;
                                        elseif (region == 6); color = colorABBBBB;
                                        else                  color = colorB;
                                        end;
                                        if (segment == 1)
                                            set(gca,'XTick',0:33.333:200);
                                            set(gca,'XTickLabel',{'a', 'aaaaab', 'aab', 'ab', 'abb', 'abbbbb' 'b'});
                                        end;
                                    end;
                                end;
                            else
                                color = colorAB;
                            end;
                            if (length(datasetDetails{dataset}.importantGaussians{chromosome,segment}) == 1)
                                area(1:200,smoothed(1:200),'FaceColor',color,'EdgeColor',color);
                            else
                                if (region_ == 1)
                                    area(        1:                                                                     ...
                                                 round(datasetDetails{dataset}.cutoffs{chromosome,segment}(region_  )*200), ...
                                        smoothed(1:                                                                     ...
                                                 round(datasetDetails{dataset}.cutoffs{chromosome,segment}(region_  )*200)),'FaceColor',color,'EdgeColor',color);
                                elseif (region_ == length(datasetDetails{dataset}.importantGaussians{chromosome,segment}))
                                    area(        round(datasetDetails{dataset}.cutoffs{chromosome,segment}(region_-1)*200): ...
                                                 200                                                                  , ...
                                        smoothed(round(datasetDetails{dataset}.cutoffs{chromosome,segment}(region_-1)*200): ...
                                                 200                                                                  ),'FaceColor',color,'EdgeColor',color);
                                else
                                    area(        round(datasetDetails{dataset}.cutoffs{chromosome,segment}(region_-1)*200): ...
                                                 round(datasetDetails{dataset}.cutoffs{chromosome,segment}(region_  )*200), ...
                                        smoothed(round(datasetDetails{dataset}.cutoffs{chromosome,segment}(region_-1)*200): ...
                                                 round(datasetDetails{dataset}.cutoffs{chromosome,segment}(region_  )*200)),'FaceColor',color,'EdgeColor',color);
                                end;
                            end;
                        end;
                        for peak = 1:length(datasetDetails{dataset}.peaks{chromosome,segment})
                            plot([datasetDetails{dataset}.peaks{chromosome,segment}(peak)*200; datasetDetails{dataset}.peaks{chromosome,segment}(peak)*200],[0; 1],'color',colorPeak);
                        end;
                        for cutoff = 1:length(datasetDetails{dataset}.cutoffs{chromosome,segment})
                            plot([datasetDetails{dataset}.cutoffs{chromosome,segment}(cutoff)*200; datasetDetails{dataset}.cutoffs{chromosome,segment}(cutoff)*200],[0; 1],'color',colorCutoff);
                        end;
                        set(gca,'FontSize',8);
                        hold off;
                        set(gca,'YTick',[]);
                        if (segment ~= 1)
                            set(gca,'XTick',[]);
                        end;
                        xlim([0,200]);
                        ylim([0,1]);
                    end;
                end;

                % make standard chromosome cartoons.
                width = Chr_max_width*chr_size(chromosome)/max(chr_size);
                if (chromosome == 8)
                    bottom = base + (chromosome)*(height+vertical_margin);
                else
                    bottom = base + (8-chromosome)*(height+vertical_margin);
                end;
                subplot('Position',[left bottom width height]);
                hold on;
                if (colorBars == true)
                    c_prev = colorInit;
                    c_post = colorInit;
                    c_     = c_prev;
                    infill = zeros(1,length(usedPlot0));
                    colors = [];
                    % determines the color of each bin.
                    for i = 1:length(usedPlot0)+1;
                        if (i-1 < length(usedPlot0))
                            c_tot_post = usedPlot0(i)+usedPlot1(i)+usedPlot2(i)+usedPlot3(i)+usedPlot4(i)+usedPlot5(i)+usedPlot6(i)+usedPlot7(i) ...
                                +usedPlot8(i)+usedPlot9(i)+usedPlot10(i)+usedPlot11(i)+usedPlot12(i)+usedPlot13(i);
                            if (c_tot_post == 0)
                                c_post = colorNoData;
                                infill(i) = 1;
                            else
                                if (show_uncalibrated == true)
                                    c_post = colorUn_DHet *usedPlot0(i)/c_tot_post + ...
                                             colorUn_HomAA*usedPlot1(i)/c_tot_post + ...
                                             colorUn_Hom  *usedPlot2(i)/c_tot_post + ...
                                             colorUn_THet *usedPlot3(i)/c_tot_post + ...
                                             colorUn_THet *usedPlot4(i)/c_tot_post + ...
                                             colorUn_Hom  *usedPlot5(i)/c_tot_post;
                                else
                                    c_post = colorAB      *usedPlot0(i)/c_tot_post + ...
                                             colorA       *usedPlot1(i)/c_tot_post + ...
                                             colorB       *usedPlot2(i)/c_tot_post + ...
                                             colorAAB     *usedPlot3(i)/c_tot_post + ...
                                             colorABB     *usedPlot4(i)/c_tot_post + ...
                                             colorAAAB    *usedPlot5(i)/c_tot_post + ...
                                             colorABBB    *usedPlot6(i)/c_tot_post + ...
                                             colorAAAAB   *usedPlot7(i)/c_tot_post + ...
                                             colorAAABB   *usedPlot8(i)/c_tot_post + ...
                                             colorAABBB   *usedPlot9(i)/c_tot_post + ...
                                             colorABBBB   *usedPlot10(i)/c_tot_post + ...
                                             colorAAAAAB  *usedPlot11(i)/c_tot_post + ...
                                             colorABBBBB  *usedPlot12(i)/c_tot_post + ...
                                             colorUn_Hom  *usedPlot13(i)/c_tot_post;
                                end;
                                infill(i) = 0;
                            end;
                        else
                            c_post = colorInit;
                            infill(i) = 1;
                        end;
                        colors(i,1) = c_post(1);
                        colors(i,2) = c_post(2);
                        colors(i,3) = c_post(3);
                    end;
                    % bleeds colors over adjacent white space.
                    if (infillColorBars == true)
                        for i = 1:length(usedPlot0)+1;
                            if (infill(i) == 1)
                                %look left.
                                foundLeft = false;
                                endLeft   = false;
                                deltaLeft = 1;
                                while (foundLeft == false)
                                    if (i-deltaLeft == 0)
                                        colorLeft(1) = colorNoData(1);
                                        colorLeft(2) = colorNoData(2);
                                        colorLeft(3) = colorNoData(3);
                                        foundLeft = true;
                                        endLeft   = true;
                                    elseif (infill(i-deltaLeft) == 0)
                                        colorLeft(1) = colors(i-deltaLeft,1);
                                        colorLeft(2) = colors(i-deltaLeft,2);
                                        colorLeft(3) = colors(i-deltaLeft,3);
                                        foundLeft = true;
                                    else
                                        deltaLeft = deltaLeft+1;
                                        %foundLeft = true;
                                    end;
                                end;
                                %look right.
                                foundRight = false;
                                endRight   = false;
                                deltaRight = 1;
                                while (foundRight == false)
                                    if (i+deltaRight == length(usedPlot0)+2)
                                        colorRight(1) = colorNoData(1);
                                        colorRight(2) = colorNoData(2);
                                        colorRight(3) = colorNoData(3);
                                        foundRight = true;
                                        endRight   = true;
                                    elseif (infill(i+deltaRight) == 0)
                                        colorRight(1) = colors(i+deltaRight,1);
                                        colorRight(2) = colors(i+deltaRight,2);
                                        colorRight(3) = colors(i+deltaRight,3);
                                        foundRight = true;
                                    else
                                        deltaRight = deltaRight+1;
                                        %foundRight = true;
                                    end;
                                end;
                                %make this the closer color.
                                if (endLeft == true)
                                    colors(i,1) = (colorRight(1)+3)/4;
                                    colors(i,2) = (colorRight(2)+3)/4;
                                    colors(i,3) = (colorRight(3)+3)/4;
                                elseif (endRight == true)
                                    colors(i,1) = (colorLeft(1)+3)/4;
                                    colors(i,2) = (colorLeft(2)+3)/4;
                                    colors(i,3) = (colorLeft(3)+3)/4;
                                elseif (deltaLeft < deltaRight)
                                    colors(i,1) = (colorLeft(1)+3)/4;
                                    colors(i,2) = (colorLeft(2)+3)/4;
                                    colors(i,3) = (colorLeft(3)+3)/4;
                                else
                                    colors(i,1) = (colorRight(1)+3)/4;
                                    colors(i,2) = (colorRight(2)+3)/4;
                                    colors(i,3) = (colorRight(3)+3)/4;
                                end;
                            end;
                        end;
                    end;
                    % draw colorbars.
                    for i = 1:length(usedPlot0)+1;
                        x_ = [i i i-1 i-1];
                        y_ = [0 maxY maxY 0];
                        c_post(1) = colors(i,1);
                        c_post(2) = colors(i,2);
                        c_post(3) = colors(i,3);
                        % makes a colorBar for each bin, using local smoothing
                        if (c_(1) > 1); c_(1) = 1; end;
                        if (c_(2) > 1); c_(2) = 1; end;
                        if (c_(3) > 1); c_(3) = 1; end;
                        if (blendColorBars == false)
                            f = fill(x_,y_,c_);
                        else
                            f = fill(x_,y_,c_/2+c_prev/4+c_post/4);
                        end;
                        c_prev = c_;
                        c_     = c_post;
                        set(f,'linestyle','none');
                    end;
                end;
                if (CGH_Genomic_display == true)
                    % cgh plot section.
                    c_ = [0 0 0];
                    for i = 1:length(usedPlotCGH);
                        x_ = [i i i-1 i-1];
                        if (strcmp(scale_type,'Ratio') == 1)
                            if (flow_ploidy{dataset} == 0) % no value; no scaling.
                                y_ = [maxY/2 usedPlotCGH(i) usedPlotCGH(i) maxY/2];
                            else
                                if (Yscale_nearest_even_ploidy == true)
                                    if (abs(2-flow_ploidy{dataset}) < abs(4-flow_ploidy{dataset})) % nearer to diploid
                                        y_ = [maxY/2 usedPlotCGH(i)*flow_ploidy{dataset}/2 usedPlotCGH(i)*flow_ploidy{dataset}/2 maxY/2];
                                    else
                                        y_ = [maxY/2 usedPlotCGH(i)*flow_ploidy{dataset}/4 usedPlotCGH(i)*flow_ploidy{dataset}/4 maxY/2];
                                    end;
                                else
                                    y_ = [maxY/2 usedPlotCGH(i)*flow_ploidy{dataset}/2 usedPlotCGH(i)*flow_ploidy{dataset}/2 maxY/2];
                                end;
                            end;
                        else % scale_type = Log2Ratio.
                            if (flow_ploidy{dataset} == 0) % no value; no scaling.
                                y_ = [maxY/2 usedPlotCGH(i)+maxY/2 usedPlotCGH(i)+maxY/2 maxY/2];
                            else
                                %y_ = [maxY/2 usedPlotCGH(i)+maxY/2 usedPlotCGH(i)+maxY/2 maxY/2];
                                y_ = [maxY/2 log2(pow2(usedPlotCGH(i))*flow_ploidy{dataset}/2)+maxY/2 ...
                                    log2(pow2(usedPlotCGH(i))*flow_ploidy{dataset}/2)+maxY/2 maxY/2];
                            end;
                        end;
                        % makes a blackbar for each bin.
                        f = fill(x_,y_,c_);
                        set(f,'linestyle','none');
                    end;
                    x2 = chr_size(chromosome)*chr_length_scale_multiplier;
                    if (strcmp(scale_type,'Ratio') == 1)
                        plot([0; x2], [maxY/2;maxY/2],'color',[0 0 0]);  % 2n line.
                        patch([0 x2], [maxY/4*3 maxY/4*3], 'b','Edgecolor',[0 0 0],'EdgeAlpha',0.15);  % top line.
                        patch([0 x2], [maxY/4 maxY/4], 'b','Edgecolor',[0 0 0],'EdgeAlpha',0.15);  % bottom line.
                    else
                        plot([0; x2], [maxY/2;maxY/2],'color',[0 0 0]);  % 2n line.
                        patch([0 x2], [maxY/2*log2(3) maxY/2*log2(3)], 'b','Edgecolor',[0 0 0],'EdgeAlpha',0.15);  % top cen.
                    end;
                    %end cgh plot section.
                end;
                %show centromere.
                x1 = cen_start(chromosome)*chr_length_scale_multiplier;
                x2 = cen_end(chromosome)*chr_length_scale_multiplier;
                leftEnd  = 0.5*(5000/bases_per_bin);
                rightEnd = chr_size(chromosome)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
                if (Centromere_format == 0)
                    h = fill([x1-maxY/2 x1 x2 x2+maxY/2], [0 maxY/4 maxY/4 0], [1 1 0]);
                    set(h,'FaceAlpha',0.5);
                    h = fill([x1-maxY/2 x1 x2 x2+maxY/2], [maxY maxY*3/4 maxY*3/4 maxY], [1 1 0]);
                    set(h,'FaceAlpha',0.5);
                elseif (Centromere_format == 1)
                    h = fill([x1-maxY/2 x1 x2 x2+maxY/2 x2 x1], [maxY/2 maxY*3/4 maxY*3/4 maxY/2 maxY/4 maxY/4], [1 1 0]);
                    set(h,'FaceAlpha',0.5);
                elseif (Centromere_format == 2)
                    dx = 5*(5000/bases_per_bin);
                    dy = 2;
                    patch([(x1-dx) (x1-dx) leftEnd      leftEnd      rightEnd      rightEnd;...
                        x1      x1      leftEnd      leftEnd      rightEnd      rightEnd;...
                        x2      x2      (leftEnd+dx) (leftEnd+dx) (rightEnd-dx) (rightEnd-dx);...
                        (x2+dx) (x2+dx) (leftEnd+dx) (leftEnd+dx) (rightEnd-dx) (rightEnd-dx)],...
                        [maxY      0  maxY      0  maxY      0;...
                        (maxY-dy) dy (maxY-dy) dy (maxY-dy) dy;...
                        (maxY-dy) dy maxY      0  maxY      0;...
                        maxY      0  maxY      0  maxY      0],...
                        'w','Edgecolor',[1 1 1]);  % top cen.
                    patch( ...
                        [leftEnd leftEnd (leftEnd+dx) x1-dx x1      x2      x2+dx rightEnd-dx rightEnd rightEnd rightEnd-dx x2+dx x2 x1 x1-dx leftEnd+dx], ...
                        [dy      maxY-dy maxY         maxY  maxY-dy maxY-dy maxY  maxY        maxY-dy  dy       0           0     dy dy 0     0         ], ...
                        [1 1 1],'FaceAlpha',0);
                end;
                %end show centromere.
                
                %show MRS locations
                if (show_MRS)
                    plot([leftEnd rightEnd], [-1.5 -1.5],'color',[0 0 0]);
                    hold on;
                    for i = 1:length(MRS_location)
                        if (MRS_chr(i) == chromosome)
                            MRSloc = MRS_location(i)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
                            plot(MRSloc,-1.5,'k:o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5)
                        end;
                    end;
                    hold off;
                end;
                %end show MRS locations.
                
                hold off;
                xlim([0,chr_size(chromosome)*chr_length_scale_multiplier]);
                if (show_MRS == true)
                    ylim([-1.5,maxY]);
                else
                    ylim([0,maxY]);
                end;
                set(gca,'YTick',[]);
                set(gca,'TickLength',[(TickSize*chr_size(1)/chr_size(chromosome)) 0]); %ensures same tick size on all subfigs.
                ylabel(chr_labels(chromosome), 'Rotation', 90, 'HorizontalAlign', 'center', 'VerticalAlign', 'bottom');
                set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
                if (show_Xlabels == true)
                    set(gca,'XTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8','2.0','2.2','2.4','2.6','2.8','3.0','3.2'});
                else
                    set(gca,'XTickLabel',[]);
                end;
                if (CGH_Genomic_display == true)
                    if (strcmp(scale_type,'Ratio') == 1)
                        set(gca,'YTick',[maxY/4 maxY/2 maxY/4*3 maxY]);
                        if (Yscale_nearest_even_ploidy == true)
                            if (abs(2-flow_ploidy{dataset}) < abs(4-flow_ploidy{dataset})) % nearer to diploid
                                set(gca,'YTickLabel',{'1','2','3','4'});
                            else % nearer to tetraploid
                                set(gca,'YTickLabel',{'2','4','6','8'});
                            end;
                        else
                            set(gca,'YTickLabel',{'1','2','3','4'});
                        end;
                    else
                        set(gca,'YTick',[0 (maxY/2) maxY/2*log2(3) maxY]);
                        set(gca,'YTickLabel',{'1','2','3','4'});
                    end;
                end;
                set(gca,'FontSize',8);
                if (chromosome == 8)
                    title(names{dataset},'Interpreter','none','FontSize',12);
                end;
            end;
            
            % Main figure %LOH and ploidy value display.
            subplot('Position',[0.65 0.0 0.1 0.1]);
            axis off square;
            xlim([0,1]);
            ylim([-0.1,1.2]);
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
            if (show_unnassigned == false)
                if (Show_Genomic_LOH_fraction == true)
                    text(0.35,2.5,['LOH_{snp/cgh} = ' num2str(100*SNPs_hom/SNPs_total) '%']);
                    if (flow_ploidy{dataset} ~= 0) % ploidy value given.
                        text(0.35,2.0,['Ploidy_{' flow_ploidy_type{dataset} '} = ' num2str(flow_ploidy{dataset}) 'n']);
                    end;
                else
                    if (flow_ploidy{dataset} ~= 0) % ploidy value given.
                        text(0.35,2.5,['Ploidy_{' flow_ploidy_type{dataset} '} = ' num2str(flow_ploidy{dataset}) 'n']);
                    end;
                end;
            else
                if (Show_Genomic_LOH_fraction == true)
                    text(0.35,2.0,['LOH_{snp/cgh} = ' num2str(100*SNPs_hom/SNPs_total) '%']);
                    if (flow_ploidy{dataset} ~= 0) % ploidy value given.
                        text(0.35,1.5,['Ploidy_{' flow_ploidy_type{dataset} '} = ' num2str(flow_ploidy{dataset}) 'n']);
                    end;
                else
                    if (flow_ploidy{dataset} ~= 0) % ploidy value given.
                        text(0.35,1.5,['Ploidy_{' flow_ploidy_type{dataset} '} = ' num2str(flow_ploidy{dataset}) 'n']);
                    end;
                end;
            end;
            
            % Main figure colors key.
            subplot('Position',[0.65 0.2 0.2 0.4]);
            axis off square;
            xlim([-0.1,1]);
            ylim([-0.1,1.6]);
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
            hold on;
            if (show_uncalibrated == true)
                patch([0 0.2 0.2 0], [1.4 1.4 1.5 1.5], [1 1 1]);      text(0.3,1.45,'No data');
                patch([0 0.2 0.2 0], [1.2 1.2 1.3 1.3], [0 0 0]);      text(0.3,1.25,'CGH ratios');
                patch([0 0.2 0.2 0], [1.0 1.0 1.1 1.1], colorUn_Hom);  text(0.3,1.05,'Homozygous unassigned');
                patch([0 0.2 0.2 0], [0.8 0.8 0.9 0.9], colorUn_DHet); text(0.3,0.85,'Heterozygous disomy');
                patch([0 0.2 0.2 0], [0.6 0.6 0.7 0.7], colorUn_THet); text(0.3,0.65,'Heterozygous trisomy unassigned');
                if (show_MRS == true)
                    plot([0 0.2],[0.45 0.45],'color',[0 0 0]);
                    plot(0.1, 0.45, 'k:o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5)
                    text(0.3,0.45,'MRS locations');
                end;
            else
                patch([0 0.2 0.2 0], [1.4 1.4 1.5 1.5], [1 1 1]);  text(0.3,1.45,'No data');
                patch([0 0.2 0.2 0], [1.2 1.2 1.3 1.3], [0 0 0]);  text(0.3,1.25,'CGH ratios');
                patch([0 0.2 0.2 0], [1.0 1.0 1.1 1.1], colorA);   text(0.3,1.05,'Homozygous homolog a');
                patch([0 0.2 0.2 0], [0.8 0.8 0.9 0.9], colorB);   text(0.3,0.85,'Homozygous homolog b');
                patch([0 0.2 0.2 0], [0.6 0.6 0.7 0.7], colorAB);  text(0.3,0.65,'Heterozygous disomy');
                patch([0 0.2 0.2 0], [0.4 0.4 0.5 0.5], colorAAB); text(0.3,0.45,'Heterozygous trisomy aab');
                patch([0 0.2 0.2 0], [0.2 0.2 0.3 0.3], colorABB); text(0.3,0.25,'Heterozygous trisomy abb');
                if (show_MRS == true)
                    plot([0 0.2],[0.05 0.05],'color',[0 0 0]);
                    plot(0.1, 0.05, 'k:o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5)
                    text(0.3,0.05,'MRS locations');
                end;
            end;
            hold off;
            
            % save then delete figures.
            if ispc  % Windows
                fig_dir = 'figures\snp-cgh\';
                if (isdir([file_dir 'figures\snp-cgh']) == 0)
                    mkdir([file_dir 'figures\snp-cgh']);
                end;
            else     % MacOS
                fig_dir = 'figures/snp-cgh/';
                if (isdir([file_dir 'figures/snp-cgh']) == 0)
                    mkdir([file_dir 'figures/snp-cgh']);
                end;
            end;
            %saveas(fig, [file_dir fig_dir strrep(names{dataset},' ','_') '.eps'], 'epsc');
            saveas(fig, [file_dir fig_dir strrep(names{dataset},' ','_') '.png'], 'png');
            delete(fig);
        end;
    end;
    
    %% ========================================================================
    % Plotting probes across genome by interpretation catagory after accounting
    % for polarity of SNP pairs.
    %    CGH info plotted as histogram.
    %    SNP info plotted as colorbar.
    % Only plot reference data against itself.
    %==========================================================================
    
    if (SC5314_display == 1)
        % define number of and labels for chromosomes.
        chr_labels = {'Chr1','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7','ChrR'};
        num_chr = length(chr_labels);
        
        bases_per_bin = SNP_bases_per_bin;
        chr_length_scale_multiplier = 1/bases_per_bin;
        
        [CGH_data_all SNP_data_all] = Load_reference_data_1();
        
        for dataset1 = 1:length(SNP_data_all(1).probe_data)
            for dataset2 = dataset1+1:length(SNP_data_all(1).probe_data)
                % Initializes vectors used to hold number of SNPs in each
                %    interpretation catagory.
                for i = 1:8   % eight chromosomes.
                    for j = 1:4   % four SNP interpretation catagories tracked.
                        chr_SNPdata{i,j} = zeros(1,ceil(chr_size(i)/bases_per_bin));
                    end;
                    for j = 1:2   % two CGH data catagories tracked.
                        chr_CGHdata{i,j} = zeros(1,ceil(chr_size(i)/bases_per_bin));
                    end;
                end;
                
                h = waitbar(0.0,'Please wait...','CreateCancelBtn','delete(h); clear h');
                set(h,'Name','Plotting SNP data.');
                counter = 0;
                SNPs_hom = 0;
                SNPs_total = 0;
                
                for i = 1:SNP_probeset_length
                    % probeset1 consists of sequential pairs, representing each allele for each SNP.
                    %    The odd numbered probes are the first of each pair.
                    if (mod(i,2) == 1)
                        counter = counter+1;
                        if (counter == 500)
                            waitbar(i/SNP_probeset_length,h,[num2str(i) '/' num2str(SNP_probeset_length)]);
                            counter = 0;
                        end;
                        
                        % determines if probe pair is useful; both probes have data.
                        if (isnan(SNP_data_all(i).probe_data(1)) == 0) && (isnan(SNP_data_all(i+1).probe_data(2)) == 0)
                            % assigns interpretations to SNP pairs.
                            % 0 : heterozygous.
                            % 1 : homozygous, "b:cyan".
                            % 2 : homozygous, "a:magenta".
                            % 3 : homozygous, unassigned.
                            if (SNP_data_all(i).probe_data(dataset1)/SNP_data_all(i).probe_data(dataset2) < SNP_data_all(i+1).probe_data(dataset1)/SNP_data_all(i+1).probe_data(dataset2)*yy/4)
                                SNP_data_all(i).probe_assignment   = 1;
                                SNP_data_all(i+1).probe_assignment = 1;
                                if (Show_Genomic_LOH_fraction == true)
                                    SNPs_hom = SNPs_hom+1;
                                    SNPs_total = SNPs_total+1;
                                end;
                            elseif (SNP_data_all(i+1).probe_data(dataset1)/SNP_data_all(i+1).probe_data(dataset2) < SNP_data_all(i).probe_data(dataset1)/SNP_data_all(i).probe_data(dataset2)*yy/4)
                                SNP_data_all(i).probe_assignment   = 2;
                                SNP_data_all(i+1).probe_assignment = 2;
                                if (Show_Genomic_LOH_fraction == true)
                                    SNPs_hom = SNPs_hom+1;
                                    SNPs_total = SNPs_total+1;
                                end;
                            else  % heterozygous regions.
                                SNP_data_all(i).probe_assignment   = 0;
                                SNP_data_all(i+1).probe_assignment = 0;
                                if (Show_Genomic_LOH_fraction == true)
                                    SNPs_total = SNPs_total+1;
                                end;
                            end;
                            
                            % val is the genomic location of the probePair.
                            val = ceil(SNP_data_all(i).probe_location/bases_per_bin);
                            
                            % Determines distribution of SNP pairs in four catagories.
                            % 0 : heterozygous.
                            % 1 : homozygous, "a:cyan".
                            % 2 : homozygous, "b:magenta".
                            % 3 : homozygous, unassigned.
                            % increment appropriate category.
                            chr_SNPdata{SNP_data_all(i).probe_chromosome,SNP_data_all(i).probe_assignment+1}(val) = chr_SNPdata{SNP_data_all(i).probe_chromosome,SNP_data_all(i).probe_assignment+1}(val)+1;
                        end;
                    end;
                end;
                if exist('h','var'); delete(h); clear h; end;
                h = waitbar(0.0,'Please wait...','CreateCancelBtn','delete(h); clear h');
                set(h,'Name','Plotting CGH data.');
                counter = 0;
                for i = 1:CGH_probeset_length
                    counter = counter+1;
                    if (counter == 500)
                        waitbar(i/CGH_probeset_length,h,[num2str(i) '/' num2str(CGH_probeset_length)]);
                        counter = 0;
                    end;
                    
                    % val is the genomic location of the probePair.
                    % val2 is the ratio of reference to experiment for the CGH probe.
                    val  = ceil((CGH_data_all(i).probe_location)/bases_per_bin);
                    if (strcmp(scale_type,'Ratio') == 1)
                        val2 = CGH_data_all(i).probe_data(dataset1)/CGH_data_all(i).probe_data(dataset2);
                    else
                        val2 = 0;
                    end;
                    
                    % Determines distribution of CGH data.
                    if (isnan(val2) ~= 1)
                        % count of data points at CGH locus.
                        chr_CGHdata{CGH_data_all(i).probe_chromosome,1}(val) = chr_CGHdata{CGH_data_all(i).probe_chromosome,1}(val)+1;
                        % total data at locus.
                        chr_CGHdata{CGH_data_all(i).probe_chromosome,2}(val) = chr_CGHdata{CGH_data_all(i).probe_chromosome,2}(val)+val2;
                    end;
                end;
                if exist('h','var'); delete(h); clear h; end;
                % divide total ave values by number of points per bin.
                for i = 1:8
                    for j = 1:length(chr_CGHdata{i,2})
                        if (chr_CGHdata{i,1}(j) == 0)
                            chr_CGHdata{i,2}(j) = 1;
                        else
                            %chr_CGHdata{i,2}(j) = chr_CGHdata{i,2}(j)/chr_CGHdata{i,1}(j);
                            chr_CGHdata{i,2}(j) = 1;
                        end;
                    end;
                end;
                
                % basic plot parameters.
                left            = 0.15;
                height          = 0.5/num_chr;
                base            = 0.1;
                vertical_margin = 0.3/num_chr;
                TickSize        = -0.005;  %negative for outside, percentage of longest chr figure.
                
                % define Y-max, for scaling plots.
                maxY = 10;
                
                %define colors for colorBars plot
                colorNoData = [1.0   1.0   1.0  ]; %used when no data is available for the bin.
                colorInit   = [0.5   0.5   0.5  ]; %external; used in blending at ends of chr.
                %define colors for angleplot
                [colorA,colorB, colorAB, colorAAB,colorABB, colorAAAB,colorABBB, colorAAAAB,colorAAABB,...
                    colorAABBB,colorABBBB, colorAAAAAB,colorABBBBB, colorPeak,colorCutoff ] = DefineColors();
                if (show_unnassigned == true) || (show_uncalibrated == true)
                    colorUn_Hom  = [0.5   0.5   1.0  ]; % homozygous unassigned.
                else
                    colorUn_Hom  = [1.0   1.0   1.0  ]; % homozygous unassigned.
                end;
                
                fig = figure;
                set(gcf, 'Position', [0 70 1024 600]);
                for chromosome = [8 1:7]
                    usedPlot0   = chr_SNPdata{chromosome,1}; % het.
                    usedPlot1   = chr_SNPdata{chromosome,2}; % hom a:magenta.
                    usedPlot2   = chr_SNPdata{chromosome,3}; % hom b:cyan.
                    usedPlot3   = chr_SNPdata{chromosome,4}; % hom unassigned.
                    usedPlotCGH = chr_CGHdata{chromosome,2}*maxY/2;
                    
                    if (length(find(chromosomes_to_analyze==chromosome)) == 1)
                        if (AnglePlot == true)
                            width = 0.075;
                            if (chromosome == 8)
                                bottom = base + (chromosome)*(height+vertical_margin);
                            else
                                bottom = base + (8-chromosome)*(height+vertical_margin);
                            end;
                            subplot('Position',[0.03 bottom width height]);
                            histAll = [];
                            for i = SNP_probeset_length:-1:1
                                if (mod(i,2) == 1)
                                    if (probeset1(i).probe_chromosome == chromosome)
                                        if (isnan(SNP_data_all(i).probe_data(dataset1)/SNP_data_all(i).probe_data(dataset2)) == 0) && (isnan(SNP_data_all(i+1).probe_data(dataset1)/SNP_data_all(i+1).probe_data(dataset2)) == 0)
                                            if (probeset1(i).probe_polarity == 1)
                                                histAll(i) = atan2(  SNP_data_all(i+1).probe_data(dataset1)/SNP_data_all(i+1).probe_data(dataset2) , SNP_data_all(i).probe_data(dataset1)/SNP_data_all(i).probe_data(dataset2)  );
                                            elseif (probeset1(i).probe_polarity == 2)
                                                histAll(i) = atan2(  SNP_data_all(i).probe_data(dataset1)/SNP_data_all(i).probe_data(dataset2) , SNP_data_all(i+1).probe_data(dataset1)/SNP_data_all(i+1).probe_data(dataset2)  );
                                            elseif (probeset1(i).probe_polarity == 0)
                                                if (show_unnassigned == true)
                                                    histAll(i) = atan2(  SNP_data_all(i+1).probe_data(dataset1)/SNP_data_all(i+1).probe_data(dataset2) , SNP_data_all(i).probe_data(dataset1)/SNP_data_all(i).probe_data(dataset2)  );
                                                else
                                                    histAll(i) = 0;
                                                end;
                                            else
                                                % null action for when (probe_polarity == 4)
                                                % due to probe design error; probes are identical.
                                            end;
                                        end;
                                    end;
                                end;
                            end;
                            
                            % make a histogram of SNP allele ratio angles, then
                            % smooth it for display.
                            histAll(histAll==0) = [];
                            histAll(length(histAll)+1) = 0;
                            histAll(length(histAll)+1) = pi/2;
                            n = 200;
                            h2 = hist(histAll,n);
                            f = fft(h2);
                            m = 160;
                            f(n/2+1-m/2:n/2+m/2) = zeros(m,1);
                            smoothed = real(ifft(f));
                            smoothed = smoothed/max(smoothed);
                            plot([100; 100],[0; 1],'color',[0.75 0.75 0.75]);
                            hold on;
                            if (datasetDetails{dataset}.chrCopyNum{chromosome}(1) == 3)
                                plot([realHomozygousAngle/22.5*50;     realHomozygousAngle/22.5*50    ],[0; 1],'color',[0.75 0.75 0.75]);
                                plot([200-realHomozygousAngle/22.5*50; 200-realHomozygousAngle/22.5*50],[0; 1],'color',[0.75 0.75 0.75]);
                                TrisomyAngle = (realHomozygousAngle+90)/3;
                                plot([TrisomyAngle/22.5*50;          TrisomyAngle/22.5*50         ],[0; 1],'color',[0.75 0.75 0.75]);
                                plot([200-TrisomyAngle/22.5*50;      200-TrisomyAngle/22.5*50     ],[0; 1],'color',[0.75 0.75 0.75]);
                                SNPcutoffAngle2 = ((realHomozygousAngle+90)/3+realHomozygousAngle)/2;
                                plot([SNPcutoffAngle2/22.5*50;          SNPcutoffAngle2/22.5*50         ],[0; 1],'color',[1.00 0.00 0.00]);
                                plot([200-SNPcutoffAngle2/22.5*50;      200-SNPcutoffAngle2/22.5*50     ],[0; 1],'color',[1.00 0.00 0.00]);
                                
                                %plot(smoothed,'color',[0 0 0]);
                                area(1:round(SNPcutoffAngle2/22.5*50),smoothed(1:round(SNPcutoffAngle2/22.5*50))            ,'FaceColor',colorA,'EdgeColor',colorA);
                                area(round(SNPcutoffAngle2/22.5*50):100,smoothed(round(SNPcutoffAngle2/22.5*50):100)        ,'FaceColor',colorAAB,'EdgeColor',colorAAB);
                                area(101:200-round(SNPcutoffAngle2/22.5*50),smoothed(101:200-round(SNPcutoffAngle2/22.5*50)),'FaceColor',colorABB,'EdgeColor',colorABB);
                                area(200-round(SNPcutoffAngle2/22.5*50):200,smoothed(200-round(SNPcutoffAngle2/22.5*50):200),'FaceColor',colorB,'EdgeColor',colorB);
                                set(gca,'XTick',[0 round(TrisomyAngle/22.5*50) 200-round(TrisomyAngle/22.5*50) 200]);
                                set(gca,'XTickLabel',{'a','aab','abb','b'});
                            else % datasetDetails{dataset}.chrCopyNum{chromosome}(1) == [1,2,4,etc]
                                plot([realHomozygousAngle/22.5*50;     realHomozygousAngle/22.5*50    ],[0; 1],'color',[0.75 0.75 0.75]);
                                plot([200-realHomozygousAngle/22.5*50; 200-realHomozygousAngle/22.5*50],[0; 1],'color',[0.75 0.75 0.75]);
                                SNPcutoffAngle = (realHomozygousAngle+45)/2;
                                plot([SNPcutoffAngle/22.5*50;          SNPcutoffAngle/22.5*50         ],[0; 1],'color',[0.75 0.75 0.75]);
                                plot([200-SNPcutoffAngle/22.5*50;      200-SNPcutoffAngle/22.5*50     ],[0; 1],'color',[0.75 0.75 0.75]);
                                
                                %plot(smoothed,'color',[0 0 0]);
                                area(1:round(SNPcutoffAngle/22.5*50),smoothed(1:round(SNPcutoffAngle/22.5*50)),                                                        'FaceColor',colorA,'EdgeColor',colorA);
                                area(round(SNPcutoffAngle/22.5*50):round(200-SNPcutoffAngle/22.5*50),smoothed(round(SNPcutoffAngle/22.5*50):round(200-SNPcutoffAngle/22.5*50)),'FaceColor',colorAB, 'EdgeColor',colorAB);
                                area(round(200-SNPcutoffAngle/22.5*50):200,smoothed(round(200-SNPcutoffAngle/22.5*50):200),                                            'FaceColor',colorB,'EdgeColor',colorB);
                                set(gca,'XTick',0:100:200);
                                set(gca,'XTickLabel',{'a','ab','b'});
                            end;
                            hold off;
                            set(gca,'YTick',[]);
                            set(gca,'FontSize',8);
                            xlim([0,n]);
                            ylim([0,1]);
                        end;
                        width = Chr_max_width*chr_size(chromosome)/max(chr_size);
                        if (chromosome == 8)
                            bottom = base + (chromosome)*(height+vertical_margin);
                        else
                            bottom = base + (8-chromosome)*(height+vertical_margin);
                        end;
                        subplot('Position',[left bottom width height]);
                        hold on;
                        if (colorBars == true)
                            c_prev = colorInit;
                            c_post = colorInit;
                            c_     = c_prev;
                            infill = zeros(1,length(usedPlot0));
                            colors = [];
                            for i = 1:length(usedPlot0)+1;
                                % determines the color of the next bin.
                                if (i-1 < length(usedPlot0))
                                    c_tot_post = usedPlot0(i)+usedPlot1(i)+usedPlot2(i)+usedPlot3(i);
                                    if (c_tot_post == 0)
                                        c_post = colorNoData;
                                        infill(i) = 1;
                                    else
                                        c_post = colorAB*usedPlot0(i)/c_tot_post + ...
                                            colorA*usedPlot1(i)/c_tot_post + ...
                                            colorB*usedPlot2(i)/c_tot_post + ...
                                            colorUn_Hom*usedPlot3(i)/c_tot_post;
                                        infill(i) = 0;
                                    end;
                                else
                                    c_post = colorInit;
                                    infill(i) = 1;
                                end;
                                colors(i,1) = c_post(1);
                                colors(i,2) = c_post(2);
                                colors(i,3) = c_post(3);
                            end;
                            if (infillColorBars == true)
                                for i = 1:length(usedPlot0)+1;
                                    if (infill(i) == 1)
                                        %look left.
                                        foundLeft = false;
                                        endLeft   = false;
                                        deltaLeft = 1;
                                        while (foundLeft == false)
                                            if (i-deltaLeft == 0)
                                                colorLeft(1) = colorNoData(1);
                                                colorLeft(2) = colorNoData(2);
                                                colorLeft(3) = colorNoData(3);
                                                foundLeft = true;
                                                endLeft   = true;
                                            elseif (infill(i-deltaLeft) == 0)
                                                colorLeft(1) = colors(i-deltaLeft,1);
                                                colorLeft(2) = colors(i-deltaLeft,2);
                                                colorLeft(3) = colors(i-deltaLeft,3);
                                                foundLeft = true;
                                            else
                                                deltaLeft = deltaLeft+1;
                                                %foundLeft = true;
                                            end;
                                        end;
                                        %look right.
                                        foundRight = false;
                                        endRight   = false;
                                        deltaRight = 1;
                                        while (foundRight == false)
                                            if (i+deltaRight == length(usedPlot0)+2)
                                                colorRight(1) = colorNoData(1);
                                                colorRight(2) = colorNoData(2);
                                                colorRight(3) = colorNoData(3);
                                                foundRight = true;
                                                endRight   = true;
                                            elseif (infill(i+deltaRight) == 0)
                                                colorRight(1) = colors(i+deltaRight,1);
                                                colorRight(2) = colors(i+deltaRight,2);
                                                colorRight(3) = colors(i+deltaRight,3);
                                                foundRight = true;
                                            else
                                                deltaRight = deltaRight+1;
                                                %foundRight = true;
                                            end;
                                        end;
                                        %make this the closer color.
                                        if (endLeft == true)
                                            colors(i,1) = (colorRight(1)+3)/4;
                                            colors(i,2) = (colorRight(2)+3)/4;
                                            colors(i,3) = (colorRight(3)+3)/4;
                                        elseif (endRight == true)
                                            colors(i,1) = (colorLeft(1)+3)/4;
                                            colors(i,2) = (colorLeft(2)+3)/4;
                                            colors(i,3) = (colorLeft(3)+3)/4;
                                        elseif (deltaLeft < deltaRight)
                                            colors(i,1) = (colorLeft(1)+3)/4;
                                            colors(i,2) = (colorLeft(2)+3)/4;
                                            colors(i,3) = (colorLeft(3)+3)/4;
                                        else
                                            colors(i,1) = (colorRight(1)+3)/4;
                                            colors(i,2) = (colorRight(2)+3)/4;
                                            colors(i,3) = (colorRight(3)+3)/4;
                                        end;
                                    end;
                                end;
                            end;
                            for i = 1:length(usedPlot0)+1;
                                x_ = [i i i-1 i-1];
                                y_ = [0 maxY maxY 0];
                                c_post(1) = colors(i,1);
                                c_post(2) = colors(i,2);
                                c_post(3) = colors(i,3);
                                % makes a colorBar for each bin, using local smoothing
                                if (blendColorBars == false)
                                    f = fill(x_,y_,c_);
                                else
                                    f = fill(x_,y_,c_/2+c_prev/4+c_post/4);
                                end;
                                c_prev = c_;
                                c_     = c_post;
                                set(f,'linestyle','none');
                            end;
                        end;
                        if (CGH_Genomic_display == true)
                            % cgh plot section.
                            c_ = [0 0 0];
                            for i = 1:length(usedPlotCGH);
                                x_ = [i i i-1 i-1];
                                if (strcmp(scale_type,'Ratio') == 1)
                                    y_ = [maxY/2 usedPlotCGH(i) usedPlotCGH(i) maxY/2];
                                else % scale_type = Log2Ratio.
                                    %y_ = [maxY/2 usedPlotCGH(i)+maxY/2 usedPlotCGH(i)+maxY/2 maxY/2];
                                    y_ = [maxY/2 log2(pow2(usedPlotCGH(i)))+maxY/2 ...
                                        log2(pow2(usedPlotCGH(i)))+maxY/2 maxY/2];
                                end;
                                % makes a blackbar for each bin.
                                f = fill(x_,y_,c_);
                                set(f,'linestyle','none');
                            end;
                            x2 = chr_size(chromosome)*chr_length_scale_multiplier;
                            if (strcmp(scale_type,'Ratio') == 1)
                                plot([0; x2],[maxY/2;maxY/2],'color',[0 0 0]);  % 2n line.
                                patch([0 x2], [maxY/4*3 maxY/4*3], 'b','Edgecolor',[0 0 0],'EdgeAlpha',0.15);  % top line.
                                patch([0 x2], [maxY/4 maxY/4], 'b','Edgecolor',[0 0 0],'EdgeAlpha',0.15);  % bottom line.
                            else
                                plot([0; x2],[maxY/2;maxY/2],'color',[0 0 0]);  % 2n line.
                                patch([0 x2], [maxY/2*log2(3) maxY/2*log2(3)], 'b','Edgecolor',[0 0 0],'EdgeAlpha',0.15);  % top cen.
                            end;
                            %end cgh plot section.
                        end;
                        %show centromere.
                        x1 = cen_start(chromosome)*chr_length_scale_multiplier;
                        x2 = cen_end(chromosome)*chr_length_scale_multiplier;
                        leftEnd  = 0.5*(5000/bases_per_bin);
                        rightEnd = chr_size(chromosome)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
                        if (Centromere_format == 0)
                            h = fill([x1-maxY/2 x1 x2 x2+maxY/2], [0 maxY/4 maxY/4 0], [1 1 0]);
                            set(h,'FaceAlpha',0.5);
                            h = fill([x1-maxY/2 x1 x2 x2+maxY/2], [maxY maxY*3/4 maxY*3/4 maxY], [1 1 0]);
                            set(h,'FaceAlpha',0.5);
                        elseif (Centromere_format == 1)
                            h = fill([x1-maxY/2 x1 x2 x2+maxY/2 x2 x1], [maxY/2 maxY*3/4 maxY*3/4 maxY/2 maxY/4 maxY/4], [1 1 0]);
                            set(h,'FaceAlpha',0.5);
                        elseif (Centromere_format == 2)
                            dx = 5*(5000/bases_per_bin);
                            dy = 2;
                            patch([(x1-dx) (x1-dx) leftEnd      leftEnd      rightEnd      rightEnd;...
                                x1      x1      leftEnd      leftEnd      rightEnd      rightEnd;...
                                x2      x2      (leftEnd+dx) (leftEnd+dx) (rightEnd-dx) (rightEnd-dx);...
                                (x2+dx) (x2+dx) (leftEnd+dx) (leftEnd+dx) (rightEnd-dx) (rightEnd-dx)],...
                                [maxY      0  maxY      0  maxY      0;...
                                (maxY-dy) dy (maxY-dy) dy (maxY-dy) dy;...
                                (maxY-dy) dy maxY      0  maxY      0;...
                                maxY      0  maxY      0  maxY      0],...
                                'w','Edgecolor',[1 1 1]);  % top cen.
                            patch( ...
                                [leftEnd leftEnd (leftEnd+dx) x1-dx x1      x2      x2+dx rightEnd-dx rightEnd rightEnd rightEnd-dx x2+dx x2 x1 x1-dx leftEnd+dx], ...
                                [dy      maxY-dy maxY         maxY  maxY-dy maxY-dy maxY  maxY        maxY-dy  dy       0           0     dy dy 0     0         ], ...
                                [1 1 1],'FaceAlpha',0);
                        end;
                        %end show centromere.
                        hold off;
                        xlim([0,chr_size(chromosome)*chr_length_scale_multiplier]);
                        ylim([0,maxY]);
                        set(gca,'YTick',[]);
                        set(gca,'TickLength',[(TickSize*chr_size(1)/chr_size(chromosome)) 0]); %ensures same tick size on all subfigs.
                        ylabel(chr_labels(chromosome), 'Rotation', 90, 'HorizontalAlign', 'center', 'VerticalAlign', 'bottom');
                        set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
                        set(gca,'XTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8','2.0','2.2','2.4','2.6','2.8','3.0','3.2'});
                        if (CGH_Genomic_display == true)
                            if (strcmp(scale_type,'Ratio') == 1)
                                set(gca,'YTick',[maxY/4 maxY/2 maxY/4*3 maxY]);
                                set(gca,'YTickLabel',{'1','2','3','4'});
                            else
                                set(gca,'YTick',[0 (maxY/2) maxY/2*log2(3) maxY]);
                                set(gca,'YTickLabel',{'1','2','3','4'});
                            end;
                        end;
                        set(gca,'FontSize',8);
                        if (chromosome == 8)
                            title(['SC5314-' num2str(dataset1) ' vs. SC5314-' num2str(dataset2) ],'Interpreter','none','FontSize',12);
                        end;
                    end;
                end;
                
                % Main figure %LOH and ploidy value display.
                subplot('Position',[0.7 0.1 0.1 0.1]);
                axis off square;
                xlim([0,1]);
                ylim([-0.1,1.2]);
                set(gca,'XTick',[]);
                set(gca,'YTick',[]);
                if (show_unnassigned == false)
                    if (Show_Genomic_LOH_fraction == true)
                        text(0.35,2.5,['LOH_s_n_p_/_c_g_h = ' num2str(100*SNPs_hom/SNPs_total) '%']);
                        text(0.35,2.0,'Ploidy_{reference} = 2n');
                    else
                        text(0.35,2.5,'Ploidy_{reference} = 2n');
                    end;
                else
                    if (Show_Genomic_LOH_fraction == true)
                        text(0.35,2.0,['LOH_s_n_p_/_c_g_h = ' num2str(100*SNPs_hom/SNPs_total) '%']);
                        text(0.35,1.5,'Ploidy_{reference} = 2n');
                    else
                        text(0.35,1.5,'Ploidy_{reference} = 2n');
                    end;
                end;
                
                % Main figure colors key.
                subplot('Position',[0.7 0.3 0.2 0.2]);
                axis off square;
                xlim([0,1]);
                ylim([-0.1,1.2]);
                set(gca,'XTick',[]);
                set(gca,'YTick',[]);
                if (show_uncalibrated == true)
                    patch([0 0.2 0.2 0], [1.0 1.0 1.1 1.1], [1   0   1  ]);   text(0.3,1.05,'Homozygous homolog a');
                    patch([0 0.2 0.2 0], [0.8 0.8 0.9 0.9], [0.5 0.5 1  ]);   text(0.3,0.85,'Homozygous unassigned');
                    patch([0 0.2 0.2 0], [0.6 0.6 0.7 0.7], [0   1   1  ]);   text(0.3,0.65,'Homozygous homolog b');
                    patch([0 0.2 0.2 0], [0.4 0.4 0.5 0.5], [2/3 2/3 2/3]);   text(0.3,0.45,'Heterozygous');
                    patch([0 0.2 0.2 0], [0.2 0.2 0.3 0.3], [1   1   1  ]);   text(0.3,0.25,'No data');
                    patch([0 0.2 0.2 0], [0.0 0.0 0.1 0.1], [0   0   0  ]);   text(0.3,0.05,'CGH ratios');
                else
                    patch([0 0.2 0.2 0], [1.0 1.0 1.1 1.1], [1   0   1  ]);   text(0.3,1.05,'Homozygous homolog a');
                    patch([0 0.2 0.2 0], [0.8 0.8 0.9 0.9], [0   1   1  ]);   text(0.3,0.85,'Homozygous homolog b');
                    patch([0 0.2 0.2 0], [0.6 0.6 0.7 0.7], [2/3 2/3 2/3]);   text(0.3,0.65,'Heterozygous');
                    patch([0 0.2 0.2 0], [0.4 0.4 0.5 0.5], [1   1   1  ]);   text(0.3,0.45,'No data');
                    patch([0 0.2 0.2 0], [0.2 0.2 0.3 0.3], [0   0   0  ]);   text(0.3,0.25,'CGH ratios');
                end;
                
                % save then delete figures.
                if ispc  % Windows
                    fig_dir = [pwd '\figures\'];
                    if (isdir([pwd '\figures\']) == 0)
                        mkdir([pwd '\figures\']);
                    end;
                else     % MacOS
                    fig_dir = [pwd '/figures/'];
                    if (isdir([pwd '/figures/']) == 0)
                        mkdir([pwd '/figures/']);
                    end;
                end;
                %saveas(fig, [fig_dir 'SC5314_' num2str(dataset1) '-' num2str(dataset2) '.eps'], 'epsc');
                saveas(fig, [fig_dir 'SC5314_' num2str(dataset1) '-' num2str(dataset2) '.png'], 'png');
                delete(fig);
            end;
        end;
        
        
        % clean up variables not needed further.
        clear plot* bases_per_bin num_chr base bottom height i ...
            vertical_margin width left maxY val counter x1 x2 ...
            data1* data2* data3* c_*;
    end;
    
    %% ========================================================================
    % Plotting probes across genome by interpretation catagory after accounting
    % for polarity of SNP pairs.
    %    CGH info plotted either as histogram.
    %    SNP info plotted as colorbar.
    % Zoomed in close-up view to examine actual data points.
    %==========================================================================
    if (Zoom_view == 1)
        for Zoom_chr = 1:8
            bases_per_bin = SNP_bases_per_bin;
            chr_length_scale_multiplier = 1/bases_per_bin;
            
            fig = figure;
            
            set(gcf, 'Position', [0 70 1024 600]);
            for dataset = 1:length(names);
                for i = 1:SNP_probeset_length
                    % probeset1 consists of sequential pairs, representing each allele for each SNP.
                    %    The odd numbered probes are the first of each pair.
                    if (mod(i,2) == 1)
                        % determines if probe pair is useful; both probes have data.
                        if (isnan(probeset1(i).probe_Ratio(dataset)) == 0) && (isnan(probeset1(i+1).probe_Ratio(dataset)) == 0)
                            % assigns interpretations to SNP pairs.
                            % 0 : heterozygous.
                            % 1 : homozygous, "b:cyan".
                            % 2 : homozygous, "a:magenta".
                            % 3 : homozygous, unassigned.
                            % 4 : unassigned interpretation.
                            
                            if (probeset1(i).probe_Ratio(dataset) < probeset1(i+1).probe_Ratio(dataset)*yy/4)
                                if (probeset1(i).probe_polarity == 1)
                                    probeset1(i).probe_assignment   = 2;
                                    probeset1(i+1).probe_assignment = 2;
                                elseif (probeset1(i).probe_polarity == 2)
                                    probeset1(i).probe_assignment   = 1;
                                    probeset1(i+1).probe_assignment = 1;
                                elseif (probeset1(i).probe_polarity == 0)
                                    probeset1(i).probe_assignment   = 3;
                                    probeset1(i+1).probe_assignment = 3;
                                else
                                    % null action for when (probe_polarity == 4)
                                    % due to probe design error; probes are identical.
                                end;
                            elseif (probeset1(i+1).probe_Ratio(dataset) < probeset1(i).probe_Ratio(dataset)*yy/4)
                                if (probeset1(i).probe_polarity == 1)
                                    probeset1(i).probe_assignment   = 1;
                                    probeset1(i+1).probe_assignment = 1;
                                elseif (probeset1(i).probe_polarity == 2)
                                    probeset1(i).probe_assignment   = 2;
                                    probeset1(i+1).probe_assignment = 2;
                                elseif (probeset1(i).probe_polarity == 0)
                                    probeset1(i).probe_assignment   = 3;
                                    probeset1(i+1).probe_assignment = 3;
                                else
                                    % null action for when (probe_polarity == 4)
                                    % due to probe design error; probes are identical.
                                end;
                            else  % heterozygous regions.
                                if (probeset1(i).probe_polarity == 0) || (probeset1(i).probe_polarity == 1) || (probeset1(i).probe_polarity == 2)
                                    probeset1(i).probe_assignment   = 0;
                                    probeset1(i+1).probe_assignment = 0;
                                else
                                    % null action for when (probe_polarity == 4)
                                    % due to probe design error; probes are identical.
                                end;
                            end;
                            
                            % val is the genomic location of the probePair.
                            val = ceil(probeset1(i).probe_location);
                        else
                            probeset1(i).probe_assignment   = 4;
                            probeset1(i+1).probe_assignment = 4;
                        end;
                    end;
                end;
                
                % basic plot parameters.
                left            = 0.15;
                height          = 0.5/8;
                base            = 0.1;
                vertical_margin = 0.3/8;
                TickSize        = -0.005;  %negative for outside, percentage of longest chr figure.
                
                % define Y-max, for scaling plots.
                maxY = 2;
                
                %define colors for colorBars plot
                colorNoData = [1.0   1.0   1.0  ]; %used when no data is available for the bin.
                colorInit   = [0.5   0.5   0.5  ]; %external; used in blending at ends of chr.
                %define colors for angleplot
                [colorA,colorB, colorAB, colorAAB,colorABB, colorAAAB,colorABBB, colorAAAAB,colorAAABB,...
                    colorAABBB,colorABBBB, colorAAAAAB,colorABBBBB, colorPeak,colorCutoff ] = DefineColors();
                if (show_unnassigned == true) || (show_uncalibrated == true)
                    colorUn_Hom  = [0.5   0.5   1.0  ]; % homozygous unassigned.
                else
                    colorUn_Hom  = [1.0   1.0   1.0  ]; % homozygous unassigned.
                end;
                
                width  = Chr_max_width*chr_size(Zoom_chr)/max(chr_size);
                bottom = base + (8-dataset)*(height+vertical_margin);
                subplot('Position',[left bottom width height]);
                
                data0 = [];   data0_ = [];
                data1 = [];   data1_ = [];
                data2 = [];   data2_ = [];
                for i = 1:SNP_probeset_length
                    if (probeset1(i).probe_chromosome == Zoom_chr)
                        if (probeset1(i).probe_assignment == 0)
                            data0  = [data0 probeset1(i).probe_assignment];
                            data0_ = [data0_ probeset1(i).probe_location];
                        elseif (probeset1(i).probe_assignment == 1)
                            data1  = [data1 probeset1(i).probe_assignment];
                            data1_ = [data1_ probeset1(i).probe_location];
                        elseif (probeset1(i).probe_assignment == 2)
                            data2  = [data2 probeset1(i).probe_assignment];
                            data2_ = [data2_ probeset1(i).probe_location];
                        end;
                    end;
                end;
                % 0 : heterozygous.
                % 1 : homozygous, "b:cyan".
                % 2 : homozygous, "a:magenta".
                plot(data0_/bases_per_bin,data0,'k.');
                hold on;
                plot(data1_/bases_per_bin,data1,'m.');
                plot(data2_/bases_per_bin,data2,'c.');
                hold off;
                
                xlim([0,chr_size(Zoom_chr)*chr_length_scale_multiplier]);
                ylim([0,maxY]);
                set(gca,'YTick',[]);
                set(gca,'TickLength',[(TickSize*chr_size(1)/chr_size(Zoom_chr)) 0]); %ensures same tick size on all subfigs.
                ylabel(names(dataset), 'Rotation', 0, 'HorizontalAlign', 'right');
                set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
                set(gca,'XTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8','2.0','2.2','2.4','2.6','2.8','3.0','3.2'});
                
                set(gca,'FontSize',6);
                if (dataset == 1)
                    title(['Zoom View   (chr=' num2str(Zoom_chr) '; range=[' num2str(Zoom_range(1)) '%..' num2str(Zoom_range(2)) '%])'],'Interpreter','none','FontSize',12);
                end;
            end;
            
            % save then delete figures.
            if ispc  % Windows
                fig_dir = 'figures\zoom\';
                if (isdir([file_dir 'figures\zoom']) == 0)
                    mkdir([file_dir 'figures\zoom']);
                end;
            else     % MacOS
                fig_dir = 'figures/zoom/';
                if (isdir([file_dir 'figures/zoom']) == 0)
                    mkdir([file_dir 'figures/zoom']);
                end;
            end;
            %saveas(fig, [file_dir fig_dir 'chr-' num2str(Zoom_chr) '_range-' num2str(Zoom_range(1)) '..' num2str(Zoom_range(2)) '.eps'], 'epsc');
            saveas(fig, [file_dir fig_dir 'chr-' num2str(Zoom_chr) '_range-' num2str(Zoom_range(1)) '..' num2str(Zoom_range(2)) '.png'], 'png');
            delete(fig);
            
            % clean up variables not needed further.
        end;
    end;
    
    %% ========================================================================
    % Graphing distribution of usable probe data across genome.
    %==========================================================================
    if (Display_distribution == 1)
        % define number of and labels for chromosomes.
        chr_labels = {'Chr1','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7','ChrR'};
        num_chr = length(chr_labels);
        
        bases_per_bin = SNP_bases_per_bin;
        chr_length_scale_multiplier = 1/bases_per_bin;
        
        for dataset = 1:length(names);
            chr1_SNPdata = zeros(1,ceil(chr_size(1)/bases_per_bin));
            chr2_SNPdata = zeros(1,ceil(chr_size(2)/bases_per_bin));
            chr3_SNPdata = zeros(1,ceil(chr_size(3)/bases_per_bin));
            chr4_SNPdata = zeros(1,ceil(chr_size(4)/bases_per_bin));
            chr5_SNPdata = zeros(1,ceil(chr_size(5)/bases_per_bin));
            chr6_SNPdata = zeros(1,ceil(chr_size(6)/bases_per_bin));
            chr7_SNPdata = zeros(1,ceil(chr_size(7)/bases_per_bin));
            chr8_SNPdata = zeros(1,ceil(chr_size(8)/bases_per_bin));
            
            chr1_CGHdata = zeros(1,ceil(chr_size(1)/bases_per_bin));
            chr2_CGHdata = zeros(1,ceil(chr_size(2)/bases_per_bin));
            chr3_CGHdata = zeros(1,ceil(chr_size(3)/bases_per_bin));
            chr4_CGHdata = zeros(1,ceil(chr_size(4)/bases_per_bin));
            chr5_CGHdata = zeros(1,ceil(chr_size(5)/bases_per_bin));
            chr6_CGHdata = zeros(1,ceil(chr_size(6)/bases_per_bin));
            chr7_CGHdata = zeros(1,ceil(chr_size(7)/bases_per_bin));
            chr8_CGHdata = zeros(1,ceil(chr_size(8)/bases_per_bin));
            
            h = waitbar(0.0,'Please wait...','CreateCancelBtn','delete(h); clear h');
            set(h,'Name','Graphing genomic distribution');
            counter = 0;
            for i = 1:SNP_probeset_length
                if (mod(i,2) == 1)
                    counter = counter+1;
                    if (counter == 500)
                        waitbar(i/SNP_probeset_length,h,[num2str(i) '/' num2str(SNP_probeset_length)]);
                        counter = 0;
                    end;
                    
                    % determines if probe pair is useful; both probes have data.
                    if (isnan(probeset1(i).probe_Ratio(dataset)) == 1) && (isnan(probeset1(i+1).probe_Ratio(dataset)) == 1)
                        % val is the genomic location of the probePair.
                        val = ceil(probeset1(i).probe_location/bases_per_bin);
                        
                        % Generates histograms of usable Probe pairs per chromosome.
                        if (probeset1(i).probe_chromosome == 1)
                            chr1_SNPdata(val) = chr1_SNPdata(val)+1;
                        elseif (probeset1(i).probe_chromosome == 2)
                            chr2_SNPdata(val) = chr2_SNPdata(val)+1;
                        elseif (probeset1(i).probe_chromosome == 3)
                            chr3_SNPdata(val) = chr3_SNPdata(val)+1;
                        elseif (probeset1(i).probe_chromosome == 4)
                            chr4_SNPdata(val) = chr4_SNPdata(val)+1;
                        elseif (probeset1(i).probe_chromosome == 5)
                            chr5_SNPdata(val) = chr5_SNPdata(val)+1;
                        elseif (probeset1(i).probe_chromosome == 6)
                            chr6_SNPdata(val) = chr6_SNPdata(val)+1;
                        elseif (probeset1(i).probe_chromosome == 7)
                            chr7_SNPdata(val) = chr7_SNPdata(val)+1;
                        elseif (probeset1(i).probe_chromosome == 8)
                            chr8_SNPdata(val) = chr8_SNPdata(val)+1;
                        end;
                    end;
                end;
            end;
            if exist('h','var'); delete(h); clear h; end;
            h = waitbar(0.0,'Please wait...','CreateCancelBtn','delete(h); clear h');
            set(h,'Name','Plotting CGH data.');
            counter = 0;
            for i = 1:CGH_probeset_length
                counter = counter+1;
                if (counter == 500)
                    waitbar(i/CGH_probeset_length,h,[num2str(i) '/' num2str(CGH_probeset_length)]);
                    counter = 0;
                end;
                
                % val is the genomic location of the probePair.
                val  = ceil((probeset2(i).probe_location)/bases_per_bin);
                val2 = probeset2(i).probe_Ratio(dataset);
                
                % Determines distribution of CGH data.
                if (isnan(val2) ~= 1)
                    if (probeset2(i).probe_chromosome == 1)
                        chr1_CGHdata(val) = chr1_CGHdata(val)+1;
                    elseif (probeset2(i).probe_chromosome == 2)
                        chr2_CGHdata(val) = chr2_CGHdata(val)+1;
                    elseif (probeset2(i).probe_chromosome == 3)
                        chr3_CGHdata(val) = chr3_CGHdata(val)+1;
                    elseif (probeset2(i).probe_chromosome == 4)
                        chr4_CGHdata(val) = chr4_CGHdata(val)+1;
                    elseif (probeset2(i).probe_chromosome == 5)
                        chr5_CGHdata(val) = chr5_CGHdata(val)+1;
                    elseif (probeset2(i).probe_chromosome == 6)
                        chr6_CGHdata(val) = chr6_CGHdata(val)+1;
                    elseif (probeset2(i).probe_chromosome == 7)
                        chr7_CGHdata(val) = chr7_CGHdata(val)+1;
                    elseif (probeset2(i).probe_chromosome == 8)
                        chr8_CGHdata(val) = chr8_CGHdata(val)+1;
                    end;
                end;
            end;
            if exist('h','var'); delete(h); clear h; end;
            
            %generates plots from data, log scaled on Yaxis for SNPs display.
            if (logScale == true)
                plot1 = log(chr1_SNPdata+0.01);
                plot2 = log(chr2_SNPdata+0.01);
                plot3 = log(chr3_SNPdata+0.01);
                plot4 = log(chr4_SNPdata+0.01);
                plot5 = log(chr5_SNPdata+0.01);
                plot6 = log(chr6_SNPdata+0.01);
                plot7 = log(chr7_SNPdata+0.01);
                plot8 = log(chr8_SNPdata+0.01);
                plot1_CGH = log(chr1_CGHdata+0.01);
                plot2_CGH = log(chr2_CGHdata+0.01);
                plot3_CGH = log(chr3_CGHdata+0.01);
                plot4_CGH = log(chr4_CGHdata+0.01);
                plot5_CGH = log(chr5_CGHdata+0.01);
                plot6_CGH = log(chr6_CGHdata+0.01);
                plot7_CGH = log(chr7_CGHdata+0.01);
                plot8_CGH = log(chr8_CGHdata+0.01);
            else
                plot1 = chr1_SNPdata;
                plot2 = chr2_SNPdata;
                plot3 = chr3_SNPdata;
                plot4 = chr4_SNPdata;
                plot5 = chr5_SNPdata;
                plot6 = chr6_SNPdata;
                plot7 = chr7_SNPdata;
                plot8 = chr8_SNPdata;
                plot1_CGH = chr1_CGHdata;
                plot2_CGH = chr2_CGHdata;
                plot3_CGH = chr3_CGHdata;
                plot4_CGH = chr4_CGHdata;
                plot5_CGH = chr5_CGHdata;
                plot6_CGH = chr6_CGHdata;
                plot7_CGH = chr7_CGHdata;
                plot8_CGH = chr8_CGHdata;
            end;
            
            % basic plot parameters.
            left            = 0.15;
            height          = 0.5/num_chr;
            base            = 0.1;
            vertical_margin = 0.3/num_chr;
            TickSize        = -0.005;  %negative for outside, percentage of longest chr figure.
            
            %find Y-max, for scaling plots.
            maxY = max([plot1 plot2 plot3 plot4 plot5 plot6 plot7 plot8]);
            
            figure;
            set(gcf, 'Position', [0 70 1024 600]);
            for chromosome = [8 1:7]
                if (chromosome == 8)
                    usedPlot1 = plot8;
                    usedPlot2 = plot8_CGH;
                elseif (chromosome == 1)
                    usedPlot1 = plot1;
                    usedPlot2 = plot1_CGH;
                elseif (chromosome == 2)
                    usedPlot1 = plot2;
                    usedPlot2 = plot2_CGH;
                elseif (chromosome == 3)
                    usedPlot1 = plot3;
                    usedPlot2 = plot3_CGH;
                elseif (chromosome == 4)
                    usedPlot1 = plot4;
                    usedPlot2 = plot4_CGH;
                elseif (chromosome == 5)
                    usedPlot1 = plot5;
                    usedPlot2 = plot5_CGH;
                elseif (chromosome == 6)
                    usedPlot1 = plot6;
                    usedPlot2 = plot6_CGH;
                elseif (chromosome == 7)
                    usedPlot1 = plot7;
                    usedPlot2 = plot7_CGH;
                end;
                
                width  = Chr_max_width*chr_size(chromosome)/max(chr_size);
                if (chromosome == 8)
                    bottom = base + (chromosome)*(height+vertical_margin);
                else
                    bottom = base + (8-chromosome)*(height+vertical_margin);
                end;
                subplot('Position',[left bottom width height]);
                %show centromere.
                x1 = cen_start(chromosome)*chr_length_scale_multiplier;
                x2 = cen_end(chromosome)*chr_length_scale_multiplier;
                plot([x1; x1],[0;maxY],'g');
                hold on;
                plot([x2; x2],[0;maxY],'g');
                stairs(usedPlot1,'color',[1 0 0]);
                stairs(usedPlot2,'color',[0 0 0]);
                hold off;
                
                xlim([0,chr_size(chromosome)*chr_length_scale_multiplier]);
                ylim([0,maxY]);
                set(gca,'TickLength',[(TickSize*chr_size(1)/chr_size(chromosome)) 0]); %ensures saame tick size on all subfigs.
                ylabel(chr_labels(chromosome), 'Rotation', 90, 'HorizontalAlign', 'center', 'VerticalAlign', 'bottom');
                set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
                set(gca,'XTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8','2.0','2.2','2.4','2.6','2.8','3.0','3.2'});
                if (chromosome == 8)
                    title('Distribution of CGH/SNP probes with usable data.');
                end;
                set(gca,'FontSize',8);
            end;
            
            % Dataset label in top right of figure.
            subplot('Position',[0.8 0.9 0.01 0.01]);
            text(0,5,['Array ' names{dataset}]);
            axis off;
            % Main figure axes key.
            subplot('Position',[0.775 0.6 0.1 0.1]);
            axis on square;
            xlim([0,1]);
            ylim([0,1]);
            set(gca,'FontSize',8);
            xlabel(['Mb (' num2Str(bases_per_bin) ' bases/bin.)']);
            ylabel('Probe count');
            set(gca,'XTick',[0 0.25 0.5 0.75 1]);
            set(gca,'YTick',[0 0.25 0.5 0.75 1]);
            set(gca,'TickLength',[0.05 0]);
            set(gca,'XTickLabel','');
            set(gca,'YTickLabel','');
            % Main figure colors key.
            subplot('Position',[0.7 0.3 0.2 0.2]);
            axis off square;
            xlim([0,1]);
            ylim([-0.1,1]);
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
            patch([0 0.2 0.2 0], [0.8 0.8 0.9 0.9], [0   0   0  ]);
            text(0.3,0.85,'CGH probe distribution');
            patch([0 0.2 0.2 0], [0.6 0.6 0.7 0.7], [1   0   0  ]);
            text(0.3,0.65,'SNP probe distribution');
            patch([0 0.2 0.2 0], [0.4 0.4 0.5 0.5], [0 1 0]);
            text(0.3,0.45,'Centromere location');
        end;
        
        % clean up variables not needed further.
        clear plot* bases_per_bin num_chr base bottom height i vertical_margin width left maxY val counter x1 x2;
        clear *_SNPdata chr_labels chr_length_scale_multiplier dataset;
    end;
    
    %% ========================================================================
    % Scatterplot of ratios from SNP probe pairs.
    %==========================================================================
    if (Scatter_display == 1)
        fprintf(['\nGenerating scatterplot figures from: ' strrep(file_dir,'\','\\')]);
        %define colors for angleplot
        [colorA,colorB, colorAB, colorAAB,colorABB, colorAAAB,colorABBB, colorAAAAB,colorAAABB,...
            colorAABBB,colorABBBB, colorAAAAAB,colorABBBBB, colorPeak,colorCutoff ] = DefineColors();
        maxY = 10;
        
        for dataset = 1:length(names);
            % Initializes vectors used to hold number of SNPs in each
            bases_per_bin = SNP_bases_per_bin;
            %    interpretation catagory.
            for i = 1:8   % eight chromosomes.
                for j = 1:14   % 14 SNP interpretation catagories tracked.
                    chr_SNPdata{i,j} = zeros(1,ceil(chr_size(i)/bases_per_bin));
                end;
                for j = 1:2   % two CGH data catagories tracked.
                    chr_CGHdata{i,j} = zeros(1,ceil(chr_size(i)/bases_per_bin));
                end;
            end;
            h = waitbar(0.0,'Please wait...','CreateCancelBtn','delete(h); clear h');
            set(h,'Name','Plotting CGH data.');
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
            for i = 1:8
                for j = 1:length(chr_CGHdata{i,2})
                    if (chr_CGHdata{i,1}(j) == 0)
                        chr_CGHdata{i,2}(j) = 1;
                    else
                        chr_CGHdata{i,2}(j) = chr_CGHdata{i,2}(j)/chr_CGHdata{i,1}(j);
                    end;
                end;
            end;
            % precalculation of chromosome copy numbers. dddd
            [chr_breaks, chrCopyNum] = ...
                FindChrSizes(segmental_aneuploidy,CGH_probeset_length,probeset2,dataset,chr_size,flow_ploidy);
            % Find homozygous peaks, by fitting to all data from whole experiment.
            [realHomozygous_peak, fit] = ...
                FindRealHomozygousPeaks_2(chrCopyNum,SNP_probeset_length,probeset1,dataset,chr_breaks,chr_size,show_unnassigned,DataTypeToUse,show_fitting);
            %% Definitions of peak locations.
            [monosomy_peak,disomy_peak,trisomy_peak,tetrasomy_peak,pentasomy_peak,hexasomy_peak ] = ...
                FindPeaks(realHomozygous_peak);
            %% Definitions of initial cutoff locations.
            [monosomy_cutoff,disomy_cutoff,trisomy_cutoff,tetrasomy_cutoff,pentasomy_cutoff,hexasomy_cutoff ] = ...
                FindCutoffs(monosomy_peak,disomy_peak,trisomy_peak,tetrasomy_peak,pentasomy_peak,hexasomy_peak);
            
            %%=================================================================
            % Initializes SNP data holders.
            data1_0 = [];    data2_0 = [];
            data1_1 = [];    data2_1 = [];
            data1_2 = [];    data2_2 = [];
            % Gathers all SNP data for this dataset.
            for i = 1:2:SNP_probeset_length
                % determines if probe pair is useful; both probes have data.
                if (isnan(probeset1(i).probe_Ratio(dataset)) == 0) && (isnan(probeset1(i+1).probe_Ratio(dataset)) == 0)
                    % swaps probes if polarity switch indicates to do so.
                    val1 = probeset1(i).probe_Ratio(dataset);
                    val2 = probeset1(i+1).probe_Ratio(dataset);
                    if (probeset1(i).probe_polarity == 1)       % [no-swap]
                        data1_1(length(data1_1)+1) = val1;
                        data2_1(length(data2_1)+1) = val2;
                    elseif (probeset1(i).probe_polarity == 2)   % [swap]
                        data1_2(length(data1_2)+1) = val2;
                        data2_2(length(data2_2)+1) = val1;
                    elseif (probeset1(i).probe_polarity == 0)   % unassigned probe pairs.   [no-swap]
                        data1_0(length(data1_0)+1) = val1;
                        data2_0(length(data2_0)+1) = val2;
                    else
                        % null action for when (probe_polarity == 4)
                        % due to probe design error; probes are identical.
                    end;
                end;
            end;
            %%=================================================================
            %% Precalculation of distance cutoff.
            % Calculates SNP magnitudes for determining distance cutoffs.
            for iii = 1;
                histAll_2a = [];
                for i = 1:length(data1_0)
                    histAll_2a = [histAll_2a sqrt(data1_0(i)^2+data2_0(i)^2)];
                end;
                for i = 1:length(data1_1)
                    histAll_2a = [histAll_2a sqrt(data1_1(i)^2+data2_1(i)^2)];
                end;
                for i = 1:length(data1_2)
                    histAll_2a = [histAll_2a sqrt(data1_2(i)^2+data2_2(i)^2)];
                end;
                % Places bounds on SNP magnitude histogram.
                histAll_2a(histAll_2a==0) = [];
                histAll_2a(histAll_2a>sqrt(32)) = [];
                histAll_2a(length(histAll_2a)+1) = 0;
                histAll_2a(length(histAll_2a)+1) = sqrt(32);
                % Smooths histogram.
                smoothed_2a = smooth_gaussian(hist(histAll_2a,200),5,20);
                % Makes a smoothed version of just the endpoints used to ensure histogram bounds.
                histAll_2b = [];
                histAll_2b(1) = 0;
                histAll_2b(2) = sqrt(32);
                smoothed_2b = smooth_gaussian(hist(histAll_2b,200),5,20);
                % Subtracts the smoothed endpoints from the histogram to remove the influence of the added endpoints.
                smoothed_2a = (smoothed_2a-smoothed_2b);
                % Scales the histogram to the histogram max.
                smoothed_2a = smoothed_2a/max(smoothed_2a);
                % Finds the location of the histogram peak.
                histMax = find(smoothed_2a == max(smoothed_2a));
                % Defines magnitude cutoff from SNP magnitude histogram peak location.
                realDistCutoff = histMax/200*sqrt(32)/2;
            end;
            %%=================================================================
            %% Precalculation of allelic-fraction cutoffs.
            % collect SNP data into raw and smoothed histograms.
            histAll_1a = [];
            if (distanceCutoff == 1)
                for i = 1:length(data1_0)
                    if (sqrt((data2_0(i))^2+(data1_0(i))^2) > realDistCutoff) && (sqrt((data2_0(i))^2+(data1_0(i))^2) < realDistCutoff*4)
                        histAll_1a = [histAll_1a data2_0(i)/(data2_0(i)+data1_0(i))];
                    end;
                end;
                for i = 1:length(data1_1)
                    if (sqrt((data2_1(i))^2+(data1_1(i))^2) > realDistCutoff) && (sqrt((data2_1(i))^2+(data1_1(i))^2) < realDistCutoff*4)
                        histAll_1a = [histAll_1a data2_1(i)/(data2_1(i)+data1_1(i))];
                    end;
                end;
                for i = 1:length(data1_2)
                    if (sqrt((data2_2(i))^2+(data1_2(i))^2) > realDistCutoff) && (sqrt((data2_2(i))^2+(data1_2(i))^2) < realDistCutoff*4)
                        histAll_1a = [histAll_1a data2_2(i)/(data2_2(i)+data1_2(i))];
                    end;
                end;
            else
                for i = 1:length(data1_0)
                    histAll_1a = [histAll_1a data2_0(i)/(data2_0(i)+data1_0(i))];
                end;
                for i = 1:length(data1_1)
                    histAll_1a = [histAll_1a data2_1(i)/(data2_1(i)+data1_1(i))];
                end;
                for i = 1:length(data1_2)
                    histAll_1a = [histAll_1a data2_2(i)/(data2_2(i)+data1_2(i))];
                end;
            end;
            histAll_1a(histAll_1a==0)        = [];
            histAll_1a(histAll_1a>1)         = [];
            histAll_1a(length(histAll_1a)+1) = 0;
            histAll_1a(length(histAll_1a)+1) = 1;
            % make a histogram of SNP allele ratio data in theoretical
            % heterozygous regions (inner quartiles only).
            hist_a   = hist(histAll_1a,200);
            hist_b   = fliplr(hist(histAll_1a,200));
            % make a summary histogram which mirrors the data.
            raw      = hist_a+hist_b;
            smoothed = smooth_gaussian(raw,2,20);
            
            fit_curve_1 = fit(1)*exp(-0.5*(((1:200)-fit(2))./fit(3)).^2);
            fit_curve_2 = fit(4)*exp(-0.5*(((1:200)-fit(5))./fit(6)).^2);
            fit_curve_3 = fit(7)*exp(-0.5*(((1:200)-fit(8))./fit(9)).^2);
            fit_curve_tot = (fit_curve_1+fit_curve_2+fit_curve_3);
            fig = figure;
            hold on;
            % SD of het peak.
            area(max([1 round(fit(2)-fit(3))]):min([200 round(fit(2)+fit(3))]),fit_curve_1(max([1 round(fit(2)-fit(3))]):min([200 round(fit(2)+fit(3))])),'FaceColor',[1.0 0.50 0.50]);
            % SD of left hom peak.
            area(max([1 round(fit(5)-fit(6))]):min([200 round(fit(5)+fit(6))]),fit_curve_2(max([1 round(fit(5)-fit(6))]):min([200 round(fit(5)+fit(6))])),'FaceColor',[0.50 1.0 0.50]);
            % SD of right hom peak.
            area(max([1 round(fit(8)-fit(9))]):min([200 round(fit(8)+fit(9))]),fit_curve_3(max([1 round(fit(8)-fit(9))]):min([200 round(fit(8)+fit(9))])),'FaceColor',[0.50 0.50 1.0]);
            plot(fit_curve_tot);
            plot(raw      ,'color',[0.5 0.5 1.0],'linestyle','-','linewidth',1);
            plot(smoothed ,'color',[0.0 0.0 1.0],'linestyle','-','linewidth',2);
            text(50,max(raw),names{dataset},'HorizontalAlign','center','VerticalAlign','middle');
            hold off;
            xlim([1,200]);
            ylim([0,max(raw)]);
            
            % save then delete figures of whole dataset allelic fraction histogram and fit curves.
            if ispc  % Windows
                fig_dir = 'figures\scatterHist_perDataset\';
                if (isdir([file_dir 'figures\scatterHist_perDataset']) == 0)
                    mkdir([file_dir 'figures\scatterHist_perDataset']);
                end;
            else     % MacOS
                fig_dir = 'figures/scatterHist_perDataset/';
                if (isdir([file_dir 'figures/scatterHist_perDataset']) == 0)
                    mkdir([file_dir 'figures/scatterHist_perDataset']);
                end;
            end;
            %saveas(fig, [file_dir fig_dir strrep(names{dataset},' ','_') '.eps'], 'epsc');
            saveas(fig, [file_dir fig_dir strrep(names{dataset},' ','_') '.png'], 'png');
            delete(fig);

            %% Generates scatter plot figures for chromosome segments.
            for analyzeChr = chromosomes_to_analyze;
                for segment = 1:length(datasetDetails{dataset}.chrCopyNum{analyzeChr})
                    for i = 1
                        if (main_scatter == true)
                            % main scatterplot figure.
                            fig = figure;
                            scatter(data2_0,data1_0,5,'k');
                            hold on;
                            scatter(data2_1,data1_1,5,'k');
                            scatter(data2_2,data1_2,5,'k');
                            ylabel('SNP_a  Ratio', 'Rotation', 90,'HorizontalAlign', 'center', 'VerticalAlign', 'middle');
                            xlabel('SNP_b  Ratio', 'Rotation', 0, 'HorizontalAlign', 'center', 'VerticalAlign', 'middle');
                            if (analyzeChr == 8)
                                title(['ChrR SNP probe Ratios : ' names{dataset}]);
                            else
                                title(['Chr' num2str(analyzeChr) ' SNP probe Ratios : ' names{dataset}]);
                            end;
                            axis([0,4,0,4],'square');
                            
                            % Scatter plot region labels.
                            if (datasetDetails{dataset}.chrCopyNum{analyzeChr}(1) == 1)
                                text(monosomy_ratio_peak(1)*0.9,4*0.9,'(a)','rotation',-monosomy_peak(1),'horizontalalignment','center','verticalalignment','middle','color',colorA);
                                text(4*0.9,monosomy_ratio_peak(1)*0.9,'(b)','rotation',-monosomy_peak(2),'horizontalalignment','center','verticalalignment','middle','color',colorB);
                            elseif (datasetDetails{dataset}.chrCopyNum{analyzeChr}(1) == 2)
                                text(disomy_ratio_peak(1)*0.9,4*0.9,'(aa)','rotation',-disomy_peak(1),'horizontalalignment','center','verticalalignment','middle','color',colorA);
                                text(4*0.9,4*0.9                   ,'(ab)','rotation',-disomy_peak(2),'horizontalalignment','center','verticalalignment','middle','color',colorAB);
                                text(4*0.9,disomy_ratio_peak(1)*0.9,'(bb)','rotation',-disomy_peak(3),'horizontalalignment','center','verticalalignment','middle','color',colorB);
                            elseif (datasetDetails{dataset}.chrCopyNum{analyzeChr}(1) == 3)
                                text(trisomy_ratio_peak(1)*0.9,4*0.9,'(aaa)','rotation',-trisomy_peak(1),'horizontalalignment','center','verticalalignment','middle','color',colorA);
                                text(trisomy_ratio_peak(2)*0.9,4*0.9,'(aab)','rotation',-trisomy_peak(2),'horizontalalignment','center','verticalalignment','middle','color',colorAAB);
                                text(4*0.9,trisomy_ratio_peak(2)*0.9,'(abb)','rotation',-trisomy_peak(3),'horizontalalignment','center','verticalalignment','middle','color',colorABB);
                                text(4*0.9,trisomy_ratio_peak(1)*0.9,'(bbb)','rotation',-trisomy_peak(4),'horizontalalignment','center','verticalalignment','middle','color',colorB);
                            elseif (datasetDetails{dataset}.chrCopyNum{analyzeChr}(1) == 4)
                                text(tetrasomy_ratio_peak(1)*0.9,4*0.9,'(aaaa)','rotation',-tetrasomy_peak(1),'horizontalalignment','center','verticalalignment','middle','color',colorA);
                                text(tetrasomy_ratio_peak(2)*0.9,4*0.9,'(aaab)','rotation',-tetrasomy_peak(2),'horizontalalignment','center','verticalalignment','middle','color',colorAAAB);
                                text(4*0.9,4*0.9                      ,'(aabb)','rotation',-tetrasomy_peak(3),'horizontalalignment','center','verticalalignment','middle','color',colorAB);
                                text(4*0.9,tetrasomy_ratio_peak(2)*0.9,'(abbb)','rotation',-tetrasomy_peak(4),'horizontalalignment','center','verticalalignment','middle','color',colorABBB);
                                text(4*0.9,tetrasomy_ratio_peak(1)*0.9,'(bbbb)','rotation',-tetrasomy_peak(5),'horizontalalignment','center','verticalalignment','middle','color',colorB);
                            elseif (datasetDetails{dataset}.chrCopyNum{analyzeChr}(1) == 5)
                                text(pentasomy_ratio_peak(1)*0.9,4*0.9,'(aaaaa)','rotation',-pentasomy_peak(1),'horizontalalignment','center','verticalalignment','middle','color',colorA);
                                text(pentasomy_ratio_peak(2)*0.9,4*0.9,'(aaaab)','rotation',-pentasomy_peak(2),'horizontalalignment','center','verticalalignment','middle','color',colorAAAAB);
                                text(pentasomy_ratio_peak(3)*0.9,4*0.9,'(aaabb)','rotation',-pentasomy_peak(3),'horizontalalignment','center','verticalalignment','middle','color',colorAAABB);
                                text(4*0.9,pentasomy_ratio_peak(3)*0.9,'(aabbb)','rotation',-pentasomy_peak(4),'horizontalalignment','center','verticalalignment','middle','color',colorAABBB);
                                text(4*0.9,pentasomy_ratio_peak(2)*0.9,'(abbbb)','rotation',-pentasomy_peak(5),'horizontalalignment','center','verticalalignment','middle','color',colorABBBB);
                                text(4*0.9,pentasomy_ratio_peak(1)*0.9,'(bbbbb)','rotation',-pentasomy_peak(6),'horizontalalignment','center','verticalalignment','middle','color',colorB);
                            elseif (datasetDetails{dataset}.chrCopyNum{analyzeChr}(1) == 6)
                                text(hexasomy_ratio_peak(1)*0.9,4*0.9,'(aaaaaa)','rotation',-hexasomy_peak(1),'horizontalalignment','center','verticalalignment','middle','color',colorA);
                                text(hexasomy_ratio_peak(2)*0.9,4*0.9,'(aaaaab)','rotation',-hexasomy_peak(2),'horizontalalignment','center','verticalalignment','middle','color',colorAAAAAB);
                                text(hexasomy_ratio_peak(3)*0.9,4*0.9,'(aaaabb)','rotation',-hexasomy_peak(3),'horizontalalignment','center','verticalalignment','middle','color',colorAAB);
                                text(4*0.9,4*0.9                     ,'(aaabbb)','rotation',-hexasomy_peak(4),'horizontalalignment','center','verticalalignment','middle','color',colorAB);
                                text(4*0.9,hexasomy_ratio_peak(3)*0.9,'(aabbbb)','rotation',-hexasomy_peak(5),'horizontalalignment','center','verticalalignment','middle','color',colorABB);
                                text(4*0.9,hexasomy_ratio_peak(2)*0.9,'(abbbbb)','rotation',-hexasomy_peak(6),'horizontalalignment','center','verticalalignment','middle','color',colorABBBBB);
                                text(4*0.9,hexasomy_ratio_peak(1)*0.9,'(bbbbbb)','rotation',-hexasomy_peak(7),'horizontalalignment','center','verticalalignment','middle','color',colorB);
                            else % datasetDetails{dataset}.chrCopyNum{analyzeChr}(1) == [1,2,4,etc]
                                text(disomy_ratio_peak(1)*0.9,4*0.9,'(aa)','rotation',-disomy_peak(1),'horizontalalignment','center','verticalalignment','middle','color',colorA);
                                text(4*0.9,4*0.9                   ,'(ab)','rotation',-disomy_peak(2),'horizontalalignment','center','verticalalignment','middle','color',colorAB);
                                text(4*0.9,disomy_ratio_peak(1)*0.9,'(bb)','rotation',-disomy_peak(3),'horizontalalignment','center','verticalalignment','middle','color',colorB);
                            end;
                            hold on;
                            
                            % Scatter plot lines.
                            peaks   = [];
                            cutoffs = [];
                            if (datasetDetails{dataset}.chrCopyNum{analyzeChr}(1) == 1)
                                peaks   = monosomy_peak;
                                cutoffs = monosomy_cutoff;
                                count   = 1;
                            elseif (datasetDetails{dataset}.chrCopyNum{analyzeChr}(1) == 2)
                                peaks   = disomy_peak;
                                cutoffs = disomy_cutoff;
                                count   = 2;
                            elseif (datasetDetails{dataset}.chrCopyNum{analyzeChr}(1) == 3)
                                peaks   = trisomy_peak;
                                cutoffs = trisomy_cutoff;
                                count   = 3;
                            elseif (datasetDetails{dataset}.chrCopyNum{analyzeChr}(1) == 4)
                                peaks   = tetrasomy_peak;
                                cutoffs = tetrasomy_cutoff;
                                count   = 4;
                            elseif (datasetDetails{dataset}.chrCopyNum{analyzeChr}(1) == 5)
                                peaks   = pentasomy_peak;
                                cutoffs = pentasomy_cutoff;
                                count   = 5;
                            elseif (datasetDetails{dataset}.chrCopyNum{analyzeChr}(1) == 6)
                                peaks   = hexasomy_peak;
                                cutoffs = hexasomy_cutoff;
                                count   = 6;
                            else
                                peaks   = disomy_peak;
                                cutoffs = disomy_cutoff;
                                count   = 2;
                            end;
                            
                            % draw peak lines, using rotation.
                            for i = 1:count+1
                                center = [0,0,0];
                                Point = [0,7,0];
                                T1 = makehgtform('translate',center);
                                T_rot = makehgtform('zrotate',-peaks(i)/180*pi);
                                NewPoint = T1 * T_rot * inv(T1) * [Point,1]';
                                NewPoint = NewPoint(1:2)';
                                plot([0 NewPoint(1)], [0 NewPoint(2)],'color',colorPeak,'linestyle','-','linewidth',2);
                            end;
                            % draw cutoff lines, using rotation.
                            for i = 1:count
                                if (use_cutoff(i) == true)
                                    center = [0,0,0];
                                    Point = [0,7,0];
                                    T1 = makehgtform('translate',center);
                                    T_rot = makehgtform('zrotate',-cutoffs(i)/180*pi);
                                    NewPoint = T1 * T_rot * inv(T1) * [Point,1]';
                                    NewPoint = NewPoint(1:2)';
                                    plot([0 NewPoint(1)], [0 NewPoint(2)],'color',colorCutoff,'linestyle','-','linewidth',2);
                                else
                                    center = [0,0,0];
                                    Point = [0,7,0];
                                    T1 = makehgtform('translate',center);
                                    T_rot = makehgtform('zrotate',-cutoffs(i)/180*pi);
                                    NewPoint = T1 * T_rot * inv(T1) * [Point,1]';
                                    NewPoint = NewPoint(1:2)';
                                    plot([0 NewPoint(1)], [0 NewPoint(2)],'color',colorCutoff,'linestyle',':','linewidth',2);
                                end;
                            end;
                            
                            % distance cutoff.
                            if (distanceCutoff == 1)
                                theta = 1:5:95;
                                r = realDistCutoff;
                                plot(r*cosd(theta),r*sind(theta),'color',colorPeak,'linestyle','-','linewidth',1);
                                plot(r*4*cosd(theta),r*4*sind(theta),'color',colorPeak,'linestyle','-','linewidth',1);
                            end;
                            hold off;
                            
                            % angleplot.
                            width           = 0.1;
                            bottom          = 0.5;
                            height          = 0.5/6;
                            subplot('Position',[0.03 bottom width height]);
                            %define colors for angleplot
                            [colorA,colorB, colorAB, colorAAB,colorABB, colorAAAB,colorABBB, colorAAAAB,colorAAABB,...
                                colorAABBB,colorABBBB, colorAAAAAB,colorABBBBB, colorPeak,colorCutoff ] = DefineColors();
                            hold on;
                            % Angleplot lines.
                            if (datasetDetails{dataset}.chrCopyNum{analyzeChr}(1) == 1)
                                area(1                                :round(monosomy_cutoff(1)/22.5*50),smoothed(1                                :round(monosomy_cutoff(1)/22.5*50)),'FaceColor',colorA, 'EdgeColor',colorA);
                                area(round(monosomy_cutoff(1)/22.5*50):200                              ,smoothed(round(monosomy_cutoff(1)/22.5*50):200                              ),'FaceColor',colorB, 'EdgeColor',colorB);
                                for i = 1:2
                                    plot([monosomy_peak(i)/22.5*50;   monosomy_peak(i)/22.5*50  ] ,[0; 1],'color',colorPeak);
                                end;
                                for i = 1
                                    plot([monosomy_cutoff(i)/22.5*50; monosomy_cutoff(i)/22.5*50] ,[0; 1],'color',colorCutoff);
                                end;
                                set(gca,'XTick',0:100:200);
                                set(gca,'XTickLabel',{'a','ab','b'});
                            elseif (datasetDetails{dataset}.chrCopyNum{analyzeChr}(1) == 2)
                                if (disomy_cutoff(1) < 0)
                                    area(1                              :round(disomy_cutoff(2)/22.5*50),smoothed(1                              :round(disomy_cutoff(2)/22.5*50)),'FaceColor',colorA, 'EdgeColor',colorAB);
                                    area(round(disomy_cutoff(2)/22.5*50):200                            ,smoothed(round(disomy_cutoff(2)/22.5*50):200                            ),'FaceColor',colorB, 'EdgeColor',colorB);
                                else
                                    area(1                              :round(disomy_cutoff(1)/22.5*50),smoothed(1                              :round(disomy_cutoff(1)/22.5*50)),'FaceColor',colorA, 'EdgeColor',colorA);
                                    area(round(disomy_cutoff(1)/22.5*50):round(disomy_cutoff(2)/22.5*50),smoothed(round(disomy_cutoff(1)/22.5*50):round(disomy_cutoff(2)/22.5*50)),'FaceColor',colorAB,'EdgeColor',colorAB);
                                    area(round(disomy_cutoff(2)/22.5*50):200                            ,smoothed(round(disomy_cutoff(2)/22.5*50):200                            ),'FaceColor',colorB, 'EdgeColor',colorB);
                                end;
                                for i = 1:3
                                    plot([disomy_peak(i)/22.5*50;   disomy_peak(i)/22.5*50  ] ,[0; 1],'color',colorPeak);
                                end;
                                for i = 1:2
                                    if (disomy_cutoff(i) > 0)
                                        plot([disomy_cutoff(i)/22.5*50; disomy_cutoff(i)/22.5*50] ,[0; 1],'color',colorCutoff);
                                    end;
                                end;
                                set(gca,'XTick',0:100:200);
                                set(gca,'XTickLabel',{'a','ab','b'});
                            elseif (datasetDetails{dataset}.chrCopyNum{analyzeChr}(1) == 3)
                                area(1                               :round(trisomy_cutoff(1)/22.5*50),smoothed(1                               :round(trisomy_cutoff(1)/22.5*50)),'FaceColor',colorA,  'EdgeColor',colorA);
                                area(round(trisomy_cutoff(1)/22.5*50):round(trisomy_cutoff(2)/22.5*50),smoothed(round(trisomy_cutoff(1)/22.5*50):round(trisomy_cutoff(2)/22.5*50)),'FaceColor',colorAAB,'EdgeColor',colorAAB);
                                area(round(trisomy_cutoff(2)/22.5*50):round(trisomy_cutoff(3)/22.5*50),smoothed(round(trisomy_cutoff(2)/22.5*50):round(trisomy_cutoff(3)/22.5*50)),'FaceColor',colorABB,'EdgeColor',colorABB);
                                area(round(trisomy_cutoff(3)/22.5*50):200                             ,smoothed(round(trisomy_cutoff(3)/22.5*50):200                             ),'FaceColor',colorB,  'EdgeColor',colorB);
                                for i = 1:4
                                    plot([trisomy_peak(i)/22.5*50;   trisomy_peak(i)/22.5*50  ] ,[0; 1],'color',colorPeak);
                                end;
                                for i = 1:3
                                    plot([trisomy_cutoff(i)/22.5*50; trisomy_cutoff(i)/22.5*50] ,[0; 1],'color',colorCutoff);
                                end;
                                set(gca,'XTick',[0 round(trisomy_peak(2)/22.5*50) round(trisomy_peak(3)/22.5*50) 200]);
                                set(gca,'XTickLabel',{'a','aab','abb','b'});
                            elseif (datasetDetails{dataset}.chrCopyNum{analyzeChr}(1) == 4)
                                area(1                                 :round(tetrasomy_cutoff(1)/22.5*50),smoothed(1                                 :round(tetrasomy_cutoff(1)/22.5*50)),'FaceColor',colorA,   'EdgeColor',colorA);
                                area(round(tetrasomy_cutoff(1)/22.5*50):round(tetrasomy_cutoff(2)/22.5*50),smoothed(round(tetrasomy_cutoff(1)/22.5*50):round(tetrasomy_cutoff(2)/22.5*50)),'FaceColor',colorAAAB,'EdgeColor',colorAAAB);
                                area(round(tetrasomy_cutoff(2)/22.5*50):round(tetrasomy_cutoff(3)/22.5*50),smoothed(round(tetrasomy_cutoff(2)/22.5*50):round(tetrasomy_cutoff(3)/22.5*50)),'FaceColor',colorAB,  'EdgeColor',colorAB);
                                area(round(tetrasomy_cutoff(3)/22.5*50):round(tetrasomy_cutoff(4)/22.5*50),smoothed(round(tetrasomy_cutoff(3)/22.5*50):round(tetrasomy_cutoff(4)/22.5*50)),'FaceColor',colorABBB,'EdgeColor',colorABBB);
                                area(round(tetrasomy_cutoff(4)/22.5*50):200                               ,smoothed(round(tetrasomy_cutoff(4)/22.5*50):200                               ),'FaceColor',colorB,   'EdgeColor',colorB);
                                for i = 1:5
                                    plot([tetrasomy_peak(i)/22.5*50;   tetrasomy_peak(i)/22.5*50  ] ,[0; 1],'color',colorPeak);
                                end;
                                for i = 1:4
                                    plot([tetrasomy_cutoff(i)/22.5*50; tetrasomy_cutoff(i)/22.5*50] ,[0; 1],'color',colorCutoff);
                                end;
                                set(gca,'XTick',[0 round(trisomy_peak(2)/22.5*50) round(trisomy_peak(3)/22.5*50) 200]);
                                set(gca,'XTickLabel',{'a','aab','abb','b'});
                            elseif (datasetDetails{dataset}.chrCopyNum{analyzeChr}(1) == 5)
                                area(1                                 :round(pentasomy_cutoff(1)/22.5*50),smoothed(1                                 :round(pentasomy_cutoff(1)/22.5*50)),'FaceColor',colorA,    'EdgeColor',colorA);
                                area(round(pentasomy_cutoff(1)/22.5*50):round(pentasomy_cutoff(2)/22.5*50),smootheda(round(pentasomy_cutoff(1)/22.5*50):round(tetrasomy_cutoff(2)/22.5*50)),'FaceColor',colorAAAAB,'EdgeColor',colorAAAAB);
                                area(round(pentasomy_cutoff(2)/22.5*50):round(pentasomy_cutoff(3)/22.5*50),smoothed(round(pentasomy_cutoff(2)/22.5*50):round(tetrasomy_cutoff(3)/22.5*50)),'FaceColor',colorAAABB,'EdgeColor',colorAAABB);
                                area(round(pentasomy_cutoff(3)/22.5*50):round(pentasomy_cutoff(4)/22.5*50),smoothed(round(pentasomy_cutoff(3)/22.5*50):round(tetrasomy_cutoff(4)/22.5*50)),'FaceColor',colorAABBB,'EdgeColor',colorAABBB);
                                area(round(pentasomy_cutoff(4)/22.5*50):round(pentasomy_cutoff(5)/22.5*50),smoothed(round(pentasomy_cutoff(4)/22.5*50):round(tetrasomy_cutoff(5)/22.5*50)),'FaceColor',colorABBBB,'EdgeColor',colorABBBB);
                                area(round(pentasomy_cutoff(5)/22.5*50):200                               ,smoothed(round(pentasomy_cutoff(5)/22.5*50):200                               ),'FaceColor',colorB,    'EdgeColor',colorB);
                                for i = 1:6
                                    plot([pentasomy_peak(i)/22.5*50;   pentasomy_peak(i)/22.5*50  ] ,[0; 1],'color',colorPeak);
                                end;
                                for i = 1:5
                                    plot([pentasomy_cutoff(i)/22.5*50; pentasomy_cutoff(i)/22.5*50] ,[0; 1],'color',colorCutoff);
                                end;
                                set(gca,'XTick',[0 round(trisomy_peak(2)/22.5*50) round(trisomy_peak(3)/22.5*50) 200]);
                                set(gca,'XTickLabel',{'a','aab','abb','b'});
                            elseif (datasetDetails{dataset}.chrCopyNum{analyzeChr}(1) == 6)
                                if (round(hexasomy_cutoff(1)/22.5*50) >= 1) && (round(hexasomy_cutoff(1)/22.5*50) <= 200)
                                    area(1                                :round(hexasomy_cutoff(1)/22.5*50),smoothed(1                                :round(hexasomy_cutoff(1)/22.5*50)),'FaceColor',colorA,     'EdgeColor',colorA);
                                end;
                                if (round(hexasomy_cutoff(1)/22.5*50) >= 1) && (round(hexasomy_cutoff(1)/22.5*50) <= 200) && ...
                                        (round(hexasomy_cutoff(2)/22.5*50) >= 1) && (round(hexasomy_cutoff(2)/22.5*50) <= 200)
                                    area(round(hexasomy_cutoff(1)/22.5*50):round(hexasomy_cutoff(2)/22.5*50),smoothed(round(hexasomy_cutoff(1)/22.5*50):round(hexasomy_cutoff(2)/22.5*50)),'FaceColor',colorAAAAAB,'EdgeColor',colorAAAAAB);
                                end;
                                if (round(hexasomy_cutoff(2)/22.5*50) >= 1) && (round(hexasomy_cutoff(2)/22.5*50) <= 200) && ...
                                        (round(hexasomy_cutoff(3)/22.5*50) >= 1) && (round(hexasomy_cutoff(3)/22.5*50) <= 200)
                                    area(round(hexasomy_cutoff(2)/22.5*50):round(hexasomy_cutoff(3)/22.5*50),smoothed(round(hexasomy_cutoff(2)/22.5*50):round(hexasomy_cutoff(3)/22.5*50)),'FaceColor',colorAAB,   'EdgeColor',colorAAB);
                                end;
                                if (round(hexasomy_cutoff(3)/22.5*50) >= 1) && (round(hexasomy_cutoff(3)/22.5*50) <= 200) && ...
                                        (round(hexasomy_cutoff(4)/22.5*50) >= 1) && (round(hexasomy_cutoff(4)/22.5*50) <= 200)
                                    area(round(hexasomy_cutoff(3)/22.5*50):round(hexasomy_cutoff(4)/22.5*50),smoothed(round(hexasomy_cutoff(3)/22.5*50):round(hexasomy_cutoff(4)/22.5*50)),'FaceColor',colorAB,    'EdgeColor',colorAB);
                                end;
                                if (round(hexasomy_cutoff(4)/22.5*50) >= 1) && (round(hexasomy_cutoff(4)/22.5*50) <= 200) && ...
                                        (round(hexasomy_cutoff(5)/22.5*50) >= 1) && (round(hexasomy_cutoff(5)/22.5*50) <= 200)
                                    area(round(hexasomy_cutoff(4)/22.5*50):round(hexasomy_cutoff(5)/22.5*50),smoothed(round(hexasomy_cutoff(4)/22.5*50):round(hexasomy_cutoff(5)/22.5*50)),'FaceColor',colorABB,   'EdgeColor',colorABB);
                                end;
                                if (round(hexasomy_cutoff(5)/22.5*50) >= 1) && (round(hexasomy_cutoff(5)/22.5*50) <= 200) && ...
                                        (round(hexasomy_cutoff(6)/22.5*50) >= 1) && (round(hexasomy_cutoff(6)/22.5*50) <= 200)
                                    area(round(hexasomy_cutoff(5)/22.5*50):round(hexasomy_cutoff(6)/22.5*50),smoothed(round(hexasomy_cutoff(5)/22.5*50):round(hexasomy_cutoff(6)/22.5*50)),'FaceColor',colorABBBBB,'EdgeColor',colorABBBBB);
                                end;
                                if (round(hexasomy_cutoff(6)/22.5*50) >= 1) && (round(hexasomy_cutoff(6)/22.5*50) <= 200)
                                    area(round(hexasomy_cutoff(6)/22.5*50):200                              ,smoothed(round(hexasomy_cutoff(6)/22.5*50):200                              ),'FaceColor',colorB,     'EdgeColor',colorB);
                                end;
                                for i = 1:7
                                    plot([hexasomy_peak(i)/22.5*50;   hexasomy_peak(i)/22.5*50  ] ,[0; 1],'color',colorPeak);
                                end;
                                for i = 1:6
                                    plot([hexasomy_cutoff(i)/22.5*50; hexasomy_cutoff(i)/22.5*50] ,[0; 1],'color',colorCutoff);
                                end;
                                set(gca,'XTick',[0 round(trisomy_peak(2)/22.5*50) round(trisomy_peak(3)/22.5*50) 200]);
                                set(gca,'XTickLabel',{'a','aab','abb','b'});
                            else % datasetDetails{dataset}.chrCopyNum{analyzeChr}(1) == [4,etc]
                                area(1                              :round(disomy_cutoff(1)/22.5*50),smoothed(1                              :round(disomy_cutoff(1)/22.5*50)),'FaceColor',colorA, 'EdgeColor',colorA);
                                area(round(disomy_cutoff(1)/22.5*50):round(disomy_cutoff(2)/22.5*50),smoothed(round(disomy_cutoff(1)/22.5*50):round(disomy_cutoff(2)/22.5*50)),'FaceColor',colorAB,'EdgeColor',colorAB);
                                area(round(disomy_cutoff(2)/22.5*50):200                            ,smoothed(round(disomy_cutoff(2)/22.5*50):200                            ),'FaceColor',colorB, 'EdgeColor',colorB);
                                for i = 1:3
                                    plot([disomy_peak(i)/22.5*50;   disomy_peak(i)/22.5*50  ] ,[0; 1],'color',colorPeak);
                                end;
                                for i = 1:2
                                    plot([disomy_cutoff(i)/22.5*50; disomy_cutoff(i)/22.5*50] ,[0; 1],'color',colorCutoff);
                                end;
                                set(gca,'XTick',0:100:200);
                                set(gca,'XTickLabel',{'a','ab','b'});
                            end;
                            hold off;
                            set(gca,'YTick',[]);
                            set(gca,'FontSize',8);
                            xlim([0,200]);
                            ylim([0,1]);
                            
                            % distanceplot
                            width           = 0.1;
                            bottom          = 0.25;
                            height          = 0.5/6;
                            subplot('Position',[0.03 bottom width height]);
                            colorBlack  = [0.8 0.8 0.8];
                            hold on;
                            area(1:200,smoothed_2a,'FaceColor',colorBlack,'EdgeColor',colorBlack);
                            
                            plot([realDistCutoff/sqrt(32)*200;   realDistCutoff/sqrt(32)*200  ],[0; 1],'color',colorPeak);
                            plot([realDistCutoff/sqrt(32)*200*2; realDistCutoff/sqrt(32)*200*2],[0; 1],'color',colorPeak);
                            hold off;
                            set(gca,'XTick',[0.000 35.355 70.711 106.066 141.421]);
                            set(gca,'XTickLabel',{'0','1','2','3','4'});
                            set(gca,'YTick',[]);
                            set(gca,'FontSize',8);
                            xlim([0,141.421]);
                            ylim([0,1]);
                            
                            % save then delete figures.
                            if ispc  % Windows
                                fig_dir = 'figures\scatter\';
                                if (isdir([file_dir 'figures\scatter']) == 0)
                                    mkdir([file_dir 'figures\scatter']);
                                end;
                            else     % MacOS
                                fig_dir = 'figures/scatter/';
                                if (isdir([file_dir 'figures/scatter']) == 0)
                                    mkdir([file_dir 'figures/scatter']);
                                end;
                            end;
                            %saveas(fig, [file_dir fig_dir 'dat-' num2str(dataset) '_chr-' num2str(analyzeChr) '_seg-' num2str(segment) '.eps'], 'epsc');
                            saveas(fig, [file_dir fig_dir 'dat-' num2str(dataset) '_chr-' num2str(analyzeChr) '_seg-' num2str(segment) '.png'], 'png');
                            delete(fig);
                        end;
                    end;
                end;
            end;
            % clean up variables not needed further.
            clear plot* bases_per_bin num_chr base bottom height i vertical_margin width left val counter x1 x2;
        end;
    end;
    
    %% ========================================================================
    % Scatterplot of ratios from SNP probe pairs.
    %   Showing probe length as circle radius.
    %==========================================================================
    if (Scatter_probeLength_display == 1)
        for analyzeChr = chromosomes_to_analyze;
            for dataset = 1:length(names);
                % initialize vectors for scattergram analysis
                data1 = [];
                data2 = [];
                data3 = [];
                
                h = waitbar(0.0,'Please wait...','CreateCancelBtn','delete(h); clear h');
                set(h,'Name','Plotting scattergrams of SNP ratios.');
                counter = 0;
                for i = 1:SNP_probeset_length
                    if (mod(i,2) == 1)
                        counter = counter+1;
                        if (counter == 500)
                            waitbar(i/SNP_probeset_length,h,[num2str(i) '/' num2str(SNP_probeset_length)]);
                            counter = 0;
                        end;
                        
                        % determines if probe pair is useful; both probes have data.
                        if (isnan(probeset1(i).probe_Ratio(dataset)) == 0) && (isnan(probeset1(i+1).probe_Ratio(dataset)) == 0)
                            if (probeset1(i).probe_chromosome == analyzeChr)
                                val1 = probeset1(i).probe_Ratio(dataset);
                                val2 = probeset1(i+1).probe_Ratio(dataset);
                                val3 = (probeset_2(1).probe_length + probeset_2(i+1).probe_length)/2;
                                %val3 = (probeset_2(1).probe_Tm + probeset_2(i+1).probe_Tm)/2;
                                
                                data1(length(data1)+1) = val1;
                                data2(length(data2)+1) = val2;
                                data3(length(data3)+1) = val3;
                            end;
                        end;
                    end;
                end;
                if exist('h','var'); delete(h); clear h; end;
                
                figure;
                scatter(data1,data2,(data3-min(data3))*10+10,'b');   %show probe length as radius.
                %scatter(data1,data2,(data3-min(data3))*80+10,'b');   %show probe Tm as radius.
                ylabel('SNP_\alpha  Ratio', 'Rotation', 90, 'HorizontalAlign', 'center', 'VerticalAlign', 'middle');
                xlabel('SNP_\beta  Ratio', 'Rotation', 0, 'HorizontalAlign', 'center', 'VerticalAlign', 'middle');
                title(['Chr2 SNP probe Ratios.   (array ' names{dataset} ')']);
                text(2,3.85,'radius = probe T_m.');
                axis([0,4,0,4]);
                axis square;
                
                % het line.
                line([0 4 ],[0 4 ],'color',[0.8 0.8 0.8],'linestyle','-');
                % hom lines.
                line([0 4     ],[0 yy_hom],'color',[0.9 0.5 0.5],'linestyle','-');
                line([0 yy_hom],[0 4     ],'color',[0.9 0.5 0.5],'linestyle','-');
                % hom lines.
                line([0 4            ],[0 yy_dis_cutoff],'color',[0.9 0.5 0.5],'linestyle','-');
                line([0 yy_dis_cutoff],[0 4            ],'color',[0.9 0.5 0.5],'linestyle','-');
            end;
        end;
        
        % clean up variables not needed further.
        clear plot* bases_per_bin num_chr base bottom height i vertical_margin width left maxY val counter x1 x2 data*;
    end;
    
    %% ========================================================================
    % Allelic fraction display as done in Forche-2005
    %
    % Forche A, May G, Magee P T.   (2005)   Demonstration of Loss of
    % Heterozygosity by Single-Nucleotide Polymorphism Microarray Analysis and
    % Alterations in Strain Morphology in Candida albicans Strains during
    % Infection.   Eukaryotic Cell 4(1): 156-165.
    %
    % Forche-2005 plots allelic fraction against the square root of total
    % intensity.
    %
    % Other publications using Allelic fraction tend to plot it against the log
    % of the total intensity.
    %==========================================================================
    if (AllelicFraction_display == 1)
        for analyzeChr = chromosomes_to_analyze;
            for dataset = 1:length(names);
                % initialize vectors for scattergram analysis
                data1_0 = [];
                data2_0 = [];
                data1_1 = [];
                data2_1 = [];
                data1_2 = [];
                data2_2 = [];
                
                h = waitbar(0.0,'Please wait...','CreateCancelBtn','delete(h); clear h');
                set(h,'Name','Plotting scattergrams of SNP ratios.');
                counter = 0;
                for i = 1:SNP_probeset_length
                    if (mod(i,2) == 1)
                        counter = counter+1;
                        if (counter == 500)
                            waitbar(i/SNP_probeset_length,h,[num2str(i) '/' num2str(SNP_probeset_length)]);
                            counter = 0;
                        end;
                        
                        % determines if probe pair is useful; both probes have data.
                        if (isnan(probeset1(i).probe_Ratio(dataset)) == 0) && (isnan(probeset1(i+1).probe_Ratio(dataset)) == 0)
                            
                            % swaps probes if polarity switch indicates to do so.
                            if (probeset1(i).probe_chromosome == analyzeChr)
                                val1 = probeset1(i).probe_Ratio(dataset)/(probeset1(i).probe_Ratio(dataset) + probeset1(i+1).probe_Ratio(dataset));
                                val2 = probeset1(i+1).probe_Ratio(dataset)/(probeset1(i).probe_Ratio(dataset) + probeset1(i+1).probe_Ratio(dataset));
                                val3 = sqrt(probeset1(i).probe_ch1(dataset) + probeset1(i+1).probe_ch2(dataset));
                                if (probeset1(i).probe_polarity == 1)
                                    data1_1(length(data1_1)+1) = val1;
                                    data2_1(length(data2_1)+1) = val3;
                                elseif (probeset1(i).probe_polarity == 2)
                                    data1_2(length(data1_2)+1) = val2;
                                    data2_2(length(data2_2)+1) = val3;
                                elseif (probeset1(i).probe_polarity == 0)
                                    data1_0(length(data1_0)+1) = val1;
                                    data2_0(length(data2_0)+1) = val3;
                                else
                                    % null action for when (probe_polarity == 4)
                                    % due to probe design error; probes are identical.
                                end;
                            end;
                        end;
                    end;
                end;
                if exist('h','var'); delete(h); clear h; end;
                
                figure;
                scatter(data1_0,data2_0,5,'b');
                hold on;
                scatter(data1_1,data2_1,5,'b');
                scatter(data1_2,data2_2,5,'b');
                
                ylabel('(total fluorescence)^1^/^2', 'Rotation', 90, 'HorizontalAlign', 'center', 'VerticalAlign', 'middle');
                xlabel('Allelic Fraction', 'Rotation', 0, 'HorizontalAlign', 'center', 'VerticalAlign', 'middle');
                if (analyzeChr == 8)
                    title(['ChrR SNP probe Ratios.   (array ' names{dataset} ')']);
                else
                    title(['Chr' num2str(analyzeChr) ' SNP probe Ratios.   (array ' names{dataset} ')']);
                end;
                % y-axis in paper was from 0-6, though in this dataset it goes from 0-hundreds
                yMax = 200;
                axis([0,1,0,yMax]);
                %text(3.75,3.75,'(Het)','rotation',-45   ,'horizontalalignment','center','verticalalignment','middle','color',[0.75 0.75 0.75]);
                axis square;
                
                line([0.4 0.4],[0 yMax],'color',[0.75 0.75 0.75],'linestyle','-');
                line([0.6 0.6],[0 yMax],'color',[0.75 0.75 0.75],'linestyle','-');
            end;
        end;
        
        % clean up variables not needed further.
        clear plot* bases_per_bin num_chr base bottom height i vertical_margin width left maxY val counter x1 x2;
    end;
    
    %% ========================================================================
    % Allelic fraction vs log(Signal ratio); from O'Meara-2002, and
    % Lovmar-2003,Fan-2000.
    %==========================================================================
    if (LogFraction_display == 1)
        for analyzeChr = chromosomes_to_analyze;
            for dataset = 1:length(names);
                % initialize vectors for scattergram analysis
                data1_0 = [];
                data2_0 = [];
                data1_1 = [];
                data2_1 = [];
                data1_2 = [];
                data2_2 = [];
                
                h = waitbar(0.0,'Please wait...','CreateCancelBtn','delete(h); clear h');
                set(h,'Name','Plotting scattergrams of SNP ratios.');
                counter = 0;
                for i = 1:SNP_probeset_length
                    if (mod(i,2) == 1)
                        counter = counter+1;
                        if (counter == 500)
                            waitbar(i/SNP_probeset_length,h,[num2str(i) '/' num2str(SNP_probeset_length)]);
                            counter = 0;
                        end;
                        
                        % determines if probe pair is useful; both probes have data.
                        if (isnan(probeset1(i).probe_Ratio(dataset)) == 0) && (isnan(probeset1(i+1).probe_Ratio(dataset)) == 0)
                            
                            % swaps probes if polarity switch indicates to do so.
                            if (probeset1(i).probe_chromosome == analyzeChr)
                                val1 = log(probeset1(i).probe_Ratio(dataset)/probeset1(i+1).probe_Ratio(dataset));
                                val2 = probeset1(i).probe_Ratio(dataset)/(probeset1(i).probe_Ratio(dataset) + probeset1(i+1).probe_Ratio(dataset));
                                
                                data1_0(length(data1_0)+1) = val1;
                                data2_0(length(data2_0)+1) = val2;
                            end;
                        end;
                    end;
                end;
                if exist('h','var'); delete(h); clear h; end;
                
                figure;
                scatter(data1_0,data2_0,5,'b');
                hold on;
                scatter(data1_1,data2_1,5,'b');
                scatter(data1_2,data2_2,5,'b');
                
                ylabel('log(A/B)', 'Rotation', 90, 'HorizontalAlign', 'center', 'VerticalAlign', 'middle');
                xlabel('log(A+B)', 'Rotation', 0, 'HorizontalAlign', 'center', 'VerticalAlign', 'middle');
                if (analyzeChr == 8)
                    title(['ChrR SNP probe Ratios.   (array ' names{dataset} ')']);
                else
                    title(['Chr' num2str(analyzeChr) ' SNP probe Ratios.   (array ' names{dataset} ')']);
                end;
                % y-axis in paper was from 0-6, though in this dataset it goes from 0-hundreds
                axis([-4,4,5,12]);
                %text(3.75,3.75,'(Het)','rotation',-45   ,'horizontalalignment','center','verticalalignment','middle','color',[0.75 0.75 0.75]);
                axis square;
            end;
        end;
        
        % clean up variables not needed further.
        clear plot* bases_per_bin num_chr base bottom height i vertical_margin width left maxY val counter x1 x2;
    end;
    
    %% ========================================================================
    % Signal intensity vs Signal ratio; from Pastinen-2000.
    %==========================================================================
    if (Intensity_v_ratio_display == 1)
        for analyzeChr = chromosomes_to_analyze;
            for dataset = 1:length(names);
                % initialize vectors for scattergram analysis
                data1_0 = [];
                data2_0 = [];
                data1_1 = [];
                data2_1 = [];
                data1_2 = [];
                data2_2 = [];
                
                h = waitbar(0.0,'Please wait...','CreateCancelBtn','delete(h); clear h');
                set(h,'Name','Plotting scattergrams of SNP ratios.');
                counter = 0;
                for i = 1:SNP_probeset_length
                    if (mod(i,2) == 1)
                        counter = counter+1;
                        if (counter == 500)
                            waitbar(i/SNP_probeset_length,h,[num2str(i) '/' num2str(SNP_probeset_length)]);
                            counter = 0;
                        end;
                        
                        % determines if probe pair is useful; both probes have data.
                        if (isnan(probeset1(i).probe_Ratio(dataset)) == 0) && (isnan(probeset1(i+1).probe_Ratio(dataset)) == 0)
                            
                            % swaps probes if polarity switch indicates to do so.
                            if (probeset1(i).probe_chromosome == analyzeChr)
                                val1 = log(probeset1(i).probe_Ratio(dataset)/probeset1(i+1).probe_Ratio(dataset));
                                val2 = probeset1(i).probe_ch1(dataset) + probeset1(i+1).probe_ch2(dataset);
                                
                                data1_0(length(data1_0)+1) = val2;
                                data2_0(length(data2_0)+1) = val1;
                            end;
                        end;
                    end;
                end;
                if exist('h','var'); delete(h); clear h; end;
                
                figure;
                scatter(data1_0,data2_0,5,'b');
                hold on;
                scatter(data1_1,data2_1,5,'b');
                scatter(data1_2,data2_2,5,'b');
                
                ylabel('Log(signal ratio)', 'Rotation', 90, 'HorizontalAlign', 'center', 'VerticalAlign', 'middle');
                xlabel('Signal Intensity', 'Rotation', 0, 'HorizontalAlign', 'center', 'VerticalAlign', 'middle');
                if (analyzeChr == 8)
                    title(['ChrR SNP probe Ratios.   (array ' names{dataset} ')']);
                else
                    title(['Chr' num2str(analyzeChr) ' SNP probe Ratios.   (array ' names{dataset} ')']);
                end;
                % y-axis in paper was from 0-6, though in this dataset it goes from 0-hundreds
                %axis([1,1000,0.01,100]);
                %text(3.75,3.75,'(Het)','rotation',-45   ,'horizontalalignment','center','verticalalignment','middle','color',[0.75 0.75 0.75]);
                axis square;
            end;
        end;
        
        % clean up variables not needed further.
        clear plot* bases_per_bin num_chr base bottom height i vertical_margin width left maxY val counter x1 x2;
    end;
    
    %% ========================================================================
    % end stuff
    %==========================================================================
    fprintf('\n');
end;

%% save workspace variables for use with Run_analysis_1.m
save('Propeset_analysis_21-workspace.mat');

% delete any hanging waitbar dialogs or already saved figures.
set(0,'ShowHiddenHandles','on');
delete(get(0,'Children'));