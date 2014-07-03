function [result_image_location, archive_data_location] = ...
    process_main(microarray_design, data_file, header_rows, probeName_col, data_col_ch1, ...
                 data_col_ch2, data_col_ratio, data_col_log2ratio, phasing_dataset, ploidy_estimate, image_format, ...
                 experiment_name, workingDir, show_MRS_string)
% PROCESS_MAIN wrapper around microarray analysis pipeline for integration into online tools and website.
% This version does not load calibration data and generate a hapmap...   see 'calibration_setup.m' for relevant script.
%

%% Let the processing pipeline know that the analysis has started.
new_fid = fopen([workingDir '/working.txt'],'w');
fprintf(new_fid,'working');
fclose(new_fid);

num_chrs = 8;

fprintf('\n');
fprintf('[process_main.m] : Input variables:\n');
fprintf(['    1  : array design         : ' microarray_design  '\n']);
fprintf(['    2  : data file            : ' data_file          '\n']);
fprintf(['    3  : header rows          : ' header_rows        '\n']);
fprintf(['    4  : col probe name       : ' probeName_col      '\n']);
fprintf(['    5  : col ch1              : ' data_col_ch1       '\n']);
fprintf(['    6  : col ch2              : ' data_col_ch2       '\n']);
fprintf(['    7  : col ratio            : ' data_col_ratio     '\n']);
fprintf(['    8  : col log2ratio        : ' data_col_log2ratio '\n']);
fprintf(['    9  : phasing dataset      : ' phasing_dataset    '\n']);
fprintf(['    10 : ploidy estimate      : ' ploidy_estimate    '\n']);
fprintf(['    11 : image format         : ' image_format       '\n']);
fprintf(['    12 : experiment name      : ' experiment_name    '\n']);
fprintf(['    12 : working directory    : ' workingDir         '\n']);
fprintf(['    12 : show MRS annotations : ' show_MRS_string    '\n']);

% calibration_file = ['designs/' microarray_design '/' phasing_dataset '.mat'];
% designs/
%     links_dir/raw_SnpCgh_dir/data1/030810_DA_Agilent_16bit_A_output_fused.xls
% /
%     2.0
% .mat

[raw_data_dir, raw_data_file, raw_data_ext] = fileparts(data_file);

fprintf('\n');
fprintf('[process_main.m] : Data file parts:\n');
fprintf(['    dir  : ''' raw_data_dir '''\n']);
fprintf(['    file : ''' raw_data_file '''\n']);
fprintf(['    ext  : ''' raw_data_ext '''\n']);
raw_data_file = [raw_data_file raw_data_ext];

%% initial test case.
%clear;
%microarray_design  = 'design1';
%raw_data_dir       = 'input';
%raw_data_file      = 'example_data.txt';
%header_rows        = 46; % lines of header to be skipped before data lines.
%probeName_col      = 1;  % column containing probe ID/name.
%data_col_ch1       = 4;  % column containing probe channel 1 data.
%data_col_ch2       = 5;  % column containing probe channel 1 data.
%data_col_ratio     = 6;  % column containing probe channel ratio data.
%data_col_log2ratio = 7;  % column containing probe channel log2ratio data.
%phasing_dataset    = 'design1.cal_paper';
%ploidy_estimate    = 2.0;
%image_format       = 'png';
%experiment_name    = 'test case 1';

figureDir           = workingDir;
matlab_save_dir     = workingDir;

ploidy_estimate     = str2num(ploidy_estimate);

tic;
% =========================================================================
% Generate graphs involved in analysis of SNP probe data from 1' Agilent
% slide.
%==========================================================================
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
%    bases_per_bin              : Controls bin sizes for SNP/CGH fractions of plot.
%    scale_type                 : 'Ratio' or 'Log2Ratio' y-axis scaling of copy number.
%                                 'Log2Ratio' does not properly scale CGH data by ploidy.
%    Chr_max_width              : max width of chromosomes as fraction of figure width.
SNP_Genomic_display              = true;
    CGH_Genomic_display          = true;
    Linear_display               = true;
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
	if strcmp(show_MRS_string,'1')
	    show_MRS                 = true;
	else
		show_MRS				 = false;
	end;
    Centromere_format            = 2;
    scale_type                   = 'Ratio';
    Chr_max_width                = 0.8;

    %not in place yet... uninformative currently looks like unnassigned to script.
    show_uninformative           = false;

% Generate CGD annotation files for analyzed microarrays.
Output_CGD_annotations           = true;

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
chromosomes_to_analyze = 1:num_chrs;

% Defines chromosome sizes in bp. (diploid total=28,567,7888)
% Defines centromere locations in bp.
% Defines MRS locations in bp.
[centromeres, chr_sizes,MRSs] = Load_information();
for i = chromosomes_to_analyze
    chr_size(i)  = chr_sizes(i).size;
    cen_start(i) = centromeres(i).start;
    cen_end(i)   = centromeres(i).end;
end;
for i = 1:length(MRSs)
    MRS_chr(i)   = MRSs(i).chr;
    MRS_location(i) = (MRSs(i).start+MRSs(i).end)/2;
end;
clear centromeres chr_sizes MRSs;
bases_per_bin = max(chr_size)/700;

%%=========================================================================
%%=========================================================================
%%= No further control variables below. ===================================
%%=========================================================================
%%=========================================================================


switch phasing_dataset
    case 'none'
        show_uncalibrated    = true;
        calibration_set      = 0;
        no_calibration       = 1;
        files_found          = 1;
    otherwise
        show_uncalibrated    = false;
        calibration_set      = 1;
        no_calibration       = 0;
        files_found          = 1;
        calibration_file = ['designs/' microarray_design '/' phasing_dataset '.mat'];
        load(calibration_file);
end;

%% ====================================================================
% Setup for experimental data analysis using calibration data.
% Swap experimental SNP data per hapmap/calibration file.
%----------------------------------------------------------------------
%info from design considerations.
load(['designs/' microarray_design '/SNP_probeset_2.mat']);
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

% process raw data file into structures containing CGH and SNP information (if not already done).

data_file_load_online(matlab_save_dir,raw_data_dir,raw_data_file, microarray_design, experiment_name);

% load experimental values for all the probes.
load([matlab_save_dir '/' experiment_name '.' microarray_design '.SNP_data.mat']);
load([matlab_save_dir '/' experiment_name '.' microarray_design '.CGH_data.mat']);

SNP_probeset_length = length(probeset1);
CGH_probeset_length = length(probeset2);


%% ====================================================================
% Get segmental aneuploidy data for experiment from common_ChARM file.
%----------------------------------------------------------------------
% Initializes vectors used to hold CGH data.
for chromosome = 1:num_chrs   % eight chromosomes.
    for j = 1:2   % two CGH data catagories tracked.
        chr_CGHdata{chromosome,j} = zeros(1,ceil(chr_size(chromosome)/bases_per_bin));
    end;
end;
% tally CGH data into bins for main figure histogram.
counter = 0;

for i = 1:CGH_probeset_length
	counter = counter+1;
	val  = ceil((probeset2(i).probe_location)/bases_per_bin);
	if (strcmp(scale_type,'Ratio') == 1)
		val2 = probeset2(i).probe_Ratio;
	else
		val2 = probeset2(i).probe_Log2Ratio;
	end;
	% Determines distribution of CGH data.
	if (length(val2) > 0)
		% count of data points at CGH locus.
		chr_CGHdata{probeset2(i).probe_chromosome,1}(val) = chr_CGHdata{probeset2(i).probe_chromosome,1}(val)+1;
		% total data at locus.
		chr_CGHdata{probeset2(i).probe_chromosome,2}(val) = chr_CGHdata{probeset2(i).probe_chromosome,2}(val)+val2;
	end;
end;
% divide total Ratio values by number of points per bin.
for chr = 1:num_chrs
    for j = 1:length(chr_CGHdata{chr,2})
        if (chr_CGHdata{chr,1}(j) == 0)
            chr_CGHdata{chr,2}(j) = 1;
        else
            chr_CGHdata{chr,2}(j) = chr_CGHdata{chr,2}(j)/chr_CGHdata{chr,1}(j);
        end;
    end;
    CNVplot2{chr} = chr_CGHdata{chr,2};
end;
% Save presented CNV data in a file format common across data types being processed.
fprintf('\nSaving "common_CNV" data file.');
genome_CNV = 'Ca_a';
save([workingDir '/' experiment_name '.common_CNV.mat'], 'CNVplot2','genome_CNV');

% ChARM_v3(experiment_name, workingDir,figureDir);
% [segmental_aneuploidy] = Load_dataset_information_1(experiment_name,workingDir);
segmental_aneuploidy = [];

%% ====================================================================
% Apply probe polarity assignment from calibration data to experimental data.
%----------------------------------------------------------------------
fprintf(['\nPhasing experimental SNP data from: ' strrep(raw_data_dir,'\','\\')]);
for i = 1:SNP_probeset_length
    if (no_calibration == 1)
        probeset1(i).probe_polarity = 0;
    else
        probeset1(i).probe_polarity = SNP_probeset(i).probe_polarity;
    end;
end;

%% ====================================================================
% Determine chromosome copy numbers.
%----------------------------------------------------------------------
datasetDetails = [];
fprintf('\nDetermining chromsome copy numbers for microarray.');
[chr_breaks, chrCopyNum] = FindChrSizes(segmental_aneuploidy,CGH_probeset_length,probeset2,chr_size,ploidy_estimate);
datasetDetails.chr_breaks = chr_breaks;
datasetDetails.chrCopyNum = chrCopyNum;

    
%% ====================================================================
% Determine cutoffs for experimental dataset chromosome segments
%----------------------------------------------------------------------

fprintf('\nDetermining SNP interpretation cutoffs for microarray.');
% Finds initial homozygous peak locations.
fprintf(['\n chr1 breaks = ' num2str(chr_breaks{1})]);
fprintf(['\n chr2 breaks = ' num2str(chr_breaks{2})]);
fprintf(['\n chr3 breaks = ' num2str(chr_breaks{3})]);
fprintf(['\n chr4 breaks = ' num2str(chr_breaks{4})]);
fprintf(['\n chr5 breaks = ' num2str(chr_breaks{5})]);
fprintf(['\n chr6 breaks = ' num2str(chr_breaks{6})]);
fprintf(['\n chr7 breaks = ' num2str(chr_breaks{7})]);
fprintf(['\n chr8 breaks = ' num2str(chr_breaks{8})]);
fprintf(['\n']);
[realHomozygous_peak, disomy_fit,skew_factor] = ...
    FindRealHomozygousPeaks_2(chrCopyNum,SNP_probeset_length,probeset1,chr_breaks,chr_size,show_unnassigned,DataTypeToUse,show_fitting);


% Finds real peak locations.
[monosomy_peak,disomy_peak,trisomy_peak,tetrasomy_peak,pentasomy_peak,hexasomy_peak ] = ...
    FindPeaks(realHomozygous_peak);
% Determine cutoffs between peaks for each datasets:chromosome:segment.
for chromosome = chromosomes_to_analyze;
    for segment = 1:length(chrCopyNum{chromosome})
        chrCopyNum = datasetDetails.chrCopyNum;
        chr_breaks = datasetDetails.chr_breaks;
        % Determins cutoffs for a single chromosome segment, using intersections of Gaussian fit curves.
        MakeFigure = Gaussian_fit_display;
        [raw,smoothed,peaks,actual_cutoffs,mostLikelyGaussians,chrCopyNum] = ...
            FindGaussianCutoffs_2(probeset1,chrCopyNum,chr_breaks,chr_size,chromosome,...
            segment,monosomy_peak,disomy_peak,trisomy_peak,tetrasomy_peak,pentasomy_peak,...
            hexasomy_peak,skew_factor,experiment_name,raw_data_dir,Gaussian_fit_display,show_fitting,...
            DataTypeToUse);
        datasetDetails.histogram_raw{chromosome,segment}      = raw;
        datasetDetails.histogram_smooth{chromosome,segment}   = smoothed;
        datasetDetails.peaks{chromosome,segment}              = peaks/200;
        datasetDetails.cutoffs{chromosome,segment}            = actual_cutoffs/200;
        datasetDetails.importantGaussians{chromosome,segment} = mostLikelyGaussians;
        datasetDetails.chrCopyNum                             = chrCopyNum;
    end;
end;
    
%% ====================================================================
% Save datasetDetails.
%----------------------------------------------------------------------
save([matlab_save_dir '/' experiment_name '.' microarray_design '.datasetDetails.mat'], 'datasetDetails');
    
%% ========================================================================
% Plotting probes across genome by interpretation catagory after accounting
% for polarity of SNP pairs.
%    CGH info plotted as histogram.
%    SNP info plotted as colorbar.
%==========================================================================
fprintf(['\nGenerating SNP/CGH figures from: ' strrep(raw_data_dir,'\','\\')]);
% define number of and labels for chromosomes.
chr_labels = {'Chr1','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7','ChrR'};
num_chr = length(chr_labels);

bases_per_bin = max(chr_size)/700;
chr_length_scale_multiplier = 1/bases_per_bin;
maxY = 10;

% Initialize linear view figure plot parameters.
if (Linear_display == true)
    Linear_fig = figure(2);
    Linear_genome_size     = sum(chr_size);
    Linear_Chr_max_width   = 0.85;
    Linear_left_start      = 0.07;
    Linear_left_chr_gap    = 0.01;
    Linear_height          = 0.6;
    Linear_base            = 0.1;
    Linear_TickSize        = -0.01;  %negative for outside, percentage of longest chr figure.
    Linear_maxY            = 10;
end;

% Initializes vectors used to hold number of SNPs in each
%    interpretation catagory.
for chromosome = 1:num_chrs   % eight chromosomes.
    for j = 1:14   % 14 SNP interpretation catagories tracked.
        chr_SNPdata{chromosome,j} = zeros(1,ceil(chr_size(chromosome)/bases_per_bin));
    end;
    for j = 1:2   % two CGH data catagories tracked.
        chr_CGHdata{chromosome,j} = zeros(1,ceil(chr_size(chromosome)/bases_per_bin));
    end;
end;

% tally CGH data into bins for main figure histogram.
% h = waitbar(0.0,'Please wait...','CreateCancelBtn','delete(h); clear h');
% set(h,'Name','Plotting CGH data.');
counter = 0;
for i = 1:CGH_probeset_length
    counter = counter+1;
    %if (counter == 500)
    %    waitbar(i/CGH_probeset_length,h,[num2str(i) '/' num2str(CGH_probeset_length)]);
    %    counter = 0;
    %end;
    % val is the genomic location of the probePair.
    val  = ceil((probeset2(i).probe_location)/bases_per_bin);
    if (strcmp(scale_type,'Ratio') == 1)
        val2 = probeset2(i).probe_Ratio;
    else
        val2 = probeset2(i).probe_Log2Ratio;
    end;
    % Determines distribution of CGH data.
    if (length(val2) > 0)
        % count of data points at CGH locus.
        chr_CGHdata{probeset2(i).probe_chromosome,1}(val) = chr_CGHdata{probeset2(i).probe_chromosome,1}(val)+1;
        % total data at locus.
        chr_CGHdata{probeset2(i).probe_chromosome,2}(val) = chr_CGHdata{probeset2(i).probe_chromosome,2}(val)+val2;
    end;
end;
% divide total Ratio values by number of points per bin.
for chr = 1:num_chrs
    for j = 1:length(chr_CGHdata{chr,2})
        if (chr_CGHdata{chr,1}(j) == 0)
            chr_CGHdata{chr,2}(j) = 1;
        else
            chr_CGHdata{chr,2}(j) = chr_CGHdata{chr,2}(j)/chr_CGHdata{chr,1}(j);
        end;
    end;
    CNVplot2{chr} = chr_CGHdata{chr,2};
end;

% Save presented CNV data in a file format common across data types being processed.
fprintf('\nSaving "common_CNV" data file.');
genome_CNV = 'Ca_a';
save([workingDir '/' experiment_name '.common_CNV.mat'], 'CNVplot2','genome_CNV');

% Define colors for use in figures;
colorNoData = [1.0    1.0    1.0   ]; %used when no data is available for the bin.
colorInit   = [0.5    0.5    0.5   ]; %external; used in blending at ends of chr.
[colorA,colorB, colorAB, colorAAB,colorABB, colorAAAB,colorABBB, colorAAAAB,colorAAABB,...
    colorAABBB,colorABBBB, colorAAAAAB,colorABBBBB, colorPeak,colorCutoff ] = DefineColors();
if (show_unnassigned == true) || (show_uncalibrated == true)
    colorUn_Hom  = [0.5   0.5   1.0  ]; % homozygous unassigned.
    colorUn_DHet = [0.667 0.667 0.667]; % disomy homozygous unassigned.
    colorUn_THet = [0.583 0.583 0.833]; % trisomy homozygous unassigned.
else
    colorUn_Hom  = [1.0   1.0   1.0  ]; % homozygous unassigned.
end;


% Initialize CGD annotation output file.
if (Output_CGD_annotations == true)
    CGDid = fopen([workingDir '/CGD_annotations.' experiment_name '.txt'], 'w');
    fprintf(CGDid,['track name=' experiment_name ' description="SNP/CGH microarray annotation of SNPs" useScore=0 itemRGB=On\n']);
end;

% h = waitbar(0.0,'Please wait...','CreateCancelBtn','delete(h); clear h');
% set(h,'Name','Plotting SNP data.');
counter = 0;
SNPs_hom = 0;
SNPs_total = 0;
for i = 1:2:SNP_probeset_length
    % probeset1 consists of sequential pairs, representing each allele for each SNP.
    %    The odd numbered probes are the first of each pair.
    counter = counter+1;
    %if (counter == 500)
    %    waitbar(i/SNP_probeset_length,h,[num2str(i) '/' num2str(SNP_probeset_length)]);
    %    counter = 0;
    %end;
    
    % determines if probe pair is useful; both probes have data.
    if (length(probeset1(i).probe_Ratio) > 0) && (length(probeset1(i+1).probe_Ratio) > 0)
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
        for iii = 1
            % Find chromosome of probe.
            chromosome = probeset1(i).probe_chromosome;
            
            % Find which segment the probe is found on.
            if (length(datasetDetails.chr_breaks{chromosome}) > 2)
                location = probeset1(i).probe_location/chr_size(chromosome);
                for last_break = length(datasetDetails.chr_breaks{chromosome}):-1:2
                    if (location < datasetDetails.chr_breaks{chromosome}(last_break))
                        segment = last_break-1;
                    end;
                end;
            else
                segment = 1;
            end;
            
            % Calculate value of SNP probe pair.
            [UsedData] = calculateValue(probeset1,i,DataTypeToUse);
            
            % Determine which section of the allelic fraction histogram the probe pair falls into.
            done = false;
            section = datasetDetails.importantGaussians{chromosome,segment}(length(datasetDetails.importantGaussians{chromosome,segment}));
            for cutoff = 1:length(datasetDetails.cutoffs{chromosome,segment})
                if (done == false) && (UsedData <= datasetDetails.cutoffs{chromosome,segment}(cutoff))
                    done = true;
                    section = datasetDetails.importantGaussians{chromosome,segment}(cutoff);
                end;
            end;
            
            % Determines probe assignment of SNP pair.
            ChrCopyNumber = datasetDetails.chrCopyNum{probeset1(i).probe_chromosome}(segment);
            if (probeset1(i).probe_polarity == 0)
		% hom unnassigned.
                probeset1(i).probe_assignment   = 13;
                probeset1(i+1).probe_assignment = 13;
                if (show_unnassigned == true) && (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
		% start BED file output for this SNP.
		outputCGDannotationLine(CGDid,probeset1,i,Output_CGD_annotations, colorUn_Hom);
            elseif (probeset1(i).probe_polarity == 4)
                % null action for when (probe_polarity == 4) due to probe design error; probes are identical.
            elseif (ChrCopyNumber <= 1)     %monosomy
                % Assigns probe pairs to specific interpretations, collects data for each interpretation.
                if (section == 1)   % monosomy:'a'
                    probeset1(i).probe_assignment   = 1;
                    probeset1(i+1).probe_assignment = 1;
                    if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
		    outputCGDannotationLine(CGDid,probeset1,i,Output_CGD_annotations, colorA);
                else % (section == 2)   % monosomy:'b'
                    probeset1(i).probe_assignment   = 2;
                    probeset1(i+1).probe_assignment = 2;
                    if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
		    outputCGDannotationLine(CGDid,probeset1,i,Output_CGD_annotations, colorB);
                end;
            elseif (ChrCopyNumber <= 2) %disomy
                % Assigns probe pairs to specific interpretations, collects data for each interpretation.
                if (section == 1)   % disomy:'aa'
                    probeset1(i).probe_assignment   = 1;
                    probeset1(i+1).probe_assignment = 1;
                    if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
		    outputCGDannotationLine(CGDid,probeset1,i,Output_CGD_annotations, colorA);
                elseif (section == 3)   % disomy:'bb'
                    probeset1(i).probe_assignment   = 2;
                    probeset1(i+1).probe_assignment = 2;
                    if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
		    outputCGDannotationLine(CGDid,probeset1,i,Output_CGD_annotations, colorB);
                else % (section == 2)   % disomy:'ab'
                    probeset1(i).probe_assignment   = 0;
                    probeset1(i+1).probe_assignment = 0;
                    if (Show_Genomic_LOH_fraction == true);   SNPs_total = SNPs_total+1;   end;
		    outputCGDannotationLine(CGDid,probeset1,i,Output_CGD_annotations, colorAB);
                end;
            elseif (ChrCopyNumber <= 3) %trisomy
                % Assigns probe pairs to specific interpretations, collects data for each interpretation.
                if (section == 1)   % trisomy:'aaa'
                    probeset1(i).probe_assignment   = 1;
                    probeset1(i+1).probe_assignment = 1;
                    if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
		    outputCGDannotationLine(CGDid,probeset1,i,Output_CGD_annotations, colorA);
                elseif (section == 4)   % trisomy:3 'bbb'
                    probeset1(i).probe_assignment   = 2;
                    probeset1(i+1).probe_assignment = 2;
                    if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
		    outputCGDannotationLine(CGDid,probeset1,i,Output_CGD_annotations, colorB);
                elseif (section == 2)   % trisomy:'aab'
                    probeset1(i).probe_assignment   = 3;
                    probeset1(i+1).probe_assignment = 3;
                    if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
		    outputCGDannotationLine(CGDid,probeset1,i,Output_CGD_annotations, colorAAB);
                else % (section == 3)   % trisomy:'abb'
                    probeset1(i).probe_assignment   = 4;
                    probeset1(i+1).probe_assignment = 4;
                    if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
		    outputCGDannotationLine(CGDid,probeset1,i,Output_CGD_annotations, colorABB);
                end;
            elseif (ChrCopyNumber <= 4) %tetrasomy
                % Assigns probe pairs to specific interpretations, collects data for each interpretation.
                if (section == 1)   % tetrasomy:'aaaa'
                    probeset1(i).probe_assignment   = 1;
                    probeset1(i+1).probe_assignment = 1;
                    if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
		    outputCGDannotationLine(CGDid,probeset1,i,Output_CGD_annotations, colorA);
                elseif (section == 5)   % tetrasomy:'bbbb'
                    probeset1(i).probe_assignment   = 2;
                    probeset1(i+1).probe_assignment = 2;
                    if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
		    outputCGDannotationLine(CGDid,probeset1,i,Output_CGD_annotations, colorB);
                elseif (section == 2)   % tetrasomy:'aaab'
                    probeset1(i).probe_assignment   = 5;
                    probeset1(i+1).probe_assignment = 5;
                    if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
		    outputCGDannotationLine(CGDid,probeset1,i,Output_CGD_annotations, colorAAAB);
                elseif (section == 4)   % tetrasomy:'abbb'
                    probeset1(i).probe_assignment   = 6;
                    probeset1(i+1).probe_assignment = 6;
                    if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
		    outputCGDannotationLine(CGDid,probeset1,i,Output_CGD_annotations, colorABBB);
                else % (section == 3)   % tetrasomy:'aabb'
                    probeset1(i).probe_assignment   = 0;
                    probeset1(i+1).probe_assignment = 0;
                    if (Show_Genomic_LOH_fraction == true);   SNPs_total = SNPs_total+1;   end;
		    outputCGDannotationLine(CGDid,probeset1,i,Output_CGD_annotations, colorAB);
                end;
            elseif (ChrCopyNumber <= 5) %pentasomy
                % Assigns probe pairs to specific interpretations, collects data for each interpretation.
                if (section == 1)   % pentasomy:'aaaaa'
                    probeset1(i).probe_assignment   = 1;
                    probeset1(i+1).probe_assignment = 1;
                    if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
		    outputCGDannotationLine(CGDid,probeset1,i,Output_CGD_annotations, colorA);
                elseif (section == 6)   % pentasomy:'bbbbb'
                    probeset1(i).probe_assignment   = 2;
                    probeset1(i+1).probe_assignment = 2;
                    if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
		    outputCGDannotationLine(CGDid,probeset1,i,Output_CGD_annotations, colorB);
                elseif (section == 2)   % pentasomy:'aaaab'
                    probeset1(i).probe_assignment   = 7;
                    probeset1(i+1).probe_assignment = 7;
                    if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
		    outputCGDannotationLine(CGDid,probeset1,i,Output_CGD_annotations, colorAAAAB);
                elseif (section == 5)   % pentasomy:'abbbb'
                    probeset1(i).probe_assignment   = 10;
                    probeset1(i+1).probe_assignment = 10;
                    if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
		    outputCGDannotationLine(CGDid,probeset1,i,Output_CGD_annotations, colorABBBB);
                elseif (section == 3)   % pentasomy:'aaabb'
                    probeset1(i).probe_assignment   = 8;
                    probeset1(i+1).probe_assignment = 8;
                    if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
		    outputCGDannotationLine(CGDid,probeset1,i,Output_CGD_annotations, colorAAABB);
                else % (section == 4)   % pentasomy:'aabbb'
                    probeset1(i).probe_assignment   = 9;
                    probeset1(i+1).probe_assignment = 9;
                    if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
		    outputCGDannotationLine(CGDid,probeset1,i,Output_CGD_annotations, colorAABBB);
                end;
            else %treat as hexasomy: if (ChrCopyNumber <= 6) %hexasomy
                % Assigns probe pairs to specific interpretations, collects data for each interpretation.
                if (section == 1)   % hexasomy:'aaaaaa'
                    probeset1(i).probe_assignment   = 1;
                    probeset1(i+1).probe_assignment = 1;
                    if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
		    outputCGDannotationLine(CGDid,probeset1,i,Output_CGD_annotations, colorA);
                elseif (section == 7)   % hexasomy:'bbbbbb'
                    probeset1(i).probe_assignment   = 2;
                    probeset1(i+1).probe_assignment = 2;
                    if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
		    outputCGDannotationLine(CGDid,probeset1,i,Output_CGD_annotations, colorB);
                elseif (section == 2)   % hexasomy:'aaaaab'
                    probeset1(i).probe_assignment   = 11;
                    probeset1(i+1).probe_assignment = 11;
                    if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
		    outputCGDannotationLine(CGDid,probeset1,i,Output_CGD_annotations, colorAAAAAB);
                elseif (section == 6)   % hexasomy:'abbbbb'
                    probeset1(i).probe_assignment   = 12;
                    probeset1(i+1).probe_assignment = 12;
                    if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
		    outputCGDannotationLine(CGDid,probeset1,i,Output_CGD_annotations, colorABBBBB);
                elseif (section == 3)   % hexasomy:'aaaabb'
                    probeset1(i).probe_assignment   = 3;
                    probeset1(i+1).probe_assignment = 3;
                    if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
		    outputCGDannotationLine(CGDid,probeset1,i,Output_CGD_annotations, colorAAB);
                elseif (section == 5)   % hexasomy:'aabbbb'
                    probeset1(i).probe_assignment   = 4;
                    probeset1(i+1).probe_assignment = 4;
                    if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
		    outputCGDannotationLine(CGDid,probeset1,i,Output_CGD_annotations, colorABB);
                else % (section == 4)   % hexasomy:'aaabbb'
                    probeset1(i).probe_assignment   = 0;
                    probeset1(i+1).probe_assignment = 0;
                    if (Show_Genomic_LOH_fraction == true);   SNPs_hom = SNPs_hom+1;   SNPs_total = SNPs_total+1;   end;
		    outputCGDannotationLine(CGDid,probeset1,i,Output_CGD_annotations, colorAB);
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

% basic plot parameters.
left            = 0.15;
height          = 0.5/num_chr;
base            = 0.1;
vertical_margin = 0.3/num_chr;
TickSize        = -0.005;  %negative for outside, percentage of longest chr figure.

if (Linear_display == true)
    Linear_left = Linear_left_start;
end;

% main figure calculation and generation.
fig = figure(1);
set(gcf, 'Position', [0 70 1024 600]);
for chromosome = 1:num_chrs
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
	numSegments            = length(datasetDetails.chrCopyNum{chromosome});
	histAll                = [];
	histAll2               = [];
	smoothed               = [];
	smoothed2              = [];
	histAll{numSegments}   = [];
	histAll2{numSegments}  = [];
	smoothed{numSegments}  = [];
	smoothed2{numSegments} = [];
        for segment = 1:numSegments
            subplot('Position',[(0.15 + Chr_max_width*chr_size(chromosome)/max(chr_size) + 0.005)+width*(segment-1) bottom width height]);
            for i = round(1+length(usedPlotCGH)*datasetDetails.chr_breaks{chromosome}(segment)):round(length(usedPlotCGH)*datasetDetails.chr_breaks{chromosome}(segment+1))
                if (strcmp(scale_type,'Ratio') == 1)
                    if (ploidy_estimate == 0) % no value; no scaling.
                        y_ = usedPlotCGH(i);
                    else
                        if (Yscale_nearest_even_ploidy == true)
                            if (abs(2-ploidy_estimate) < abs(4-ploidy_estimate)) % nearer to diploid
                                y_ = usedPlotCGH(i)*ploidy_estimate/2;
                            else
                                y_ = usedPlotCGH(i)*ploidy_estimate/4;
                            end;
                        else
                            y_ = usedPlotCGH(i)*ploidy_estimate/2;
                        end;
                    end;
                else % scale_type = Log2Ratio.
                    if (ploidy_estimate == 0) % no value; no scaling.
                        y_ = usedPlotCGH(i)+maxY/2;
                    else
                        y_ = log2(pow2(usedPlotCGH(i))*ploidy_estimate/2)+maxY/2;
                    end;
                end;
                histAll{segment}(i) = y_;
            end;
            % make a histogram of CGH data, then smooth it for display.
	    fprintf(['\nchr = ' num2str(chromosome) '; segment = ' num2str(segment) '; ']);
	    fprintf(['length(histAll{segment}) = ' num2str(length(histAll{segment})) '\n']);

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
            subplot('Position',[(0.15 + Chr_max_width*chr_size(chromosome)/max(chr_size) + 0.005 + width*(length(datasetDetails.chrCopyNum{chromosome})-1) + width+0.001) bottom width height]);
        else
            subplot('Position',[(0.15 + Chr_max_width*chr_size(chromosome)/max(chr_size) + 0.005) bottom width height]);
        end;
        axis off square;
        set(gca,'YTick',[]);
        set(gca,'XTick',[]);
        if (length(datasetDetails.chrCopyNum{chromosome}) == 1)
            chr_string = num2str(datasetDetails.chrCopyNum{chromosome}(1));
        else
            chr_string = num2str(datasetDetails.chrCopyNum{chromosome}(1));
            for i = 2:length(datasetDetails.chrCopyNum{chromosome})
                chr_string = [chr_string ',' num2str(datasetDetails.chrCopyNum{chromosome}(i))];
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
        
        for segment = 1:length(datasetDetails.chrCopyNum{chromosome})
            if (segment == 1)
                subplot('Position',[0.03 bottom width (height/length(datasetDetails.chrCopyNum{chromosome}))]);
            else
                subplot('Position',[0.03 (bottom+height/length(datasetDetails.chrCopyNum{chromosome})*(segment-1)) width (height/length(datasetDetails.chrCopyNum{chromosome}))]);
            end;
            histAll = [];
            histAll2 = [];
            %for i = round(1+length(usedPlotCGH)*chr_breaks{chromosome}(segment))
            %        round(  length(usedPlotCGH)*chr_breaks{chromosome}(segment+1))
            for i = 1:2:SNP_probeset_length
                if (probeset1(i).probe_chromosome == chromosome) && ...
                        (probeset1(i).probe_location >= round(1+chr_size(chromosome)*datasetDetails.chr_breaks{chromosome}(segment))) && ...
                        (probeset1(i).probe_location <= round(chr_size(chromosome)*datasetDetails.chr_breaks{chromosome}(segment+1)))
                    if (length(probeset1(i).probe_Ratio) > 0) && (length(probeset1(i+1).probe_Ratio) > 0)
                        % Calculate value of SNP probe pair.
                        [UsedData] = calculateValue(probeset1,i,DataTypeToUse);
                        
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

	    % ddd
            % Invert histogram values;
            histAll = 1-histAll;
	    % fprintf(['\nRaw Angleplot histogram data: \n\t' num2str(histAll) '\n']);

            smoothed = smooth_gaussian(hist(histAll,200),3,20);
            % make a smoothed version of just the endpoints used to ensure histogram bounds.
            histAll2(1) = 0;
            histAll2(2) = 1;
            smoothed2 = smooth_gaussian(hist(histAll2,200),3,20);
            % subtract the smoothed endpoints from the histogram to remove the influence of the added endpoints.
            smoothed = (smoothed-smoothed2);
            smoothed = smoothed/max(smoothed);
            
            hold on;
            copynum    = datasetDetails.chrCopyNum{chromosome}(segment);
            region_    = 0;
            for region = datasetDetails.importantGaussians{chromosome,segment}
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
                if (length(datasetDetails.importantGaussians{chromosome,segment}) == 1)
                    area(1:200,smoothed(1:200),'FaceColor',color,'EdgeColor',color);
                else
                    if (region_ == 1)
                        area(         round(200-datasetDetails.cutoffs{chromosome,segment}(region_  )*200): ...
			              200                                                                 , ...
			     smoothed(round(200-datasetDetails.cutoffs{chromosome,segment}(region_  )*200): ...
                                      200                                                                 ) ...
                             ,'FaceColor',color,'EdgeColor',color);
                    elseif (region_ == length(datasetDetails.importantGaussians{chromosome,segment}))
                        area(        1                                                                   :  ...
			             round(200-datasetDetails.cutoffs{chromosome,segment}(region_-1)*200),  ...
		 	     smoothed(1                                                                   : ...
                                      round(200-datasetDetails.cutoffs{chromosome,segment}(region_-1)*200)) ...
      		             ,'FaceColor',color,'EdgeColor',color);
                    else
                        area(         round(200-datasetDetails.cutoffs{chromosome,segment}(region_  )*200): ...
                                      round(200-datasetDetails.cutoffs{chromosome,segment}(region_-1)*200), ...
                             smoothed(round(200-datasetDetails.cutoffs{chromosome,segment}(region_  )*200): ...
                                      round(200-datasetDetails.cutoffs{chromosome,segment}(region_-1)*200)) ...
			     ,'FaceColor',color,'EdgeColor',color);
                    end;
                end;
            end;
            for peak = 1:length(datasetDetails.peaks{chromosome,segment})
                plot([200-datasetDetails.peaks{chromosome,segment}(peak)*200; 200-datasetDetails.peaks{chromosome,segment}(peak)*200],[0; 1],'color',colorPeak);
            end;
            for cutoff = 1:length(datasetDetails.cutoffs{chromosome,segment})
                plot([200-datasetDetails.cutoffs{chromosome,segment}(cutoff)*200; 200-datasetDetails.cutoffs{chromosome,segment}(cutoff)*200],[0; 1],'color',colorCutoff);
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
                if (ploidy_estimate == 0) % no value; no scaling.
                    y_ = [maxY/2 usedPlotCGH(i) usedPlotCGH(i) maxY/2];
                else
                    if (Yscale_nearest_even_ploidy == true)
                        if (abs(2-ploidy_estimate) < abs(4-ploidy_estimate)) % nearer to diploid
                            y_ = [maxY/2 usedPlotCGH(i)*ploidy_estimate/2 usedPlotCGH(i)*ploidy_estimate/2 maxY/2];
                        else
                            y_ = [maxY/2 usedPlotCGH(i)*ploidy_estimate/4 usedPlotCGH(i)*ploidy_estimate/4 maxY/2];
                        end;
                    else
                        y_ = [maxY/2 usedPlotCGH(i)*ploidy_estimate/2 usedPlotCGH(i)*ploidy_estimate/2 maxY/2];
                    end;
                end;
            else % scale_type = Log2Ratio.
                if (ploidy_estimate == 0) % no value; no scaling.
                    y_ = [maxY/2 usedPlotCGH(i)+maxY/2 usedPlotCGH(i)+maxY/2 maxY/2];
                else
                    %y_ = [maxY/2 usedPlotCGH(i)+maxY/2 usedPlotCGH(i)+maxY/2 maxY/2];
                    y_ = [maxY/2 log2(pow2(usedPlotCGH(i))*ploidy_estimate/2)+maxY/2 ...
                        log2(pow2(usedPlotCGH(i))*ploidy_estimate/2)+maxY/2 maxY/2];
                end;
            end;
            % makes a blackbar for each bin.
            f = fill(x_,y_,c_);
            set(f,'linestyle','none');
        end;
        x2 = chr_size(chromosome)*chr_length_scale_multiplier;
        if (strcmp(scale_type,'Ratio') == 1)
	    plot([0; x2], [maxY/2;maxY/2],'Color',[0 0 0]);  % 2n line.
	    line([0 x2], [maxY/4*3 maxY/4*3],'Color',[0.85 0.85 0.85]); % top line.
	    line([0 x2], [maxY/4   maxY/4],  'Color',[0.85 0.85 0.85]); % bottom line.
        else
            plot([0; x2], [maxY/2;maxY/2],'Color',[0 0 0]);  % 2n line.
	    line([0 x2], [maxY/2*log2(3) maxY/2*log2(3)], 'Color',[0.85 0.85 0.85]); % top cen.
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
		fill([leftEnd  leftEnd  (leftEnd+dx) ],[(maxY-dy) maxY       maxY], colorNoData, 'EdgeColor', colorNoData);
		fill([rightEnd rightEnd (rightEnd-dx)],[(maxY-dy) maxY       maxY], colorNoData, 'EdgeColor', colorNoData);
		fill([leftEnd  leftEnd  (leftEnd+dx) ],[(0+dy)    0          0   ], colorNoData, 'EdgeColor', colorNoData);
		fill([rightEnd rightEnd (rightEnd-dx)],[(0+dy)    0          0   ], colorNoData, 'EdgeColor', colorNoData);
		fill([(x1-dx)  x1       (x1+dx)      ],[maxY      (maxY-dy)  maxY], colorNoData, 'EdgeColor', colorNoData);
		fill([(x1-dx)  x1       (x1+dx)      ],[0         dy         0   ], colorNoData, 'EdgeColor', colorNoData);
		plot([leftEnd leftEnd (leftEnd+dx) x1-dx x1      x2      x2+dx rightEnd-dx rightEnd rightEnd rightEnd-dx x2+dx x2 x1 x1-dx leftEnd+dx leftEnd], ...
			 [dy      maxY-dy maxY         maxY  maxY-dy maxY-dy maxY  maxY        maxY-dy  dy       0           0     dy dy 0     0          dy     ], ...
			 'Color',[0 0 0]);
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
                if (abs(2-ploidy_estimate) < abs(4-ploidy_estimate)) % nearer to diploid
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
        title(experiment_name,'Interpreter','none','FontSize',12);
    end;

    %% Linear figure draw section
    if (Linear_display == true)
	figure(Linear_fig);
	Linear_width = Linear_Chr_max_width*chr_size(chromosome)/Linear_genome_size;
	subplot('Position',[Linear_left Linear_base Linear_width Linear_height]);
	Linear_left = Linear_left + Linear_width + Linear_left_chr_gap;
	hold on;
	title(chr_labels(chromosome),'Interpreter','none','FontSize',10);

	% display color bars for this chromosome.
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
	% end color bars section.

	% cgh plot section.
	c_ = [0 0 0];
	for i = 1:length(usedPlotCGH);
	    x_ = [i i i-1 i-1];
	    if (strcmp(scale_type,'Ratio') == 1)
		if (ploidy_estimate == 0) % no value; no scaling.
		    y_ = [maxY/2 usedPlotCGH(i) usedPlotCGH(i) maxY/2];
		else
		    if (Yscale_nearest_even_ploidy == true)
			if (abs(2-ploidy_estimate) < abs(4-ploidy_estimate)) % nearer to diploid
			    y_ = [maxY/2 usedPlotCGH(i)*ploidy_estimate/2 usedPlotCGH(i)*ploidy_estimate/2 maxY/2];
			else
			    y_ = [maxY/2 usedPlotCGH(i)*ploidy_estimate/4 usedPlotCGH(i)*ploidy_estimate/4 maxY/2];
			end;
		    else
			y_ = [maxY/2 usedPlotCGH(i)*ploidy_estimate/2 usedPlotCGH(i)*ploidy_estimate/2 maxY/2];
		    end;
		end;
	    else % scale_type = Log2Ratio.
		if (ploidy_estimate == 0) % no value; no scaling.
		    y_ = [maxY/2 usedPlotCGH(i)+maxY/2 usedPlotCGH(i)+maxY/2 maxY/2];
		else
		    %y_ = [maxY/2 usedPlotCGH(i)+maxY/2 usedPlotCGH(i)+maxY/2 maxY/2];
		    y_ = [maxY/2 log2(pow2(usedPlotCGH(i))*ploidy_estimate/2)+maxY/2 ...
			log2(pow2(usedPlotCGH(i))*ploidy_estimate/2)+maxY/2 maxY/2];
		end;
	    end;
	    % makes a blackbar for each bin.
	    f = fill(x_,y_,c_);
	    set(f,'linestyle','none');
	end;
	x2 = chr_size(chromosome)*chr_length_scale_multiplier;
	if (strcmp(scale_type,'Ratio') == 1)
            plot([0; x2], [maxY/2;maxY/2],'Color',[0 0 0]);  % 2n line.
            line([0 x2], [maxY/4*3 maxY/4*3],'Color',[0.85 0.85 0.85]); % top line.
            line([0 x2], [maxY/4   maxY/4],  'Color',[0.85 0.85 0.85]); % bottom line.
	else
	    plot([0; x2], [maxY/2;maxY/2],'Color',[0 0 0]);  % 2n line.
	    line([0 x2], [maxY/2*log2(3) maxY/2*log2(3)], 'Color',[0.85 0.85 0.85]); % top cen.
	end;
	%end cgh plot section.

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
		fill([leftEnd  leftEnd  (leftEnd+dx) ],[(maxY-dy) maxY       maxY], colorNoData, 'EdgeColor', colorNoData);
		fill([rightEnd rightEnd (rightEnd-dx)],[(maxY-dy) maxY       maxY], colorNoData, 'EdgeColor', colorNoData);
		fill([leftEnd  leftEnd  (leftEnd+dx) ],[(0+dy)    0          0   ], colorNoData, 'EdgeColor', colorNoData);
		fill([rightEnd rightEnd (rightEnd-dx)],[(0+dy)    0          0   ], colorNoData, 'EdgeColor', colorNoData);
		fill([(x1-dx)  x1       (x1+dx)      ],[maxY      (maxY-dy)  maxY], colorNoData, 'EdgeColor', colorNoData);
		fill([(x1-dx)  x1       (x1+dx)      ],[0         dy         0   ], colorNoData, 'EdgeColor', colorNoData);
		plot([leftEnd leftEnd (leftEnd+dx) x1-dx x1      x2      x2+dx rightEnd-dx rightEnd rightEnd rightEnd-dx x2+dx x2 x1 x1-dx leftEnd+dx leftEnd], ...
			 [dy      maxY-dy maxY         maxY  maxY-dy maxY-dy maxY  maxY        maxY-dy  dy       0           0     dy dy 0     0          dy     ], ...
			 'Color',[0 0 0]);
	end;
	%end show centromere.

	%show MRS locations
	if (show_MRS)
	    plot([leftEnd rightEnd], [-1.5 -1.5],'color',[0 0 0]);
	    hold on;
	    for i = 1:length(MRS_location)
		if (MRS_chr(i) == chromosome)
		    MRSloc = MRS_location(i)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
		    plot(MRSloc,-1.5,'k:o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3)
		end;
	    end;
	    hold off;
	end;
	%end show MRS locations.

	%final reformatting.
	xlim([0,chr_size(chromosome)*chr_length_scale_multiplier]);
	if (show_MRS == true)
	    ylim([-1.5,maxY]);
	else
	    ylim([0,maxY]);
	end;
	set(gca,'YTick',[]);
	set(gca,'TickLength',[(Linear_TickSize*chr_size(1)/chr_size(chromosome)) 0]); %ensures same tick size on all subfigs.
	if (chromosome == 1)
	    ylabel(experiment_name, 'Rotation', 0, 'HorizontalAlign', 'right', 'VerticalAlign', 'bottom','Interpreter','none','FontSize',5);
	end;
	set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
	if (show_Xlabels == true)
	    set(gca,'XTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8','2.0','2.2','2.4','2.6','2.8','3.0','3.2'});
	else
	    set(gca,'XTickLabel',[]);
	end;
        if (chromosome == 1)
            if (CGH_Genomic_display == true)
	        if (strcmp(scale_type,'Ratio') == 1)
		    set(gca,'YTick',[Linear_maxY/4 Linear_maxY/2 Linear_maxY/4*3 Linear_maxY]);
		    if (Yscale_nearest_even_ploidy == true)
		        if (abs(2-ploidy_estimate) < abs(4-ploidy_estimate)) % nearer to diploid
			    set(gca,'YTickLabel',{'1','2','3','4'});
                        else % nearer to tetraploid
			    set(gca,'YTickLabel',{'2','4','6','8'});
                        end;
		    else
			set(gca,'YTickLabel',{'1','2','3','4'});
		    end;
	        else
		    set(gca,'YTick',[0 (Linear_maxY/2) Linear_maxY/2*log2(3) Linear_maxY]);
		    set(gca,'YTickLabel',{'1','2','3','4'});
		end;
	    end;
	end;
	set(gca,'FontSize',8);
	%end final reformatting.

	% shift back to main figure generation.
	figure(fig);
	hold on;
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
        if (ploidy_estimate ~= 0) % ploidy value given.
            text(0.35,2.0,['Ploidy = ' num2str(ploidy_estimate) 'n']);
        end;
    else
        if (ploidy_estimate ~= 0) % ploidy value given.
            text(0.35,2.5,['Ploidy = ' num2str(ploidy_estimate) 'n']);
        end;
    end;
else
    if (Show_Genomic_LOH_fraction == true)
        text(0.35,2.0,['LOH_{snp/cgh} = ' num2str(100*SNPs_hom/SNPs_total) '%']);
        if (ploidy_estimate ~= 0) % ploidy value given.
            text(0.35,1.5,['Ploidy = ' num2str(ploidy_estimate) 'n']);
        end;
    else
        if (ploidy_estimate ~= 0) % ploidy value given.
            text(0.35,1.5,['Ploidy = ' num2str(ploidy_estimate) 'n']);
        end;
    end;
end;

hold off;

fprintf('\n');
% save then close main genome figure.
set(fig,'PaperPosition',[0 0 8 6]);
saveas(fig, [workingDir '/fig_1.' experiment_name '.eps'], 'epsc');
saveas(fig, [workingDir '/fig_1.' experiment_name '.png'], 'png');
delete(fig);

% save then delete linear genome figure, after resizing image output.
set(Linear_fig,'PaperPosition',[0 0 8 0.62222222]*2);
saveas(Linear_fig, [workingDir '/fig_2.' experiment_name '.eps'], 'epsc');
saveas(Linear_fig, [workingDir '/fig_2.' experiment_name '.png'], 'png');
delete(Linear_fig);

%% ========================================================================
% end stuff
%==========================================================================
fprintf('\n');

% delete any hanging waitbar dialogs or already saved figures.
%set(0,'ShowHiddenHandles','on');
%delete(get(0,'Children'));

result_image_location1 = [workingDir '/fig_1.' experiment_name '.png'];
result_image_location2 = [workingDir '/fig_2.' experiment_name '.png'];
archive_data_location = {[workingDir '/' microarray_design '.' strrep(experiment_name,' ','_') '.SNP_data.mat'], ...
                         [workingDir '/' microarray_design '.' strrep(experiment_name,' ','_') '.CGH_data.mat'], ...
                         [workingDir '/' microarray_design '.' strrep(experiment_name,' ','_') '.datasetDetails.mat']};

new_fid = fopen([workingDir '/output2.txt'],'w');
fprintf(new_fid,'Input:\n');
fprintf(new_fid,['\tMicroarray design  = \t' microarray_design '\n']);
fprintf(new_fid,['\tRaw data dir       = \t' raw_data_dir '\n']);
fprintf(new_fid,['\tRaw data file      = \t' raw_data_file '\n']);
fprintf(new_fid,['\tHeader rows        = \t' num2str(header_rows) '\n']);
fprintf(new_fid,['\tProbe name col     = \t' num2str(probeName_col) '\n']);
fprintf(new_fid,['\tData col ch1       = \t' num2str(data_col_ch1) '\n']);
fprintf(new_fid,['\tData col ch2       = \t' num2str(data_col_ch2) '\n']);
fprintf(new_fid,['\tData col ratio     = \t' num2str(data_col_ratio) '\n']);
fprintf(new_fid,['\tData col log2ratio = \t' num2str(data_col_log2ratio) '\n']);
fprintf(new_fid,['\tPhasing dataset    = \t' phasing_dataset '\n']);
fprintf(new_fid,['\tPloidy estimate    = \t' num2str(ploidy_estimate) '\n']);
fprintf(new_fid,['\tImage format       = \t' image_format '\n']);
fprintf(new_fid,['\tExperiment name    = \t' experiment_name '\n']);
fprintf(new_fid,'Output:\n');
fprintf(new_fid,['\tImage output 1     = \t' result_image_location1 '\n']);
fprintf(new_fid,['\tImage output 2     = \t' result_image_location2 '\n']);
fprintf(new_fid,['\tArchvive file #1   = \t' archive_data_location{1} '\n']);
fprintf(new_fid,['\tArchvive file #2   = \t' archive_data_location{2} '\n']);
fprintf(new_fid,['\tArchvive file #3   = \t' archive_data_location{3} '\n']);
time_to_complete = toc;
fprintf(new_fid,['\tTime to process    = \t' num2str(time_to_complete) ' s']);
fclose(new_fid);

%% Let the processing pipeline know that the analysis has completed.
new_fid = fopen([workingDir '/complete.txt'],'w');
fprintf(new_fid,'complete');
fclose(new_fid);

end
