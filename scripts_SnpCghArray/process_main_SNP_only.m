function [result_image_location, archive_data_location] = ...
    process_main_SNP_only(microarray_design, data_file, header_rows, probeName_col, data_col_ch1, ...
                 data_col_ch2, data_col_ratio, data_col_log2ratio, phasing_dataset, ploidy_estimate, ploidyBase, image_format, ...
                 experiment_name, workingDir, show_MRS_string)
% PROCESS_MAIN wrapper around microarray analysis pipeline for integration into online tools and website.
% This version does not load calibration data and generate a hapmap...   see 'calibration_setup.m' for relevant script.
%

%%% Let the processing pipeline know that the analysis has started.
%new_fid = fopen([workingDir 'working.txt'],'a');
%fprintf(new_fid,'working');
%fclose(new_fid);

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
fprintf(['    11 : ploidy baseline      : ' ploidyBase        '\n']);
fprintf(['    12 : image format         : ' image_format       '\n']);
fprintf(['    13 : experiment name      : ' experiment_name    '\n']);
fprintf(['    14 : working directory    : ' workingDir         '\n']);
fprintf(['    15 : show MRS annotations : ' show_MRS_string    '\n']);

fprintf(['[process_main.m] : current folder.\n']);
	currentDir = pwd;
    fprintf(['    ' currentDir '\n']);

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
ploidyBase          = round(str2num(ploidyBase));
if (ploidyBase > 4);   ploidyBase = 4;   end;
if (ploidyBase < 1);   ploidyBase = 1;   end;
fprintf(['\nEuploid base = "' num2str(ploidyBase) '"\n']);


tic;
% =========================================================================
% Generate graphs involved in analysis of SNP probe data from 1' Agilent
% slide.
%==========================================================================
% Displays SNP and CGH data in a common figure format.
%    colorBars                  : Draws SNP information as background colors on chromosomes.
%    blendColorBars             : Blends adjacent colorBars to remove high contrast edges.
%    infillColorBars            : Fills in "no-data" areas with adjacent colored areas.
%    AnglePlot                  : Produces angle histogram of scatterplots, to left of chromosome.
%       FillColors              : Fills angle histogram with homolog identity colors.
%    ChrNum                     : Draws a large numeral of copy# estimage to right of chromosome.
%    Yscale_nearest_even_ploidy : automatically alters y-axis from [0..4] to [0..8] if ploidy is ~tetraploid.
%    Show_Genomic_LOH_fraction  : Shows an automatically determined %LOH, from how many SNP probes are homozygous.
%    show_unnassigned           : Shows unnassigned SNP pair data in histogram/chromosome plots.
%    Centromere_format          : Controls how centromeres are depicted.   [0..2]   '2' is pinched cartoon default.
%    bases_per_bin              : Controls bin sizes for SNP/CGH fractions of plot.
%    scale_type                 : 'Ratio' or 'Log2Ratio' y-axis scaling of copy number.
%                                 'Log2Ratio' does not properly scale CGH data by ploidy.
%    Chr_max_width              : max width of chromosomes as fraction of figure width.
SNP_Genomic_display                  = true;
	Linear_display               = true;
	colorBars                    = true;
	blendColorBars               = false;
	infillColorBars              = false;
	AnglePlot                    = true;
        FillColors                   = true;
	ChrNum                       = false;
	Yscale_nearest_even_ploidy   = true;
	Show_Genomic_LOH_fraction    = true;
	show_unnassigned             = false;
	if strcmp(show_MRS_string,'1')
		show_MRS             = true;
	else
		show_MRS             = false;
	end;
	Centromere_format            = 2;
	scale_type                   = 'Ratio';
	Chr_max_width                = 0.8;
	Show_ChARM_edges             = true;

    %not in place yet... uninformative currently looks like unnassigned to script.
    show_uninformative           = false;

axisLabelPosition_vert  = 0.01125;
axisLabelPosition_horiz = 0.01125;

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
% Control variables for Candida albicans SC5314.
%--------------------------------------------------------------------------
% Determines for which chromosomes to construct figures from.
chromosomes_to_analyze = 1:num_chrs;

% Defines chromosome sizes in bp. (diploid total=28,567,7888)
% Defines centromere locations in bp.
% Defines MRS locations in bp.
[centromeres, chr_sizes,MRSs,rDNA] = Load_information();
for i = chromosomes_to_analyze
	chr_size(i)  = chr_sizes(i).size;
	cen_start(i) = centromeres(i).start;
	cen_end(i)   = centromeres(i).end;
end;
for i = 1:length(MRSs)
	MRS_chr(i)      = MRSs(i).chr;
	MRS_location(i) = (MRSs(i).start+MRSs(i).end)/2;
	MRS_color_1(i)  = MRSs(i).color_fill;
	MRS_color_2(i)  = MRSs(i).color_edge;
end;
rDNA_chr            = rDNA.chr;
rDNA_location       = (rDNA.start+rDNA.end)/2;
rDNA_color_1        = rDNA.color_fill;
rDNA_color_2        = rDNA.color_edge;
clear centromeres chr_sizes MRSs rDNA;


% basic plot parameters not defined per genome.
TickSize         = -0.005;  %negative for outside, percentage of longest chr figure.
bases_per_bin    = max(chr_size)/700;
maxY             = ploidyBase*2;
cen_tel_Xindent  = 5; 
cen_tel_Yindent  = maxY/5;


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

%% ====================================================================
% Process raw data file into structures containing CGH and SNP information (if not already done).
%----------------------------------------------------------------------
data_file_load_online(matlab_save_dir,raw_data_dir,raw_data_file, microarray_design, experiment_name);
% load experimental values for all the probes.
load([matlab_save_dir '/' experiment_name '.' microarray_design '.SNP_data.mat']);
SNP_probeset_length = length(probeset1);

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
% Load cutoffs for experimental dataset chromosome segments
%----------------------------------------------------------------------
fprintf('\nLoading previously calculated SNP interpretation cutoffs for microarray.');
load([workingDir 'SNP_cutoffs.mat']);
fprintf('\nLoading previously calculated datasetDetails structure containing CNV segmentation info.');
load([matlab_save_dir '/' experiment_name '.' microarray_design '.datasetDetails.mat']);
fprintf('\nLoading previously processed SNP data.');
load([workingDir 'SNPdata_collected.mat']);

    
%% ========================================================================
% Plotting probes across genome by interpretation catagory after accounting
% for polarity of SNP pairs.
%    CGH info plotted as histogram.
%    SNP info plotted as colorbar.
%==========================================================================
fprintf(['\nGenerating SNP figures from: ' strrep(raw_data_dir,'\','\\')]);
% define number of and labels for chromosomes.
chr_labels = {'Chr1','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7','ChrR'};

bases_per_bin = max(chr_size)/700;
chr_length_scale_multiplier = 1/bases_per_bin;

% Initialize linear view figure plot parameters.
if (Linear_display == true)
	Linear_fig           = figure(2);
	Linear_genome_size   = sum(chr_size);
	Linear_Chr_max_width = 0.91;
	Linear_left_start    = 0.01;
	Linear_left_chr_gap  = 0.07/(num_chrs);
	Linear_height        = 0.6;
	Linear_base          = 0.1;
	Linear_TickSize      = -0.01;  %negative for outside, percentage of longest chr figure.
	maxY                 = ploidyBase*2;
	Linear_left          = Linear_left_start;
end;


%% ====================================================================
% Define colors for use in figures.
%----------------------------------------------------------------------
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

% basic plot parameters.
left            = 0.15;
height          = 0.5/num_chrs;
base            = 0.1;
vertical_margin = 0.3/num_chrs;


%% ====================================================================
% Main figure drawing, including CNV data and SNP data.
%----------------------------------------------------------------------

% main figure calculation and generation.
fig = figure(1);
set(gcf, 'Position', [0 70 1024 600]);
for chr = 1:num_chrs
	usedPlot0   = chr_SNPdata{chr,1};  % 'ab'
	usedPlot1   = chr_SNPdata{chr,2};  % 'a'/'aa'/'aaa'/'aaaa'/'aaaaa'/'aaaaaa'
	usedPlot2   = chr_SNPdata{chr,3};  % 'b'/'bb'/'bbb'/'bbbb'/'bbbbb'/'bbbbbb'
	usedPlot3   = chr_SNPdata{chr,4};  % 'aab'
	usedPlot4   = chr_SNPdata{chr,5};  % 'abb'
	usedPlot5   = chr_SNPdata{chr,6};  % 'aaab'
	usedPlot6   = chr_SNPdata{chr,7};  % 'abbb'
	usedPlot7   = chr_SNPdata{chr,8};  % 'aaaab'
	usedPlot8   = chr_SNPdata{chr,9};  % 'aaabb'
	usedPlot9   = chr_SNPdata{chr,10}; % 'aabbb'
	usedPlot10  = chr_SNPdata{chr,11}; % 'abbbb'
	usedPlot11  = chr_SNPdata{chr,12}; % 'aaaaab'
	usedPlot12  = chr_SNPdata{chr,13}; % 'abbbbb'
	usedPlot13  = chr_SNPdata{chr,14}; % hom unassigned.

	% places chr copy number to the right of the main chr cartoons.
	if (ChrNum == true)
		% subplot to show chromosome copy number value.
		width = 0.020;
		if (chr == 8)
			bottom = base + (chr)*(height+vertical_margin);
		else
			bottom = base + (8-chr)*(height+vertical_margin);
		end;

		subplot('Position',[(0.15 + Chr_max_width*chr_size(chr)/max(chr_size) + 0.005) bottom width height]);
		axis off square;
		set(gca,'YTick',[]);
		set(gca,'XTick',[]);
		if (length(datasetDetails.chrCopyNum{chr}) == 1)
			chr_string = num2str(datasetDetails.chrCopyNum{chr}(1));
		else
			chr_string = num2str(datasetDetails.chrCopyNum{chr}(1));
			for i = 2:length(datasetDetails.chrCopyNum{chr})
				chr_string = [chr_string ',' num2str(datasetDetails.chrCopyNum{chr}(i))];
			end;
		end;
		text(0.1,0.5, chr_string,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',24 );
	end;

	% make allelic fraction histograms to the left of the main chr cartoons.
	if (AnglePlot == true)
		width = 0.075;
		if (chr == 8);  bottom = base + (chr)*(height+vertical_margin);
		else                   bottom = base + (8-chr)*(height+vertical_margin);
		end;
		height = 0.5/num_chrs;

		for segment = 1:length(datasetDetails.chrCopyNum{chr})
			if (segment == 1)
				subplot('Position',[0.03 bottom width (height/length(datasetDetails.chrCopyNum{chr}))]);
			else
				subplot('Position',[0.03 (bottom+height/length(datasetDetails.chrCopyNum{chr})*(segment-1)) width (height/length(datasetDetails.chrCopyNum{chr}))]);
			end;
			histAll = [];
			histAll2 = [];
			for i = 1:2:SNP_probeset_length
				if (probeset1(i).probe_chromosome == chr) && ...
					(probeset1(i).probe_location >= round(1+chr_size(chr)*datasetDetails.chr_breaks{chr}(segment))) && ...
						(probeset1(i).probe_location <= round(chr_size(chr)*datasetDetails.chr_breaks{chr}(segment+1)))
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

			% Invert histogram values;
			histAll = 1-histAll;

			smoothed = smooth_gaussian(hist(histAll,200),3,20);
			% make a smoothed version of just the endpoints used to ensure histogram bounds.
			histAll2(1) = 0;
			histAll2(2) = 1;
			smoothed2 = smooth_gaussian(hist(histAll2,200),3,20);
			% subtract the smoothed endpoints from the histogram to remove the influence of the added endpoints.
			smoothed = (smoothed-smoothed2);
			smoothed = smoothed/max(smoothed);

			hold on;
			copynum    = datasetDetails.chrCopyNum{chr}(segment);
			region_    = 0;

			fprintf(['important Gaussians = ' num2str(datasetDetails.importantGaussians{chr,segment}) '\n']);
			fprintf(['Gaussian Cutoffs    = ' num2str(datasetDetails.cutoffs{chr,segment}) '\n']);

			for region = datasetDetails.importantGaussians{chr,segment}
				if (region_ < length(datasetDetails.importantGaussians{chr,segment}))
					region_ = region_+1;
				else
					break;
				end;

				if (FillColors == true)
					if (show_uncalibrated == true)
						color = colorUn_Hom;
					else
						if (copynum <= 1.5) %monosomy
							if (region == 1); color = colorA;
							else              color = colorB;
							end;
							if (segment == 1)
								set(gca,'XTick',[0 200]);
								set(gca,'XTickLabel',{'a','b'});
							end;
						elseif (copynum <= 2.5) %disomy
							if (region == 1);     color = colorA;
							elseif (region == 2); color = colorAB;
							else                  color = colorB;
							end;
							if (segment == 1)
								set(gca,'XTick',0:100:200);
								set(gca,'XTickLabel',{'a','ab','b'});
							end;
						elseif (copynum <= 3.5) %trisomy
							if (region == 1);     color = colorA;
							elseif (region == 2); color = colorAAB;
							elseif (region == 3); color = colorABB;
							else                  color = colorB;
							end;
							if (segment == 1)
								set(gca,'XTick',[0 66.667 133.333 200]);
								set(gca,'XTickLabel',{'a','aab','abb','b'});
							end;
						elseif (copynum <= 4.5) %tetrasomy
							if (region == 1);     color = colorA;
							elseif (region == 2); color = colorAAAB;
							elseif (region == 3); color = colorAB;
							elseif (region == 4); color = colorABBB;
							else                  color = colorB;
							end;
							if (segment == 1)
								set(gca,'XTick',0:50:200);
								set(gca,'XTickLabel',{'a', '3ab', 'ab', 'a3b' 'b'});
							end;
						elseif (copynum <= 5.5) %pentasomy
							if (region == 1);     color = colorA;
							elseif (region == 2); color = colorAAAAB;
							elseif (region == 3); color = colorAAABB;
							elseif (region == 4); color = colorAABBB;
							elseif (region == 5); color = colorABBBB;
							else                  color = colorB;
							end;
							if (segment == 1)
								set(gca,'XTick',0:40:200);
								set(gca,'XTickLabel',{'a', '4ab', '3abb', 'aa3b', 'a4b' 'b'});
							end;
						else % if (copynum <= 6.5) %hexasomy
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
								set(gca,'XTickLabel',{'a', '5ab', 'aab', 'ab', 'abb', 'a5b' 'b'});
							end;
						end;
					end;
				else
					color = colorAB;
				end;
				if (length(datasetDetails.importantGaussians{chr,segment}) == 1)
					area(1:200,smoothed(1:200),'FaceColor',color,'EdgeColor',color);
				else
					if (region_ > 1) && (region_ < length(datasetDetails.importantGaussians{chr,segment}))
						if (round(200-datasetDetails.cutoffs{chr,segment}(region_)*200) == 0)
							% If the difference between the position and the end is < 0, move the "region_" to the end of the list.
							last_region_ = region_-1;
							region_      = length(datasetDetails.importantGaussians{chr,segment});
						else
							last_region_ = region_-1;
						end;
					else
						last_region_ = region_-1;
					end;

					fprintf(['test :\t' num2str(region) '\t' num2str(region_) '\t']);
					if (region_ < length(datasetDetails.cutoffs{chr,segment}))
						fprintf([num2str(datasetDetails.cutoffs{chr,segment}(region_)) '\t' ...
						         num2str(round(200-datasetDetails.cutoffs{chr,segment}(region_)*200)) '\t' ...
						         num2str(smoothed(round(200-datasetDetails.cutoffs{chr,segment}(region_)*200))) '\t' ...
						         num2str(last_region_) '\n']);
					else
						fprintf([num2str(last_region_) '\n']);
					end;

					if (region_ == 1)
						area(		  round(200-datasetDetails.cutoffs{chr,segment}(region_)*200):200, ...
							 smoothed(round(200-datasetDetails.cutoffs{chr,segment}(region_)*200):200) ...
							 ,'FaceColor',color,'EdgeColor',color);
					elseif (region_ == length(datasetDetails.importantGaussians{chr,segment}))
						area(		  1:round(200-datasetDetails.cutoffs{chr,segment}(last_region_)*200), ...
							 smoothed(1:round(200-datasetDetails.cutoffs{chr,segment}(last_region_)*200)) ...
							 ,'FaceColor',color,'EdgeColor',color);
					else
						area(		  round(200-datasetDetails.cutoffs{chr,segment}(region_     )*200): ...
									  round(200-datasetDetails.cutoffs{chr,segment}(last_region_)*200), ...
							 smoothed(round(200-datasetDetails.cutoffs{chr,segment}(region_     )*200): ...
									  round(200-datasetDetails.cutoffs{chr,segment}(last_region_)*200)) ...
							,'FaceColor',color,'EdgeColor',color);
					end;
				end;
			end;
			for peak = 1:length(datasetDetails.peaks{chr,segment})
				plot([200-datasetDetails.peaks{chr,segment}(peak)*200; 200-datasetDetails.peaks{chr,segment}(peak)*200],[0; 1],'color',colorPeak);
			end;
			for cutoff = 1:length(datasetDetails.cutoffs{chr,segment})
				plot([200-datasetDetails.cutoffs{chr,segment}(cutoff)*200; 200-datasetDetails.cutoffs{chr,segment}(cutoff)*200],[0; 1],'color',colorCutoff);
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
	width = Chr_max_width*chr_size(chr)/max(chr_size);
	if (chr == 8)
		bottom = base + (chr)*(height+vertical_margin);
	else
		bottom = base + (8-chr)*(height+vertical_margin);
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
		% standard : draw colorbars.
		for chr_bin = 1:length(usedPlot0)+1;
			x_ = [chr_bin chr_bin chr_bin-1 chr_bin-1];
			y_ = [0 maxY maxY 0];
			c_post(1) = colors(chr_bin,1);
			c_post(2) = colors(chr_bin,2);
			c_post(3) = colors(chr_bin,3);
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

	% standard : show centromere.
	x1 = cen_start(chr)*chr_length_scale_multiplier;
	x2 = cen_end(chr)*chr_length_scale_multiplier;
	leftEnd  = 0.5*(5000/bases_per_bin);
	rightEnd = chr_size(chr)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
	if (Centromere_format == 0)
		h = fill([x1-maxY/2 x1 x2 x2+maxY/2], [0 maxY/4 maxY/4 0], [1 1 0]);
		set(h,'FaceAlpha',0.5);
		h = fill([x1-maxY/2 x1 x2 x2+maxY/2], [maxY maxY*3/4 maxY*3/4 maxY], [1 1 0]);
		set(h,'FaceAlpha',0.5);
	elseif (Centromere_format == 1)
		h = fill([x1-maxY/2 x1 x2 x2+maxY/2 x2 x1], [maxY/2 maxY*3/4 maxY*3/4 maxY/2 maxY/4 maxY/4], [1 1 0]);
		set(h,'FaceAlpha',0.5);
	elseif (Centromere_format == 2)
		dx = cen_tel_Xindent; %5*(5000/bases_per_bin);
		dy = cen_tel_Yindent; %maxY/5;
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
	% end show centromere.

	% standard : show MRS locations
	if (show_MRS)
		plot([leftEnd rightEnd], [-maxY/10*1.5 -maxY/10*1.5],'color',[0 0 0]);
		hold on;
		for i = 1:length(MRS_location)
			if (MRS_chr(i) == chr)
				MRSloc = MRS_location(i)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
				plot(MRSloc,-maxY/10*1.5,'k:o','MarkerEdgeColor',MRS_color_2(i),'MarkerFaceColor',MRS_color_1(i),'MarkerSize',5);
			end;
		end;
		hold off;
	end;
	% end show MRS locations.

	% standard : show rDNA location.
	if (show_MRS)
		hold on;
		if (rDNA_chr == chr)
			rDNAloc = rDNA_location*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
			plot(rDNAloc,-maxY/10*1.5,'k:o','MarkerEdgeColor',rDNA_color_2,'MarkerFaceColor',rDNA_color_1,'MarkerSize',5);
		end;
		hold off;
	end;
	% end show rDNA location.

	% standard : show ChARM edges.
	if ((Show_ChARM_edges) && (show_MRS))
		for i = 1:length(segmental_aneuploidy)
			if (segmental_aneuploidy(i).chr == chr)
				EDGEloc = segmental_aneuploidy(i).break*chr_size(chr)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
				plot([EDGEloc EDGEloc], [(-maxY/10*2.5) 0],  'Color',[1 0 0],'LineWidth',2);
			end;
		end;
	end;
	% end show ChARM edges.

	hold off;
	xlim([0,chr_size(chr)*chr_length_scale_multiplier]);
	if (show_MRS == true)
		ylim([-maxY/10*1.5,maxY]);
	else
		ylim([0,maxY]);
	end;
	set(gca,'YTick',[]);
	set(gca,'YTickLabel',[]);
	set(gca,'TickLength',[(TickSize*chr_size(1)/chr_size(chr)) 0]); %ensures same tick size on all subfigs.
	text(-50000/5000/2*3, maxY*3/2,     chr_labels(chr), 'Rotation',90, 'HorizontalAlignment','center', 'VerticalAlign','bottom', 'Fontsize',20);
	set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
	set(gca,'XTickLabel',[]);
	set(gca,'FontSize',8);
	if (chr == 8)
		title(experiment_name,'Interpreter','none','FontSize',12);
	end;

	%% Linear figure draw section
	if (Linear_display == true)
		figure(Linear_fig);
		Linear_width = Linear_Chr_max_width*chr_size(chr)/Linear_genome_size;
		subplot('Position',[Linear_left Linear_base Linear_width Linear_height]);
		Linear_left = Linear_left + Linear_width + Linear_left_chr_gap;
		hold on;
		title(chr_labels(chr),'Interpreter','none','FontSize',20);

		% linear : display color bars for this chr.
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

		% linear : show centromere.
		x1 = cen_start(chr)*chr_length_scale_multiplier;
		x2 = cen_end(chr)*chr_length_scale_multiplier;
		leftEnd  = 0.5*(5000/bases_per_bin);
		rightEnd = chr_size(chr)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
		if (Centromere_format == 0)
			h = fill([x1-maxY/2 x1 x2 x2+maxY/2], [0 maxY/4 maxY/4 0], [1 1 0]);
			set(h,'FaceAlpha',0.5);
			h = fill([x1-maxY/2 x1 x2 x2+maxY/2], [maxY maxY*3/4 maxY*3/4 maxY], [1 1 0]);
			set(h,'FaceAlpha',0.5);
		elseif (Centromere_format == 1)
			h = fill([x1-maxY/2 x1 x2 x2+maxY/2 x2 x1], [maxY/2 maxY*3/4 maxY*3/4 maxY/2 maxY/4 maxY/4], [1 1 0]);
			set(h,'FaceAlpha',0.5);
		elseif (Centromere_format == 2)
			dx = cen_tel_Xindent; %5*(5000/bases_per_bin);
			dy = cen_tel_Yindent; %maxY/10;
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
		% end show centromere.

		% linear : show MRS locations
		if (show_MRS)
			plot([leftEnd rightEnd], [-maxY/10*1.5 -maxY/10*1.5],'color',[0 0 0]);
			hold on;
			for i = 1:length(MRS_location)
				if (MRS_chr(i) == chr)
					MRSloc = MRS_location(i)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
					plot(MRSloc,-maxY/10*1.5,'k:o','MarkerEdgeColor',MRS_color_2(i),'MarkerFaceColor',MRS_color_1(i),'MarkerSize',5);
				end;
			end;
			hold off;
		end;
		% end show MRS locations.

		% linear : show rDNA location.
		if (show_MRS)
			hold on;
			if (rDNA_chr == chr)
				rDNAloc = rDNA_location*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
				plot(rDNAloc,-maxY/10*1.5,'k:o','MarkerEdgeColor',rDNA_color_2,'MarkerFaceColor',rDNA_color_1,'MarkerSize',5);
			end;
			hold off;
		end;
		% end show rDNA location.

		% linear : final reformatting.
		xlim([0,chr_size(chr)*chr_length_scale_multiplier]);
		if (show_MRS == true)
			ylim([-maxY/10*1.5,maxY]);
		else
			ylim([0,maxY]);
		end;
		set(gca,'YTick',[]);
		set(gca,'YTickLabel',[]);
		set(gca,'TickLength',[(Linear_TickSize*chr_size(1)/chr_size(chr)) 0]); %ensures same tick size on all subfigs.
		set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
		set(gca,'XTickLabel',[]);
		set(gca,'FontSize',8);
		%end final reformatting.

		% shift back to main figure generation.
		hold off;
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
set(   fig,        'PaperPosition',[0 0 8 6]*2);
saveas(fig,        [workingDir 'fig.SNP-map.1.eps'], 'epsc');
saveas(fig,        [workingDir 'fig.SNP-map.1.png'], 'png' );
delete(fig);
set(   Linear_fig, 'PaperPosition',[0 0 8 0.62222222]*2);
saveas(Linear_fig, [workingDir 'fig.SNP-map.2.eps'], 'epsc');
saveas(Linear_fig, [workingDir 'fig.SNP-map.2.png'], 'png' );
delete(Linear_fig);


%% ========================================================================
% End stuff
%==========================================================================
fprintf('\n');

% delete any hanging waitbar dialogs or already saved figures.
%set(0,'ShowHiddenHandles','on');
%delete(get(0,'Children'));

result_image_location1 = [workingDir 'fig_1.' experiment_name '.png'];
result_image_location2 = [workingDir 'fig_2.' experiment_name '.png'];
archive_data_location = {	[workingDir microarray_design '.' strrep(experiment_name,' ','_') '.SNP_data.mat'], ...
							[workingDir microarray_design '.' strrep(experiment_name,' ','_') '.CGH_data.mat'], ...
							[workingDir microarray_design '.' strrep(experiment_name,' ','_') '.datasetDetails.mat']};

%% Let the processing pipeline know that the analysis has completed.
new_fid        = fopen([workingDir 'complete.txt'],'w');
completionTime = datestr(clock, 0)
fprintf(new_fid,completionTime);
fclose(new_fid);

delete([workingDir 'working.txt']);

end
