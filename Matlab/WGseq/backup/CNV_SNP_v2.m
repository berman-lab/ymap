function [] = CNV_SNP_v2(main_dir,user,genomeUser,projectChild,projectParent,genome,ploidyEstimateString,ploidyBaseString, ...
                         CNV_verString,SNP_verString,LOH_verString,displayBREAKS, referenceCHR);
%% ========================================================================
% Generate graphs involved in analysis of CNV and SNP information from NGS
% datasets.
%==========================================================================
%    Centromere_format          : Controls how centromeres are depicted.   [0..2]   '2' is pinched cartoon default.
%    bases_per_bin              : Controls bin sizes for CGH fractions of plot.
%    scale_type                 : 'Ratio' or 'Log2Ratio' y-axis scaling of copy number.
%                                 'Log2Ratio' does not properly scale CGH data by ploidy.
%    Chr_max_width              : max width of chrs as fraction of figure width.
Centromere_format_default   = 0;
Yscale_nearest_even_ploidy  = true;
HistPlot                    = true;
ChrNum                      = true;
Chr_max_width               = 0.8;
show_annotations            = true;
analyze_rDNA                = true;
colorBars                   = true;
blendColorBars              = false;
Linear_display              = true;
Low_quality_ploidy_estimate = true;


%%=========================================================================
% Load FASTA file name from 'reference.txt' file for project.
%--------------------------------------------------------------------------
userReference    = [main_dir 'users/' user '/genomes/' genome '/reference.txt'];
defaultReference = [main_dir 'users/default/genomes/' genome '/reference.txt'];
if (exist(userReference,'file') == 0)   
	FASTA_string = strtrim(fileread(defaultReference));
else                    
	FASTA_string = strtrim(fileread(userReference));
end;
[FastaPath,FastaName,FastaExt] = fileparts(FASTA_string);


%%=========================================================================
% Control variables.
%--------------------------------------------------------------------------
projectChildDir  = [main_dir 'users/' user '/projects/' projectChild '/'];
if (exist([[main_dir 'users/default/projects/' projectParent '/']],'dir') == 7)
	projectParentDir = [main_dir 'users/default/projects/' projectParent '/'];
else  
	projectParentDir = [main_dir 'users/' user '/projects/' projectParent '/'];
end;
genomeDir        = [main_dir 'users/' genomeUser '/genomes/' genome '/'];

fprintf(['\n$$ projectParentDir : ' projectParentDir '\n$$ projectChildDir : ' projectChildDir '\n$$ genomeDir  : ' genomeDir '\n']);
[centromeres, chr_sizes, figure_details, annotations, ploidy_default] = Load_genome_information_1(projectChildDir,genomeDir, genome);
[Aneuploidy]                                                          = Load_dataset_information_1(projectChildDir, projectChild);

for i = 1:length(chr_sizes)
    chr_size(i) = 0;
end;
for i = 1:length(chr_sizes)
    chr_size(chr_sizes(i).chr)    = chr_sizes(i).size;
end;
for i = 1:length(centromeres)
    cen_start(centromeres(i).chr) = centromeres(i).start;
    cen_end(centromeres(i).chr)   = centromeres(i).end;
end;
if (length(annotations) > 0)
    fprintf(['\nAnnotations for ' genome '.\n']);
    for i = 1:length(annotations)
        annotation_chr(i)       = annotations(i).chr;
        annotation_type{i}      = annotations(i).type;
        annotation_start(i)     = annotations(i).start;
        annotation_end(i)       = annotations(i).end;
        annotation_fillcolor{i} = annotations(i).fillcolor;
        annotation_edgecolor{i} = annotations(i).edgecolor;
        annotation_size(i)      = annotations(i).size;
        fprintf(['\t[' num2str(annotations(i).chr) ':' annotations(i).type ':' num2str(annotations(i).start) ':' num2str(annotations(i).end) ':' annotations(i).fillcolor ':' annotations(i).edgecolor ':' num2str(annotations(i).size) ']\n']);
    end;
end;
for i = 1:length(figure_details)
    if (figure_details(i).chr == 0)
        key_posX   = figure_details(i).posX;
        key_posY   = figure_details(i).posY;
        key_width  = figure_details(i).width;
        key_height = figure_details(i).height;
    else
        chr_id    (figure_details(i).chr) = figure_details(i).chr;
        chr_label {figure_details(i).chr} = figure_details(i).label;
        chr_name  {figure_details(i).chr} = figure_details(i).name;
        chr_posX  (figure_details(i).chr) = figure_details(i).posX;
        chr_posY  (figure_details(i).chr) = figure_details(i).posY;
        chr_width (figure_details(i).chr) = figure_details(i).width;
        chr_height(figure_details(i).chr) = figure_details(i).height;
        chr_in_use(figure_details(i).chr) = str2num(figure_details(i).useChr);
    end;
end;
num_chrs = length(chr_size);


%%=========================================================================
%%= No further control variables below. ===================================
%%=========================================================================

% Sanitize user input of euploid state.
ploidyBase = round(str2num(ploidyBaseString));
if (ploidyBase > 4);   ploidyBase = 4;   end;
if (ploidyBase < 1);   ploidyBase = 1;   end;
fprintf(['\nEuploid base = "' num2str(ploidyBase) '"\n']);


% basic plot parameters not defined per genome.
TickSize         = -0.005;  %negative for outside, percentage of longest chr figure.
bases_per_bin    = max(chr_size)/700;
maxY             = ploidyBase*2;
cen_tel_Xindent  = 5;
cen_tel_Yindent  = maxY/5;


if (strcmp(projectParent,projectChild) == 1)
    fprintf(['\nGenerating CNV-SNP map figure from ' projectChild ' genome data.\n']);
else
    fprintf(['\nGenerating CNV-LOH map figure from ' projectParent '(parent) and ' projectChild '(child) genome data.\n']);
end;


%%=========================================================================
%% Load GC-bias corrected CGH data.
%%-------------------------------------------------------------------------
load([projectChildDir 'Common_CNV.mat']);		% 'CNVplot2','genome_CNV'


%% ====================================================================
% Load CNV edges from ChARM algorithm.
%----------------------------------------------------------------------
[segmental_aneuploidy] = Load_dataset_information_1(projectChildDir,projectChild);
segmental_aneuploidy = [];


%%================================================================================================
% Load SNP/LOH data.
%-------------------------------------------------------------------------------------------------
load([projectChildDir ['SNP_' SNP_verString '.mat']]);


%define colors for colorBars plot
colorNoData = [1.0   1.0   1.0  ]; %used when no data is available for the bin.
colorInit   = [0.5   0.5   0.5  ]; %external; used in blending at ends of chr.
%colorHET    = [0.0   0.0   0.0  ]; % near 1:1 ratio SNPs
colorHET    = [0.333 0.333 0.333];
%colorHET    = [0.667 0.667 0.667]; % near 1:1 ratio SNPs.
colorOddHET = [0.0   1.0   0.0  ]; % Het, but not near 1:1 ratio SNPs.
colorHOM    = [1.0   0.0   0.0  ]; % Hom SNPs;

%% ====================================================================
% Apply GC bias correction to the SNP data.
%   Amount of putative SNPs per standard bin vs GCbias per standard bin.
%----------------------------------------------------------------------
% Load standard bin GC_bias data from : standard_bins.GC_ratios.txt
fprintf(['standard_bins_GC_ratios_file :\n\t' main_dir 'users/' genomeUser '/genomes/' genome '/' FastaName '.GC_ratios.standard_bins.txt\n']);
standard_bins_GC_ratios_fid = fopen([main_dir 'users/' genomeUser '/genomes/' genome '/' FastaName '.GC_ratios.standard_bins.txt'], 'r');
fprintf(['\t' num2str(standard_bins_GC_ratios_fid) '\n']);
lines_analyzed = 0;
for chr = 1:num_chrs
    if (chr_in_use(chr) == 1)
        chr_GCratioData{chr} = zeros(1,ceil(chr_size(chr)/bases_per_bin));
    end;
end;
while not (feof(standard_bins_GC_ratios_fid))
    dataLine = fgetl(standard_bins_GC_ratios_fid);
    if (length(dataLine) > 0)
        if (dataLine(1) ~= '#')
            lines_analyzed = lines_analyzed+1;
            chr            = str2num(sscanf(dataLine, '%s',1));
            fragment_start = sscanf(dataLine, '%s',2);  for i = 1:size(sscanf(dataLine,'%s',1),2);      fragment_start(1) = []; end;    fragment_start = str2num(fragment_start);
            fragment_end   = sscanf(dataLine, '%s',3);  for i = 1:size(sscanf(dataLine,'%s',2),2);      fragment_end(1) = [];   end;    fragment_end   = str2num(fragment_end);
            GCratio        = sscanf(dataLine, '%s',4);  for i = 1:size(sscanf(dataLine,'%s',3),2);      GCratio(1) = [];        end;    GCratio        = str2num(GCratio);
            position       = ceil(fragment_start/bases_per_bin);
            chr_GCratioData{chr}(position) = GCratio;
        end;
    end;
end;
fclose(standard_bins_GC_ratios_fid);

% Gather SNP data into Totplot arrays;
for chr = 1:num_chrs
    if (chr_in_use(chr) == 1)
        TOTplot{chr} = chr_SNPdata{chr,1}+chr_SNPdata{chr,2}+chr_SNPdata{chr,3};
    end;
end;

% Gather SNP and GCratio data for LOWESS fitting.
SNPdata_all          = [];
GCratioData_all      = [];
for chr = 1:num_chrs
    if (chr_in_use(chr) == 1)
        SNPdata_all      = [SNPdata_all     TOTplot{chr}        ];
        GCratioData_all  = [GCratioData_all chr_GCratioData{chr}];
    end;
end;
medianRawY = median(SNPdata_all);
fprintf(['medianRawY = ' num2str(medianRawY) '\n']);

%% Clean up data by:
%%    deleting GC ratio data near zero.
%%    deleting CGH data beyond 3* the median value.  (rDNA, etc.)
SNPdata_clean                                          = SNPdata_all;
GCratioData_clean                                      = GCratioData_all;
SNPdata_clean(GCratioData_clean < 0.01)                = [];
GCratioData_clean(GCratioData_clean < 0.01)            = [];
GCratioData_clean(SNPdata_clean > max(medianRawY*3,3)) = [];
SNPdata_clean(SNPdata_clean > max(medianRawY*3,3))     = [];
GCratioData_clean(SNPdata_clean == 0)                  = [];
SNPdata_clean(SNPdata_clean == 0)                      = [];

% Perform LOWESS fitting.
rawData_X1        = GCratioData_clean;
rawData_Y1        = SNPdata_clean;
fprintf(['Lowess X:Y size : [' num2str(size(rawData_X1,1)) ',' num2str(size(rawData_X1,2)) ']:[' num2str(size(rawData_Y1,1)) ',' num2str(size(rawData_Y1,2)) ']\n']);
[fitX1, fitY1]    = optimize_mylowess_SNP(rawData_X1,rawData_Y1);

% Correct data using normalization to LOWESS fitting
Y_target = 1;
for chr = 1:num_chrs
    if (chr_in_use(chr) == 1)
        rawData_chr_X{chr}           = chr_GCratioData{chr};

        rawDataAll_chr_Y{chr}        = TOTplot{chr};
        fitDataAll_chr_Y{chr}        = interp1(fitX1,fitY1,rawData_chr_X{chr},'spline');
        normalizedDataAll_chr_Y{chr} = rawDataAll_chr_Y{chr}./fitDataAll_chr_Y{chr}*Y_target;

        rawData_chr_Y{chr,1}         = chr_SNPdata{chr,1};
        rawData_chr_Y{chr,2}         = chr_SNPdata{chr,2};
        rawData_chr_Y{chr,3}         = chr_SNPdata{chr,3};
        fitData_chr_Y{chr}           = interp1(fitX1,fitY1,rawData_chr_X{chr},'spline');
        normalizedData_chr_Y{chr,1}  = rawData_chr_Y{chr,1}./fitData_chr_Y{chr}*Y_target;
        normalizedData_chr_Y{chr,2}  = rawData_chr_Y{chr,2}./fitData_chr_Y{chr}*Y_target;
        normalizedData_chr_Y{chr,3}  = rawData_chr_Y{chr,3}./fitData_chr_Y{chr}*Y_target;
    end;
end;

% Move LOWESS-normalizd SNP data back into display pipeline.
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		chr_SNPdata{chr,1} = rawData_chr_Y{chr,1};
		chr_SNPdata{chr,2} = rawData_chr_Y{chr,2};
		chr_SNPdata{chr,3} = rawData_chr_Y{chr,3};
		TOTplot{chr}       = rawDataAll_chr_Y{chr};
	end;
end;


%% -----------------------------------------------------------------------------------------
% Setup for main figure generation.
%-------------------------------------------------------------------------------------------
ploidy = str2num(ploidyEstimateString);
[chr_breaks, chrCopyNum, ploidyAdjust] = FindChrSizes_4(Aneuploidy,CNVplot2,ploidy,num_chrs,chr_in_use)
fprintf('\n');

full_data_threshold = floor(bases_per_bin/100);

fig = figure(1);
set(gcf, 'Position', [0 70 1024 600]);
data_mode = 3;
for chr = 1:num_chrs
    if (chr_in_use(chr) == 1)
        if (data_mode == 1)			% Regenerate chr plot data if the save file does not exist.
            TOTplot{chr}                           = chr_SNPdata{chr,1}+chr_SNPdata{chr,2}+chr_SNPdata{chr,3};  % TOT data
            TOTave{chr}                            = sum(TOTplot{chr})/length(TOTplot{chr});
            TOTplot2{chr}                          = TOTplot{chr}/full_data_threshold;
            TOTplot2{chr}(TOTplot2{chr} > 1)       = 1;
            TOTave2{chr}                           = sum(TOTplot2{chr})/length(TOTplot2{chr});

			HOMplot{chr}                           = chr_SNPdata{chr,1};  % HOM data
            HOMave{chr}                            = sum(HOMplot{chr})/length(HOMplot{chr});
            HOMplot2{chr}                          = HOMplot{chr}/full_data_threshold;
            HOMplot2{chr}(HOMplot2{chr} > 1)       = 1;
            HOMave2{chr}                           = sum(HOMplot2{chr})/length(HOMplot2{chr});

            HETplot{chr}                           = chr_SNPdata{chr,2};  % HET data
            HETave{chr}                            = sum(HETplot{chr})/length(HETplot{chr});
            HETplot2{chr}                          = HETplot{chr}/full_data_threshold;
            HETplot2{chr}(HETplot2{chr} > 1)       = 1;
            HETave2{chr}                           = sum(HETplot2{chr})/length(HETplot2{chr});

            oddHETplot{chr}                        = chr_SNPdata{chr,3};  % oddHET data
            oddHETave{chr}                         = sum(oddHETplot{chr})/length(oddHETplot{chr});
            oddHETplot2{chr}                       = oddHETplot{chr}/full_data_threshold;
            oddHETplot2{chr}(oddHETplot2{chr} > 1) = 1;
            oddHETave2{chr}                        = sum(oddHETplot2{chr})/length(oddHETplot2{chr});
	    elseif (data_mode == 2)
			%% Details from LOH_v2a.m :
            % Regenerate chr plot data if the save file does not exist.
            TOTplot{chr}                                  = chr_SNPdata{chr,1}+chr_SNPdata{chr,2}+chr_SNPdata{chr,3};  % TOT data
			HOMplot{chr}                                  = chr_SNPdata{chr,1};  % HOM data
            HETplot{chr}                                  = chr_SNPdata{chr,2};  % HET data
            oddHETplot{chr}                               = chr_SNPdata{chr,3};  % oddHET data

            TOTave{chr}                                   = sum(TOTplot{chr})/length(TOTplot{chr});
            TOTplot2{chr}                                 = TOTplot{chr}/full_data_threshold;
            TOTplot2{chr}(TOTplot2{chr} > 1)              = 1;
            TOTave2{chr}                                  = sum(TOTplot2{chr})/length(TOTplot2{chr});

			HOMave{chr}                                   = sum(HOMplot{chr})/length(HOMplot{chr});
            HOMplot2{chr}                                 = HOMplot{chr}/full_data_threshold;
            HOMplot2{chr}(HOMplot2{chr} > 1)              = 1;
            HOMave2{chr}                                  = sum(HOMplot2{chr})/length(HOMplot2{chr});
            HETave{chr}                                   = sum(HETplot{chr})/length(HETplot{chr});
            HETplot2{chr}                                 = HETplot{chr}/full_data_threshold;
            HETplot2{chr}(HETplot2{chr} > 1)              = 1;
            HETave2{chr}                                  = sum(HETplot2{chr})/length(HETplot2{chr});
            oddHETave{chr}                                = sum(oddHETplot{chr})/length(oddHETplot{chr});
            oddHETplot2{chr}                              = oddHETplot{chr}/full_data_threshold;
            oddHETplot2{chr}(oddHETplot2{chr} > 1)        = 1;
            oddHETave2{chr}                               = sum(oddHETplot2{chr})/length(oddHETplot2{chr});
	    elseif (data_mode == 3)
			%% Details from LOH_v3a.m :
            % Regenerate chr plot data if the save file does not exist.
            TOTplot{chr}                                  = chr_SNPdata{chr,1}+chr_SNPdata{chr,2}+chr_SNPdata{chr,3};  % TOT data
            TOTave{chr}                                   = sum(TOTplot{chr})/length(TOTplot{chr});
            TOTplot2{chr}                                 = TOTplot{chr}/full_data_threshold;
            TOTplot2{chr}(TOTplot2{chr} > 1)              = 1;
            TOTave2{chr}                                  = sum(TOTplot2{chr})/length(TOTplot2{chr});

			HOMplot{chr}                                  = chr_SNPdata{chr,1};  % HOM data
            HOMave{chr}                                   = sum(HOMplot{chr})/length(HOMplot{chr});
            HOMplot2{chr}                                 = HOMplot{chr}/full_data_threshold;
            HOMplot2{chr}(HOMplot2{chr} > 1)              = 1;

            HETplot{chr}                                  = chr_SNPdata{chr,2};  % HET data
            HETave{chr}                                   = sum(HETplot{chr})/length(HETplot{chr});
            HETplot2{chr}                                 = HETplot{chr}/full_data_threshold;
            HETplot2{chr}(HETplot2{chr} > 1)              = 1;
            HETave2{chr}                                  = sum(HETplot2{chr})/length(HETplot2{chr});

            oddHETplot{chr}                               = chr_SNPdata{chr,3};  % oddHET data
            oddHETave{chr}                                = sum(oddHETplot{chr})/length(oddHETplot{chr});
            oddHETplot2{chr}                              = oddHETplot{chr}/full_data_threshold;
            oddHETplot2{chr}(oddHETplot2{chr} > 1)        = 1;
        elseif (data_mode == 4)
        end;
    end;
end;
fprintf('\n');
largestChr = find(chr_width == max(chr_width));

%% -----------------------------------------------------------------------------------------
% Setup for linear-view figure generation.
%-------------------------------------------------------------------------------------------
if (Linear_display == true)
	Linear_fig = figure(2);
	Linear_genome_size   = sum(chr_size);
	Linear_Chr_max_width = 0.85;
	Linear_left_start    = 0.07;
	Linear_left_chr_gap  = 0.01*7/(num_chrs-1);
	Linear_height        = 0.6;
	Linear_base          = 0.1;
	Linear_TickSize      = -0.01;  %negative for outside, percentage of longest chr figure.
	maxY                 = ploidyBase*2;
	Linear_left          = Linear_left_start;
end;

%% -----------------------------------------------------------------------------------------
% Make figures
%-------------------------------------------------------------------------------------------
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		figure(fig);
		% make standard chr cartoons.
		left   = chr_posX(chr);
		bottom = chr_posY(chr);
		width  = chr_width(chr);
		height = chr_height(chr);
		subplot('Position',[left bottom width height]);
		fprintf(['\tfigposition = [' num2str(left) ' | ' num2str(bottom) ' | ' num2str(width) ' | ' num2str(height) ']\n']);
		hold on;

		c_prev = colorInit;
		c_post = colorInit;
		c_     = c_prev;
		infill = zeros(1,length(HETplot2{chr}));
		colors = [];

		% standard : determines the color of each bin.
		for i = 1:length(TOTplot2{chr})+1;
			if (i-1 < length(TOTplot2{chr}))
				c_tot_post = TOTplot2{chr}(i)+TOTplot2{chr}(i);
				if (c_tot_post == 0)
					c_post = colorNoData;
				else
					%c_post =   colorHET*HETplot2{chr}(i) + ...
					%           colorHOM*HOMplot2{chr}(i) + ...
					%           colorNoData*(1-min([HETplot2{chr}(i)+HOMplot2{chr}(i) 1]));
					%colorMix = colorHET   *HETplot2   {chr}(i)/TOTplot2{chr}(i) + ...
					%           colorOddHET*oddHETplot2{chr}(i)/TOTplot2{chr}(i) + ...
					%           colorHOM   *HOMplot2   {chr}(i)/TOTplot2{chr}(i);
					colorMix = colorHET   *   HETplot2{chr}(i)/TOTplot2{chr}(i) + ...
					           colorOddHET*oddHETplot2{chr}(i)/TOTplot2{chr}(i) + ...
					           colorHOM   *   HOMplot2{chr}(i)/TOTplot2{chr}(i);
					c_post =   colorMix   *   min(1,TOTplot2{chr}(i)) + ...
					           colorNoData*(1-min(1,TOTplot2{chr}(i)));
					%colorNoData*(1-min([HETplot2{chr}(i)+oddHETplot2{chr}(i)+HOMplot2{chr}(i) 1]));
				end;
			else
				c_post = colorInit;
			end;
			colors(i,1) = c_post(1);
			colors(i,2) = c_post(2);
			colors(i,3) = c_post(3);
		end;

		% standard : draw colorbars.
		for i = 1:length(HETplot2{chr})+1;
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

		%% standard : cgh plot section.
		c_ = [0 0 0];
		fprintf(['main-plot : chr' num2str(chr) ':' num2str(length(CNVplot2{chr})) '\n']);
		for i = 1:length(CNVplot2{chr});
			x_ = [i i i-1 i-1];
			if (CNVplot2{chr}(i) == 0)
				CNVhistValue = 1;
			else
				CNVhistValue = CNVplot2{chr}(i);
			end;

			% The CNV-histogram values were normalized to a median value of 1.
			% The ratio of 'ploidy' to 'ploidyBase' determines where the data is displayed relative to the median line.
			startY = maxY/2;
			if (Low_quality_ploidy_estimate == true)
				endY = CNVhistValue*ploidy*ploidyAdjust;
			else
				endY = CNVhistValue*ploidy;
			end;
			y_ = [startY endY endY startY];

			% makes a blackbar for each bin.
			f = fill(x_,y_,c_);
			set(f,'linestyle','none');
		end;

		%% standard : draw lines across plots for easier interpretation of CNV regions.
		x2 = chr_size(chr)/bases_per_bin;
		plot([0; x2], [maxY/2; maxY/2],'color',[0 0 0]);  % 2n line.
		switch ploidyBase
			case 1
			case 2
				line([0 x2], [maxY/4*1 maxY/4*1],'Color',[0.85 0.85 0.85]);
				line([0 x2], [maxY/4*3 maxY/4*3],'Color',[0.85 0.85 0.85]);
			case 3
				line([0 x2], [maxY/6*1 maxY/6*1],'Color',[0.85 0.85 0.85]);
				line([0 x2], [maxY/6*2 maxY/6*2],'Color',[0.85 0.85 0.85]);
				line([0 x2], [maxY/6*4 maxY/6*4],'Color',[0.85 0.85 0.85]);
				line([0 x2], [maxY/6*5 maxY/6*5],'Color',[0.85 0.85 0.85]);
			case 4
				line([0 x2], [maxY/8*1 maxY/8*1],'Color',[0.85 0.85 0.85]);
				line([0 x2], [maxY/8*2 maxY/8*2],'Color',[0.85 0.85 0.85]);
				line([0 x2], [maxY/8*3 maxY/8*3],'Color',[0.85 0.85 0.85]);
				line([0 x2], [maxY/8*5 maxY/8*5],'Color',[0.85 0.85 0.85]);
				line([0 x2], [maxY/8*6 maxY/8*6],'Color',[0.85 0.85 0.85]);
				line([0 x2], [maxY/8*7 maxY/8*7],'Color',[0.85 0.85 0.85]);
		end;
		% standard : end cgh plot section.

		% standard : axes labels etc.
		hold off;
		xlim([0,chr_size(chr)/bases_per_bin]);

		%% standard : modify y axis limits to show annotation locations if any are provided.
		if (length(annotations) > 0)
			ylim([-maxY/10*1.5,maxY]);
		else
			ylim([0,maxY]);
		end;
		set(gca,'YTick',[]);
		set(gca,'TickLength',[(TickSize*chr_size(largestChr)/chr_size(chr)) 0]); %ensures same tick size on all subfigs.
		ylabel(chr_label{chr}, 'Rotation', 90, 'HorizontalAlign', 'center', 'VerticalAlign', 'bottom');
		set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
		set(gca,'XTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8','2.0','2.2','2.4','2.6','2.8','3.0','3.2'});
		set(gca,'YTick',[0 maxY/4 maxY/2 maxY/4*3 maxY]);
		set(gca,'YTickLabel',{'','','','',''});
		text(-50000/bases_per_bin, maxY/4,   '1','HorizontalAlignment','right','Fontsize',5);
		text(-50000/bases_per_bin, maxY/2,   '2','HorizontalAlignment','right','Fontsize',5);
		text(-50000/bases_per_bin, maxY/4*3, '3','HorizontalAlignment','right','Fontsize',5);
		text(-50000/bases_per_bin, maxY,     '4','HorizontalAlignment','right','Fontsize',5);
		set(gca,'FontSize',6);
		if (chr == find(chr_posY == max(chr_posY)))
			if (strcmp(projectParent,projectChild) == 1)
				title([ projectChild ' SNP map'],'Interpreter','none','FontSize',12);
			else
				title([ projectChild ' vs. (parent)' projectParent ' SNP/LOH map'],'Interpreter','none','FontSize',12);
			end;
		end;
		hold on;
		% standard : end axes labels etc.

		% standard : show centromere.
		if (chr_size(chr) < 100000)
			Centromere_format = 1;
		else
			Centromere_format = Centromere_format_default;
		end;
		% standard : show centromere outlines and horizontal marks.
		x1 = cen_start(chr)/bases_per_bin;
		x2 = cen_end(chr)/bases_per_bin;
		leftEnd  = 0.5*5000/bases_per_bin;
		rightEnd = (chr_size(chr) - 0.5*5000)/bases_per_bin;
		if (Centromere_format == 0)
			% standard chromosome cartoons in a way which will not cause segfaults when running via commandline.
			dx = cen_tel_Xindent;
			dy = cen_tel_Yindent;
			% draw white triangles at corners and centromere locations.
			fill([leftEnd   leftEnd   leftEnd+dx ],       [maxY-dy   maxY      maxY],         [1.0 1.0 1.0], 'EdgeColor',[1.0 1.0 1.0]);	% top left corner.
			fill([leftEnd   leftEnd   leftEnd+dx ],       [dy        0         0   ],         [1.0 1.0 1.0], 'EdgeColor',[1.0 1.0 1.0]);	% bottom left corner.
			fill([rightEnd  rightEnd  rightEnd-dx],       [maxY-dy   maxY      maxY],         [1.0 1.0 1.0], 'EdgeColor',[1.0 1.0 1.0]);	% top right corner.
			fill([rightEnd  rightEnd  rightEnd-dx],       [dy        0         0   ],         [1.0 1.0 1.0], 'EdgeColor',[1.0 1.0 1.0]);	% bottom right corner.
			fill([x1-dx     x1        x2           x2+dx],[maxY      maxY-dy   maxY-dy  maxY],[1.0 1.0 1.0], 'EdgeColor',[1.0 1.0 1.0]);	% top centromere.
			fill([x1-dx     x1        x2           x2+dx],[0         dy        dy       0   ],[1.0 1.0 1.0], 'EdgeColor',[1.0 1.0 1.0]);	% bottom centromere.
			% draw outlines of chromosome cartoon.   (drawn after horizontal lines to that cartoon edges are not interrupted by horiz lines.
			plot([leftEnd   leftEnd   leftEnd+dx   x1-dx   x1        x2        x2+dx    rightEnd-dx   rightEnd   rightEnd   rightEnd-dx   x2+dx   x2   x1   x1-dx   leftEnd+dx   leftEnd],...
				 [dy        maxY-dy   maxY         maxY    maxY-dy   maxY-dy   maxY     maxY          maxY-dy    dy         0             0       dy   dy   0       0            dy     ],...
				  'Color',[0 0 0]);
		elseif (Centromere_format == 1)
			leftEnd  = 0;
			rightEnd = chr_size(chr)/bases_per_bin;

			% Minimal outline for examining very small sequence regions, such as C.albicans MTL locus.
			plot([leftEnd   leftEnd   rightEnd   rightEnd   leftEnd], [0   maxY   maxY   0   0], 'Color',[0 0 0]);
		end;
		% standard : end show centromere.

		% standard : show annotation locations
		if (show_annotations) && (length(annotations) > 0)
			plot([leftEnd rightEnd], [-maxY/10*1.5 -maxY/10*1.5],'color',[0 0 0]);
			hold on;
			annotation_location = (annotation_start+annotation_end)./2;
			for i = 1:length(annotation_location)
				if (annotation_chr(i) == chr)
					annotationloc = annotation_location(i)/bases_per_bin-0.5*(5000/bases_per_bin);
					annotationStart = annotation_start(i)/bases_per_bin-0.5*(5000/bases_per_bin);
					annotationEnd   = annotation_end(i)/bases_per_bin-0.5*(5000/bases_per_bin);
					if (strcmp(annotation_type{i},'dot') == 1)
						plot(annotationloc,-maxY/10*1.5,'k:o','MarkerEdgeColor',annotation_edgecolor{i}, ...
						                                      'MarkerFaceColor',annotation_fillcolor{i}, ...
						                                      'MarkerSize',     annotation_size(i));
					elseif (strcmp(annotation_type{i},'block') == 1)
						fill([annotationStart annotationStart annotationEnd annotationEnd], ...
						     [-maxY/10*(1.5+0.75) -maxY/10*(1.5-0.75) -maxY/10*(1.5-0.75) -maxY/10*(1.5+0.75)], ...
						     annotation_fillcolor{i},'EdgeColor',annotation_edgecolor{i});
					end;
				end;
			end;
			hold off;
		end;
		% standard : end show annotation locations.

		% standard : make CGH histograms to the right of the main chr cartoons.
		if (HistPlot == true)
			width     = 0.020;
			height    = chr_height(chr);
			bottom    = chr_posY(chr);
			histAll   = [];
			histAll2  = [];
			smoothed  = [];
			smoothed2 = [];
			for segment = 1:length(chrCopyNum{chr})
				subplot('Position',[(left+chr_width(chr)+0.005)+width*(segment-1) bottom width height]);
				% The CNV-histogram values were normalized to a median value of 1.
				for i = round(1+length(CNVplot2{chr})*chr_breaks{chr}(segment)):round(length(CNVplot2{chr})*chr_breaks{chr}(segment+1))
					if (Low_quality_ploidy_estimate == true)
						histAll{segment}(i) = CNVplot2{chr}(i)*ploidy*ploidyAdjust;
					else
						histAll{segment}(i) = CNVplot2{chr}(i)*ploidy;
					end;
				end;
				% make a histogram of CGH data, then smooth it for display.
				histAll{segment}(histAll{segment}<=0)             = [];
				histAll{segment}(histAll{segment}>ploidyBase*2+2) = ploidyBase*2+2;
				histAll{segment}(length(histAll{segment})+1)      = 0;   % endpoints added to ensure histogram bounds.
				histAll{segment}(length(histAll{segment})+1)      = ploidyBase*2+2;
				smoothed{segment}    = smooth_gaussian(hist(histAll{segment},(ploidyBase*2+2)*50),5,20);
				% make a smoothed version of just the endpoints used to ensure histogram bounds.
				histAll2{segment}(1) = 0;
				histAll2{segment}(2) = ploidyBase*2+2;
				smoothed2{segment}   = smooth_gaussian(hist(histAll2{segment},(ploidyBase*2+2)*50),5,20)*4;
				% subtract the smoothed endpoints from the histogram to remove the influence of the added endpoints.
				smoothed{segment}    = (smoothed{segment}-smoothed2{segment});
				smoothed{segment}    = smoothed{segment}/max(smoothed{segment});
				hold on;
				for i = 1:(ploidyBase*2-1)
					plot([0; 1],[i*50; i*50],'color',[0.75 0.75 0.75]);
				end;
				area(smoothed{segment},1:length(smoothed{segment}),'FaceColor',[0 0 0]);
				hold off;
				set(gca,'YTick',[]);
				set(gca,'XTick',[]);
				xlim([0,1]);
				ylim([0,(ploidyBase*2)*50]);
			end;
		end;
		% standard : end of CGH histograms at right.

		% standard : places chr copy number to the right of the main chr cartoons.
	    if (ChrNum == true)
	        % subplot to show chr copy number value.
	        width  = 0.020;	
	        height = chr_height(chr);
	        bottom = chr_posY(chr);
	        if (HistPlot == true)
	            subplot('Position',[(left + chr_width(chr) + 0.005 + width*(length(chrCopyNum{chr})-1) + width+0.001) bottom width height]);
	        else
	            subplot('Position',[(left + chr_width(chr) + 0.005) bottom width height]);
	        end;
	        axis off square;
	        set(gca,'YTick',[]);
	        set(gca,'XTick',[]);
	        if (length(chrCopyNum{chr}) == 1)
	            chr_string = num2str(chrCopyNum{chr}(1));
	        else
	            chr_string = num2str(chrCopyNum{chr}(1));
	            for i = 2:length(chrCopyNum{chr})
	                chr_string = [chr_string ',' num2str(chrCopyNum{chr}(i))];
	            end;
	        end;
	        text(0.1,0.5, chr_string,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',12);
	    end;
		%% END standard draw section.

		%% Linear figure draw section.
        if (Linear_display == true)
            figure(Linear_fig);
            Linear_width = Linear_Chr_max_width*chr_size(chr)/Linear_genome_size;
            subplot('Position',[Linear_left Linear_base Linear_width Linear_height]);
            Linear_left = Linear_left + Linear_width + Linear_left_chr_gap;
            hold on;
            title(chr_label{chr},'Interpreter','none','FontSize',10);

            % linear : draw colorbars.
            for i = 1:length(HETplot2{chr})+1;
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

			%% linear : cgh plot section.
			c_ = [0 0 0];
			fprintf(['linear-plot : chr' num2str(chr) ':' num2str(length(CNVplot2{chr})) '\n']);
			for i = 1:length(CNVplot2{chr});
				x_ = [i i i-1 i-1];
				if (CNVplot2{chr}(i) == 0)
					CNVhistValue = 1;
				else
					CNVhistValue = CNVplot2{chr}(i);
				end;
				% The CNV-histogram values were normalized to a median value of 1.
				% The ratio of 'ploidy' to 'ploidyBase' determines where the data is displayed relative to the median line.
				startY = maxY/2;
				if (Low_quality_ploidy_estimate == true)
					endY = CNVhistValue*ploidy*ploidyAdjust;
				else
					endY = CNVhistValue*ploidy;
				end;
				y_ = [startY endY endY startY];
				% makes a blackbar for each bin.
				f = fill(x_,y_,c_);
				set(f,'linestyle','none');
			end;

			%% linear : draw lines across plots for easier interpretation of CNV regions.
			x2 = chr_size(chr)/bases_per_bin;
			plot([0; x2], [maxY/2; maxY/2],'color',[0 0 0]);  % 2n line.
            switch ploidyBase
                case 1
                case 2
                    line([0 x2], [maxY/4*1 maxY/4*1],'Color',[0.85 0.85 0.85]);
                    line([0 x2], [maxY/4*3 maxY/4*3],'Color',[0.85 0.85 0.85]);
                case 3
                    line([0 x2], [maxY/6*1 maxY/6*1],'Color',[0.85 0.85 0.85]);
                    line([0 x2], [maxY/6*2 maxY/6*2],'Color',[0.85 0.85 0.85]);
                    line([0 x2], [maxY/6*4 maxY/6*4],'Color',[0.85 0.85 0.85]);
                    line([0 x2], [maxY/6*5 maxY/6*5],'Color',[0.85 0.85 0.85]);
                case 4
                    line([0 x2], [maxY/8*1 maxY/8*1],'Color',[0.85 0.85 0.85]);
                    line([0 x2], [maxY/8*2 maxY/8*2],'Color',[0.85 0.85 0.85]);
                    line([0 x2], [maxY/8*3 maxY/8*3],'Color',[0.85 0.85 0.85]);
                    line([0 x2], [maxY/8*5 maxY/8*5],'Color',[0.85 0.85 0.85]);
                    line([0 x2], [maxY/8*6 maxY/8*6],'Color',[0.85 0.85 0.85]);
                    line([0 x2], [maxY/8*7 maxY/8*7],'Color',[0.85 0.85 0.85]);
            end;
            %% linear : end cgh plot section.

            % linear : show segmental anueploidy breakpoints.
            if (displayBREAKS == true)
                for segment = 2:length(chr_breaks{chr})-1
                    bP = chr_breaks{chr}(segment)*length(HETplot2{chr});
                    c_ = [0 0 1];
                    x_ = [bP bP bP-1 bP-1];
                    y_ = [0 maxY maxY 0];
                    f = fill(x_,y_,c_);
                    set(f,'linestyle','none');
                end;
            end;

			% linear : show centromere.
		    if (chr_size(chr) < 100000)
		        Centromere_format = 1;
		    else
		        Centromere_format = Centromere_format_default;
		    end;
            x1 = cen_start(chr)/bases_per_bin;
            x2 = cen_end(chr)/bases_per_bin;
            leftEnd  = 0.5*5000/bases_per_bin;
            rightEnd = (chr_size(chr) - 0.5*5000)/bases_per_bin;
            if (Centromere_format == 0)
                % standard chromosome cartoons in a way which will not cause segfaults when running via commandline.
                dx = cen_tel_Xindent; %5*5000/bases_per_bin;
                dy = cen_tel_Yindent; %maxY/10;
				% draw white triangles at corners and centromere locations.
				fill([leftEnd   leftEnd   leftEnd+dx ],       [maxY-dy   maxY      maxY],         [1.0 1.0 1.0], 'EdgeColor',[1.0 1.0 1.0]);  % top left corner.
				fill([leftEnd   leftEnd   leftEnd+dx ],       [dy        0         0   ],         [1.0 1.0 1.0], 'EdgeColor',[1.0 1.0 1.0]);  % bottom left corner.
				fill([rightEnd  rightEnd  rightEnd-dx],       [maxY-dy   maxY      maxY],         [1.0 1.0 1.0], 'EdgeColor',[1.0 1.0 1.0]);  % top right corner.
				fill([rightEnd  rightEnd  rightEnd-dx],       [dy        0         0   ],         [1.0 1.0 1.0], 'EdgeColor',[1.0 1.0 1.0]);  % bottom right corner.
				fill([x1-dx     x1        x2           x2+dx],[maxY      maxY-dy   maxY-dy  maxY],[1.0 1.0 1.0], 'EdgeColor',[1.0 1.0 1.0]);  % top centromere.
				fill([x1-dx     x1        x2           x2+dx],[0         dy        dy       0   ],[1.0 1.0 1.0], 'EdgeColor',[1.0 1.0 1.0]);  % bottom centromere.
                % draw outlines of chromosome cartoon.   (drawn after horizontal lines to that cartoon edges are not interrupted by horiz lines.
                plot([leftEnd   leftEnd   leftEnd+dx   x1-dx   x1        x2        x2+dx    rightEnd-dx   rightEnd   rightEnd   rightEnd-dx   x2+dx   x2   x1   x1-dx   leftEnd+dx   leftEnd],...
                     [dy        maxY-dy   maxY         maxY    maxY-dy   maxY-dy   maxY     maxY          maxY-dy    dy         0             0       dy   dy   0       0            dy],...
                      'Color',[0 0 0]);
			elseif (Centromere_format == 1) 
				leftEnd  = 0;
	            rightEnd = chr_size(chr)/bases_per_bin;

	            % Minimal outline for examining very small sequence regions, such as C.albicans MTL locus.
	            plot([leftEnd   leftEnd   rightEnd   rightEnd   leftEnd], [0   maxY   maxY   0   0], 'Color',[0 0 0]);
            end;
            % linear : end show centromere.

			% linear : show annotation locations
            if (show_annotations) && (length(annotations) > 0)
                plot([leftEnd rightEnd], [-maxY/10*1.5 -maxY/10*1.5],'color',[0 0 0]);
                hold on;
                annotation_location = (annotation_start+annotation_end)./2;
                for i = 1:length(annotation_location)
                    if (annotation_chr(i) == chr)
                        annotationloc = annotation_location(i)/bases_per_bin-0.5*(5000/bases_per_bin);
                        annotationStart = annotation_start(i)/bases_per_bin-0.5*(5000/bases_per_bin);
                        annotationEnd   = annotation_end(i)/bases_per_bin-0.5*(5000/bases_per_bin);
                        if (strcmp(annotation_type{i},'dot') == 1)
                            plot(annotationloc,-maxY/10*1.5,'k:o','MarkerEdgeColor',annotation_edgecolor{i}, ...
                                                                  'MarkerFaceColor',annotation_fillcolor{i}, ...
                                                                  'MarkerSize',     annotation_size(i));
                        elseif (strcmp(annotation_type{i},'block') == 1)
                            fill([annotationStart annotationStart annotationEnd annotationEnd], ...
                                 [-maxY/10*(1.5+0.75) -maxY/10*(1.5-0.75) -maxY/10*(1.5-0.75) -maxY/10*(1.5+0.75)], ...
                                 annotation_fillcolor{i},'EdgeColor',annotation_edgecolor{i});
                        end;
                    end;
                end;
                hold off;
            end;
            % linear : end show annotation locations.

            %% linear : Final formatting stuff.
            xlim([0,chr_size(chr)/bases_per_bin]);
            % modify y axis limits to show annotation locations if any are provided.
            if (length(annotations) > 0)
                ylim([-maxY/10*1.5,maxY]);
            else
                ylim([0,maxY]);
            end;
            set(gca,'TickLength',[(Linear_TickSize*chr_size(1)/chr_size(chr)) 0]); %ensures same tick size on all subfigs.
            set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
            set(gca,'XTickLabel',[]);
            if (chr == 1)
                if (strcmp(projectParent,projectChild) == 1)
                    ylabel(projectChild, 'Rotation', 0, 'HorizontalAlign', 'right', 'VerticalAlign', 'bottom','Interpreter','none','FontSize',5);
                else
                    ylabel({projectChild;'vs. (parent)';projectParent}, 'Rotation', 0, 'HorizontalAlign', 'right', 'VerticalAlign', 'bottom','Interpreter','none','FontSize',5);
                end;
                set(gca,'YTick',[0 maxY/4 maxY/2 maxY/4*3 maxY]);
                set(gca,'YTickLabel',{'','','','',''});
				text(-50000/bases_per_bin, maxY/4,   '1','HorizontalAlignment','right','Fontsize',5);
                text(-50000/bases_per_bin, maxY/2,   '2','HorizontalAlignment','right','Fontsize',5);
                text(-50000/bases_per_bin, maxY/4*3, '3','HorizontalAlignment','right','Fontsize',5);
                text(-50000/bases_per_bin, maxY,     '4','HorizontalAlignment','right','Fontsize',5);
            else
                set(gca,'YTick',[]);
                set(gca,'YTickLabel',[]);
            end;
            set(gca,'FontSize',6);
            % linear : end final reformatting.

            % shift back to main figure generation.
            figure(fig);
            hold on;
        end;
    end;
end;


if (strcmp(projectParent,projectChild) == 1)
	set(fig,'PaperPosition',[0 0 8 6]*2);
    saveas(fig, [projectChildDir 'fig.CNV-SNP-map.1.eps'], 'epsc');
	saveas(fig, [projectChildDir 'fig.CNV-SNP-map.1.png'], 'png');
	set(Linear_fig,'PaperPosition',[0 0 8 0.62222222]*2);
    saveas(Linear_fig, [projectChildDir 'fig.CNV-SNP-map.2.eps'], 'epsc');
	saveas(Linear_fig, [projectChildDir 'fig.CNV-SNP-map.2.png'], 'png');
else
	set(fig,'PaperPosition',[0 0 8 6]*2);
    saveas(fig, [projectChildDir 'fig.CNV-LOH-map.1.eps'], 'epsc');
	saveas(fig, [projectChildDir 'fig.CNV-LOH-map.1.png'], 'png');
	set(Linear_fig,'PaperPosition',[0 0 8 0.62222222]*2);
    saveas(Linear_fig, [projectChildDir 'fig.CNV-LOH-map.2.eps'], 'epsc');
	saveas(Linear_fig, [projectChildDir 'fig.CNV-LOH-map.2.png'], 'png');
end;
delete(fig);
delete(Linear_fig);
end
