function [] = CNV_SNP_hapmap_v1(main_dir,user,genomeUser,project,hapmap,genome,ploidyEstimateString,ploidyBaseString, ...
                                SNP_verString,LOH_verString,CNV_verString,displayBREAKS);
%% ========================================================================
%    Centromere_format          : Controls how centromeres are depicted.   [0..2]   '2' is pinched cartoon default.
%    bases_per_bin              : Controls bin sizes for SNP/CGH fractions of plot.
%    scale_type                 : 'Ratio' or 'Log2Ratio' y-axis scaling of copy number.
%                                 'Log2Ratio' does not properly scale CGH data by ploidy.
%    Chr_max_width              : max width of chrs as fraction of figure width.
Centromere_format           = 0;
Chr_max_width               = 0.8;
colorBars                   = true;
blendColorBars              = false;
show_annotations            = true;
Yscale_nearest_even_ploidy  = true;
HistPlot                    = true;
ChrNum                      = true;
Linear_display              = true;
Low_quality_ploidy_estimate = true

%% =========================================================================================
% Load FASTA file name from 'reference.txt' file for project.
%-------------------------------------------------------------------------------------------
userReference    = [main_dir 'users/' user '/genomes/' genome '/reference.txt'];
defaultReference = [main_dir 'users/default/genomes/' genome '/reference.txt'];
if (exist(userReference,'file') == 0)   
	FASTA_string = strtrim(fileread(defaultReference));
else                    
	FASTA_string = strtrim(fileread(userReference));
end;
[FastaPath,FastaName,FastaExt] = fileparts(FASTA_string);


%% =========================================================================================
% Control variables for Candida albicans SC5314.
%-------------------------------------------------------------------------------------------
projectDir  = [main_dir 'users/' user '/projects/' project '/'];

if (exist([[main_dir 'users/default/hapmaps/' hapmap '/']],'dir') == 7)
	hapmapDir = [main_dir 'users/default/hapmaps/' hapmap '/'];
else
	hapmapDir = [main_dir 'users/' user '/hapmaps/' hapmap '/'];
end;

genomeDir  = [main_dir 'users/' genomeUser '/genomes/' genome '/'];

[centromeres, chr_sizes, figure_details, annotations, ploidy_default] = Load_genome_information_1(projectDir,genomeDir, genome);
[Aneuploidy]                                                          = Load_dataset_information_1(projectDir, project);

num_chrs = length(chr_sizes);

for i = 1:length(chr_sizes)
	chr_size(i)  = 0;
	cen_start(i) = 0;
	cen_end(i)   = 0;
end;
for i = 1:length(chr_sizes)
	chr_size(chr_sizes(i).chr)    = chr_sizes(i).size;
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
		if (strcmp(figure_details(i).label,'Key') == 1)
			key_posX   = figure_details(i).posX;
			key_posY   = figure_details(i).posY;
			key_width  = figure_details(i).width;
			key_height = figure_details(i).height;
		end;
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

%% This block is normally calculated in FindChrSizes_2 in CNV analysis.
for usedChr = 1:num_chrs
	if (chr_in_use(usedChr) == 1)
		% determine where the endpoints of ploidy segments are.
		chr_breaks{usedChr}(1) = 0.0;
		break_count = 1;
		if (length(Aneuploidy) > 0)
			for i = 1:length(Aneuploidy)
				if (Aneuploidy(i).chr == usedChr)
					break_count = break_count+1;
					chr_broken = true;
					chr_breaks{usedChr}(break_count) = Aneuploidy(i).break;
				end;
			end;
		end;
		chr_breaks{usedChr}(length(chr_breaks{usedChr})+1) = 1;
	end;
end;


%% =========================================================================================
%% =========================================================================================
%% =========================================================================================
%% = No further control variables below. ===================================================
%% =========================================================================================
%% =========================================================================================
%% =========================================================================================


% Process input ploidy.
ploidy = str2num(ploidyEstimateString);

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

fprintf(['\nGenerating LOH-map figure from ''' project ''' vs. (hapmap)''' hapmap ''' data.\n']);


%% =========================================================================================
% Initializes vectors used to hold copy number and SNP data.
%-------------------------------------------------------------------------------------------
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		% 1 category tracked : average read counts per bin.
		chr_CNVdata{chr} = zeros(1,ceil(chr_size(chr)/bases_per_bin));
		% fprintf(['0|' num2str(chr) ':' num2str(length(chr_CNVdata{chr})) '\n']);
	end;
end;
for chr = 1:length(chr_sizes)
	for j = 1:2
		% 3 SNP interpretation catagories tracked.
		%	1 : phased ratio data.
		%	2 : unphased ratio data.
		chr_length = ceil(chr_size(chr)/bases_per_bin);
		chr_SNPdata{chr,j} = cell(1,chr_length);
	end;
	% fprintf(['0|' num2str(chr) ':' num2str(length(chr_SNPdata{chr,1})) '\n']);
end;


%% =========================================================================================
% Load GC-bias corrected CGH data.
%-------------------------------------------------------------------------------------------
load([projectDir 'Common_CNV.mat']);       % 'CNVplot2','genome_CNV'
[chr_breaks, chrCopyNum, ploidyAdjust] = FindChrSizes_3(Aneuploidy,CNVplot2,ploidy,num_chrs,chr_in_use);
largestChr = find(chr_width == max(chr_width));


%% =========================================================================================
% Load CNV edges from ChARM algorithm.
%-------------------------------------------------------------------------------------------
[segmental_aneuploidy] = Load_dataset_information_1(projectDir,project);
segmental_aneuploidy = [];


%% =========================================================================================
% Load SNP/LOH data.
%-------------------------------------------------------------------------------------------
load([projectDir 'SNP_' SNP_verString '.mat']);


%% =========================================================================================
% Setup for figure generation.
%-------------------------------------------------------------------------------------------
fig = figure(1);
set(gcf, 'Position', [0 70 1024 600]);

% calculate SNP bin values.
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		for chr_bin = 1:length(chr_SNPdata{chr,1})
			SNPplot{chr,1}{chr_bin} = chr_SNPdata{chr,1}{chr_bin};
			SNPplot{chr,2}{chr_bin} = chr_SNPdata{chr,2}{chr_bin};
		end;
	end;
end;

SNPdata_all = [];
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		for chr_bin = 1:length(SNPplot{chr,1})
			TOTplot{chr}{chr_bin} = length(chr_SNPdata{chr,1}{chr_bin}) + length(chr_SNPdata{chr,2}{chr_bin});   % TOT = phased+nonphased;
			SNPdata_all           = [SNPdata_all TOTplot{chr}{chr_bin}];
		end;
	end;
end;
medianRawY = median(SNPdata_all);

% Gather median-normalized SNP data for LOWESS fitting.
SNPdata_all = [];
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		for chr_bin = 1:length(SNPplot{chr,1})
			TOTplot{chr}{chr_bin} = length(chr_SNPdata{chr,1}{chr_bin}/medianRawY) + length(chr_SNPdata{chr,2}{chr_bin}/medianRawY);   % TOT = phased+nonphased;
			SNPdata_all           = [SNPdata_all TOTplot{chr}{chr_bin}];
		end;
	end;
end;


%% =========================================================================================
% Apply GC bias correction to SNP data.
%   Number of putative SNPs vs. GCbias per standard bin.
%-------------------------------------------------------------------------------------------
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

% Gather SNP and GCratio data for LOWESS fitting.
GCratioData_all = [];
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		GCratioData_all = [GCratioData_all chr_GCratioData{chr}];
	end;
end;
medianRawY = median(SNPdata_all);
fprintf(['medianRawY = ' num2str(medianRawY) '\n']);
fprintf(['SNPdata_all     => ' num2str(length(SNPdata_all)    ) '\n']);
fprintf(['GCratioData_all => ' num2str(length(GCratioData_all)) '\n']);


%% Clean up data by:
%%    deleting GC ratio data near zero.
%%    deleting CGH data beyond 3* the median value.  (rDNA, etc.)
SNPdata_clean                                            = SNPdata_all;
GCratioData_clean                                        = GCratioData_all;
SNPdata_clean(     GCratioData_clean < 0.01            ) = [];
GCratioData_clean( GCratioData_clean < 0.01            ) = [];
GCratioData_clean( SNPdata_clean > max(medianRawY*3,3) ) = [];
SNPdata_clean(     SNPdata_clean > max(medianRawY*3,3) ) = [];
GCratioData_clean( SNPdata_clean == 0                  ) = [];
SNPdata_clean(     SNPdata_clean == 0                  ) = [];

% Perform LOWESS fitting.
rawData_X1     = GCratioData_clean;
rawData_Y1     = SNPdata_clean;
fprintf(['Lowess X:Y size : [' num2str(size(rawData_X1,1)) ',' num2str(size(rawData_X1,2)) ']:[' num2str(size(rawData_Y1,1)) ',' num2str(size(rawData_Y1,2)) ']\n']);
[fitX1, fitY1] = optimize_mylowess_SNP(rawData_X1,rawData_Y1);

% Correct data using normalization to LOWESS fitting
Y_target = 1;
% Initialize data structures, to simplify debugging.   These nested loops can be deleted without problems arising.
for chr = 1:num_chrs
	for j = 1:2
		rawData_chr_Y{chr,j}        = [];
		normalizedData_chr_Y{chr,j} = [];
		rawData_chr_X{chr}          = [];
		rawDataAll_chr_Y{chr}       = [];
		fitDataAll_chr_Y{chr}       = [];
		fitData_chr_Y{chr}          = [];
	end;
end;
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		rawData_chr_X{chr}           = chr_GCratioData{chr};
		rawDataAll_chr_Y{chr}        = cell2mat(TOTplot{chr});
		fitDataAll_chr_Y{chr}        = interp1(fitX1,fitY1,rawData_chr_X{chr},'spline');
		normalizedDataAll_chr_Y{chr} = rawDataAll_chr_Y{chr}./fitDataAll_chr_Y{chr}*Y_target;
		fitData_chr_Y{chr}           = interp1(fitX1,fitY1,rawData_chr_X{chr},'spline');

 		for chr_bin = 1:length(SNPplot{chr,1})
			rawData_chr_Y{chr,1}{chr_bin}        = length(chr_SNPdata{chr,1}{chr_bin});
			rawData_chr_Y{chr,2}{chr_bin}        = length(chr_SNPdata{chr,2}{chr_bin});

			normalizedData_chr_Y{chr,1}{chr_bin} = length(rawData_chr_Y{chr,1}{chr_bin})./fitData_chr_Y{chr}*Y_target;
			normalizedData_chr_Y{chr,2}{chr_bin} = length(rawData_chr_Y{chr,2}{chr_bin})./fitData_chr_Y{chr}*Y_target;
		end;
	end;
end;

% Move LOWESS-normalizd SNP data back into display pipeline.
SNPdata_all    = [];
SNPdata_all_1  = [];
SNPdata_all_2  = [];
cSNPdata_all   = [];
cSNPdata_all_1 = [];
cSNPdata_all_2 = [];
GCdata_all     = [];
GCdata_all_1   = [];
GCdata_all_2   = [];
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		TOTplot{chr}        = rawDataAll_chr_Y{chr};
		cTOTplot{chr}       = normalizedDataAll_chr_Y{chr};
		SNPdata_all         = [SNPdata_all    TOTplot{chr}              ];
		cSNPdata_all        = [cSNPdata_all   cTOTplot{chr}             ];
		GCdata_all          = [GCdata_all     rawData_chr_X{chr}        ];

		for chr_bin = 1:length(SNPplot{chr,1})
% darren : attempt at data normalization is failing, introduces major bug due to data compaction.
%			chr_SNPdata{chr,1}{chr_bin}  = rawData_chr_Y{chr,1}{chr_bin};
%			chr_SNPdata{chr,2}{chr_bin}  = rawData_chr_Y{chr,2}{chr_bin};
			cchr_SNPdata{chr,1}{chr_bin} = normalizedData_chr_Y{chr,1}{chr_bin};
			cchr_SNPdata{chr,2}{chr_bin} = normalizedData_chr_Y{chr,2}{chr_bin};
			SNPdata_all_1                = [SNPdata_all_1  length(chr_SNPdata{chr,1}{chr_bin})  ];
			SNPdata_all_2                = [SNPdata_all_2  length(chr_SNPdata{chr,2}{chr_bin})  ];
			cSNPdata_all_1               = [cSNPdata_all_1 length(cchr_SNPdata{chr,1}{chr_bin}) ];
			cSNPdata_all_2               = [cSNPdata_all_2 length(cchr_SNPdata{chr,2}{chr_bin}) ];
		end;
	end;
end;


%% Generate figure showing subplots of LOWESS fittings.
GCfig = figure(3);
subplot(2,3,1);
    plot(GCratioData_all,SNPdata_all,'k.','markersize',1);
    hold on;	plot(fitX1,fitY1,'r','LineWidth',2);   hold off;
    xlabel('GC ratio');   ylabel('SNP data');
    xlim([0.0 1.0]);      ylim([0 max(medianRawY*5,5)]);   axis square;
subplot(2,3,2);
	plot(GCdata_all,SNPdata_all_1,'r.','markersize',1);
	hold on;    plot(fitX1,fitY1,'k','LineWidth',2);   hold off;
	xlabel('GC ratio');   ylabel('SNP data');
	xlim([0.0 1.0]);      ylim([0 max(medianRawY*5,5)]);   axis square;
subplot(2,3,3);
	plot(GCdata_all,SNPdata_all_2,'g.','markersize',1);
	hold on;    plot(fitX1,fitY1,'k','LineWidth',2);   hold off;
	xlabel('GC ratio');   ylabel('SNP data');
	xlim([0.0 1.0]);      ylim([0 max(medianRawY*5,5)]);   axis square;

subplot(2,3,4);
	plot(GCratioData_all,cSNPdata_all,'k.','markersize',1);
	hold on;   plot([min(GCratioData_all) max(GCratioData_all)],[Y_target Y_target],'r','LineWidth',2);   hold off;
	xlabel('GC ratio');   ylabel('corrected SNP data');
	xlim([0.0 1.0]);      ylim([0 5]);                    axis square;
subplot(2,3,5);
	plot(GCdata_all,cSNPdata_all_1,'r.','markersize',1);
	hold on;   plot([min(GCratioData_all) max(GCratioData_all)],[Y_target Y_target],'k','LineWidth',2);   hold off;
	xlabel('GC ratio');   ylabel('corrected SNP data');
	xlim([0.0 1.0]);      ylim([0 5]);                    axis square;
subplot(2,3,6);
	plot(GCdata_all,cSNPdata_all_1,'g.','markersize',1);
	hold on;   plot([min(GCratioData_all) max(GCratioData_all)],[Y_target Y_target],'k','LineWidth',2);   hold off;
	xlabel('GC ratio');   ylabel('corrected SNP data');
	xlim([0.0 1.0]);      ylim([0 5]);                    axis square;

saveas(GCfig, [projectDir '/fig.GCratio_vs_SNP.png'], 'png');


%% =========================================================================================
% Setup for main figure generation.
%-------------------------------------------------------------------------------------------
% threshold for full color saturation in SNP/LOH figure.
% synced to bases_per_bin as below, or defaulted to 50.
full_data_threshold = floor(bases_per_bin/100);
fig = figure(1);
set(gcf, 'Position', [0 70 1024 600]);
data_mode = 1;
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		if (data_mode == 1)
			for chr_bin = 1:length(chr_SNPdata{chr,1})
				% Regenerate chr plot data if the save file does not exist.
				SNPs_count{chr}(chr_bin)                                     = length(chr_SNPdata{chr,1}{chr_bin}) + length(chr_SNPdata{chr,2}{chr_bin});   % the number of data points in this bin.
				SNPs_to_fullData_ratio{chr}(chr_bin)                         = SNPs_count{chr}(chr_bin)/full_data_threshold;                                % divide by the threshold for full color saturation in SNP/LOH figure.
				SNPs_to_fullData_ratio{chr}(SNPs_to_fullData_ratio{chr} > 1) = 1;                                                                           % any bins with more data than the threshold for full color saturation around limited to full saturation.

				% darren
				%if (length(chr_SNPdata{chr,1}{chr_bin}) > 0)
				%	if (mean(chr_SNPdata{chr,1}{chr_bin}) < 0.5)
				%		fprintf(['chr' num2str(chr) ':chr_bin' num2str(chr_bin) ' => A\n']);
				%	elseif (mean(chr_SNPdata{chr,1}{chr_bin}) > 0.5)
				%		fprintf(['chr' num2str(chr) ':chr_bin' num2str(chr_bin) ' => B\n']);
				%	else
				%		fprintf(['chr' num2str(chr) ':chr_bin' num2str(chr_bin) ' => ab\n']);
				%	end;
				%	% fprintf(['chr' num2str(chr) ':chr_bin' num2str(chr_bin) ' => dataAve = (' num2str(mean(chr_SNPdata{chr,1}{chr_bin})) ')\n']);
				%end;
    
				HOMplot{chr}(chr_bin)                  = length(chr_SNPdata{chr,1}{chr_bin});         % phased data.
				HOMplot2{chr}(chr_bin)                 = HOMplot{chr}(chr_bin)/full_data_threshold;   %
				HOMplot2{chr}(HOMplot2{chr} > 1)       = 1;                                           %

				HETplot{chr}(chr_bin)                  = length(chr_SNPdata{chr,2}{chr_bin});         % unphased data.
				HETplot2{chr}(chr_bin)                 = HETplot{chr}(chr_bin)/full_data_threshold;   %
				HETplot2{chr}(HETplot2{chr} > 1)       = 1;                                           %
			end;
		end;
	end;
end;
fprintf('\n');
largestChr = find(chr_width == max(chr_width));


%% =========================================================================================
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


%% =========================================================================================
% Define colors for figure generation.
%-------------------------------------------------------------------------------------------
%define colors for colorBars plot
colorNoData     = [1.0   1.0   1.0  ]; %used when no data is available for the bin.
colorInit       = [0.5   0.5   0.5  ]; %external; used in blending at ends of chr.
colorHET        = [0.0   0.0   0.0  ]; % near 1:1 ratio SNPs
colorHOM        = [1.0   0.0   0.0  ]; % Hom SNPs;

%% Figure out the user of the in-use hapmap.
userHapmap    = [main_dir 'users/' user '/hapmaps/' hapmap '/'];
defaultHapmap = [main_dir 'users/default/hapmaps/' hapmap '/'];
if (exist(userHapmap,'dir') == 0)
    hapmapUser = 'default';
else
	hapmapUser = user;
end;

%% Load color names defined for hapmap;
colors_fid = fopen([main_dir 'users/' hapmapUser '/hapmaps/' hapmap '/colors.txt'], 'r');

% The swapped colors are to correct for a polarity mistake in the python preprocessing steps.
%    correcting the error there would require reprocessing all current datasets.
colorB_string = fgetl(colors_fid);
colorA_string = fgetl(colors_fid);
fclose(colors_fid);
fprintf(['\nHapmap colors:\n\tcolorA = ' colorA_string '\n\tcolorB = ' colorB_string '\n\n']);

switch colorA_string
	case 'green'
		homolog_a_color = [0.0 1.0 0.0];
	case 'blue'
		homolog_a_color = [0.0 0.0 1.0];
	case 'cyan'
		homolog_a_color = [0.0 1.0 1.0];
	case 'magenta'
		homolog_a_color = [1.0 0.0 1.0];
	case 'yellow'
		homolog_a_color = [1.0 1.0 0.0];
	otherwise
		homolog_a_color = [1/3 1/3 1/3];
end;
switch colorB_string
	case 'green'
		homolog_b_color = [0.0 1.0 0.0];
	case 'blue'
		homolog_b_color = [0.0 0.0 1.0];
	case 'cyan'
		homolog_b_color = [0.0 1.0 1.0];
	case 'magenta'
		homolog_b_color = [1.0 0.0 1.0];
	case 'yellow'
		homolog_b_color = [1.0 1.0 0.0];
	otherwise
		homolog_b_color = [1/3 1/3 1/3];
end;

% phased data colors.
	het_color       = [0.66667 0.66667 0.66667]; % heterozygous.
	% haploid colors.
	colorA          = homolog_a_color;
	colorB          = homolog_b_color;
	% diploid colors.
	colorAA         = homolog_a_color;
	colorAB         = het_color;
	colorBB         = homolog_b_color;
	% triploid colors.
	colorAAA        = homolog_a_color;
	colorAAB        = homolog_a_color*2/3 + homolog_b_color*1/3;
	colorABB        = homolog_a_color*1/3 + homolog_b_color*2/3;
	colorBBB        = homolog_b_color;
	% tetraploid colors.
	colorAAAA       = homolog_a_color;
	colorAAAB       = homolog_a_color*3/4 + homolog_b_color*1/4;
	colorAABB       = het_color;
	colorABBB       = homolog_a_color*1/4 + homolog_b_color*3/4;
	colorBBBB       = homolog_b_color;
	% pentaploid colors.
	colorAAAAA      = homolog_a_color;
	colorAAAAB      = homolog_a_color*4/5 + homolog_b_color*1/5;
	colorAAABB      = homolog_a_color*3/5 + homolog_b_color*2/5;
	colorAABBB      = homolog_a_color*2/5 + homolog_b_color*3/5;
	colorABBBB      = homolog_a_color*1/5 + homolog_b_color*4/5;
	colorBBBBB      = homolog_b_color;
	% hexaploid colors.
	colorAAAAAA     = homolog_a_color;
	colorAAAAAB     = homolog_a_color*5/6 + homolog_b_color*1/6;
	colorAAAABB     = homolog_a_color*4/6 + homolog_b_color*2/6;
	colorAAABBB     = het_color;
	colorAABBBB     = homolog_a_color*2/6 + homolog_b_color*4/6;
	colorABBBBB     = homolog_a_color*1/6 + homolog_b_color*5/6;
	colorBBBBBB     = homolog_b_color;
	% heptaploid colors.
	colorAAAAAAA    = homolog_a_color;
	colorAAAAAAB    = homolog_a_color*6/7 + homolog_b_color*1/7;
	colorAAAAABB    = homolog_a_color*5/7 + homolog_b_color*2/7;
	colorAAAABBB    = homolog_a_color*4/7 + homolog_b_color*3/7;
	colorAAABBBB    = homolog_a_color*3/7 + homolog_b_color*4/7;
	colorAABBBBB    = homolog_a_color*2/7 + homolog_b_color*5/7;
	colorABBBBBB    = homolog_a_color*1/7 + homolog_b_color*6/7;
	colorBBBBBBB    = homolog_b_color;
	% octaploid colors.
	colorAAAAAAAA   = homolog_a_color;
	colorAAAAAAAB   = homolog_a_color*7/8 + homolog_b_color*1/8;
	colorAAAAAABB   = homolog_a_color*6/8 + homolog_b_color*2/8;
	colorAAAAABBB   = homolog_a_color*5/8 + homolog_b_color*3/8;
	colorAAAABBBB   = het_color;
	colorAAABBBBB   = homolog_a_color*3/8 + homolog_b_color*5/8;
	colorAABBBBBB   = homolog_a_color*2/8 + homolog_b_color*6/8;
	colorABBBBBBB   = homolog_a_color*1/8 + homolog_b_color*7/8;
	colorBBBBBBBB   = homolog_b_color;
	% nonaploid colors.
	colorAAAAAAAAA  = homolog_a_color;
	colorAAAAAAAAB  = homolog_a_color*8/9 + homolog_b_color*1/9;
	colorAAAAAAABB  = homolog_a_color*7/9 + homolog_b_color*2/9;
	colorAAAAAABBB  = homolog_a_color*6/9 + homolog_b_color*3/9;
	colorAAAAABBBB  = homolog_a_color*5/9 + homolog_b_color*4/9;
	colorAAAABBBBB  = homolog_a_color*4/9 + homolog_b_color*5/9;
	colorAAABBBBBB  = homolog_a_color*3/9 + homolog_b_color*6/9;
	colorAABBBBBBB  = homolog_a_color*2/9 + homolog_b_color*7/9;
	colorABBBBBBBB  = homolog_a_color*1/9 + homolog_b_color*8/9;
	colorBBBBBBBBB  = homolog_b_color;

% unphased colors.
	hom_unphased_color  = [1.0     0.0     0.0    ]; % homozygous, unphased.
	het_unphased_color  = [0.66667 0.66667 0.66667]; % heterozygous.
	% haploid colors.
	unphased_color_1of1 = hom_unphased_color;
	% diploid colors.
	unphased_color_2of2 = hom_unphased_color;
	unphased_color_1of2 = het_unphased_color;
	% triploid colors.
	unphased_color_3of3 = hom_unphased_color;
	unphased_color_2of3 = hom_unphased_color*2/3 + het_unphased_color*1/3;
	% tetraploid colors.
	unphased_color_4of4 = hom_unphased_color;
	unphased_color_3of4 = hom_unphased_color*3/4 + het_unphased_color*1/4;
	unphased_color_2of4 = het_unphased_color;
	% pentaploid colors.
	unphased_color_5of5 = hom_unphased_color;
	unphased_color_4of5 = hom_unphased_color*4/5 + het_unphased_color*1/5;
	unphased_color_3of5 = hom_unphased_color*3/5 + het_unphased_color*2/5;
	% hexaploid colors.
	unphased_color_6of6 = hom_unphased_color;
	unphased_color_5of6 = hom_unphased_color*5/6 + het_unphased_color*1/6;
	unphased_color_4of6 = hom_unphased_color*4/6 + het_unphased_color*2/6;
	unphased_color_3of6 = het_unphased_color;
	% heptaploid colors.
	unphased_color_7of7 = hom_unphased_color;
	unphased_color_6of7 = hom_unphased_color*6/7 + het_unphased_color*1/7;
	unphased_color_5of7 = hom_unphased_color*5/7 + het_unphased_color*2/7;
	unphased_color_4of7 = hom_unphased_color*4/7 + het_unphased_color*3/7;
	% octaploid colors.
	unphased_color_8of8 = hom_unphased_color;
	unphased_color_7of8 = hom_unphased_color*7/8 + het_unphased_color*1/8;
	unphased_color_6of8 = hom_unphased_color*6/8 + het_unphased_color*2/8;
	unphased_color_5of8 = hom_unphased_color*5/8 + het_unphased_color*3/8;
	unphased_color_4of8 = het_unphased_color;
	% nonaploid colors.
	unphased_color_9of9 = hom_unphased_color;
	unphased_color_8of9 = hom_unphased_color*8/9 + het_unphased_color*1/9;
	unphased_color_7of9 = hom_unphased_color*7/9 + het_unphased_color*2/9;
	unphased_color_6of9 = hom_unphased_color*6/9 + het_unphased_color*3/9;
	unphased_color_5of9 = hom_unphased_color*5/9 + het_unphased_color*4/9;

%% =========================================================================================
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
	    for chr_bin = 1:length(SNPs_to_fullData_ratio{chr})+1;
	        if (chr_bin-1 < length(SNPs_to_fullData_ratio{chr}))
	            c_tot_post = SNPs_to_fullData_ratio{chr}(chr_bin)+SNPs_to_fullData_ratio{chr}(chr_bin);
	            if (c_tot_post == 0)
	                c_post = colorNoData;
	            else
					% darren
					% Define colorMix using localized copy number estimate to define SNP cutoff thresholds,
					%     then the ratio of SNP data in each SNP ratio bin.
					% For testing, consider all loci haploid, so only two ratio bins.
					ratioData_phased   = chr_SNPdata{chr,1}{chr_bin};
					ratioData_unphased = chr_SNPdata{chr,2}{chr_bin};
					localCopyEstimate  = round(CNVplot2{chr}(chr_bin));
					if (localCopyEstimate <= 0)
					elseif (localCopyEstimate == 1)
						binCounts_phased   = zeros(1,2);	
						binCounts_unphased = 0;
						if (length(ratioData_phased) > 0)
							for i = 1:length(ratioData_phased)
								    if (ratioData_phased(i) < 1/2);   binCounts_phased(1) = binCounts_phased(1)+1;
								elseif (ratioData_phased(i) > 1/2);   binCounts_phased(2) = binCounts_phased(2)+1;
								else
									binCounts_phased(1) = binCounts_phased(1)+0.5;
									binCounts_phased(2) = binCounts_phased(2)+0.5;
								end;
							end;
						end;
						if (length(ratioData_unphased) > 0)
							for i = 1:length(ratioData_unphased)
								binCounts_unphased = binCounts_unphased+1;
							end;
						end;
						colorMix = colorA              * binCounts_phased(1)/(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorB              * binCounts_phased(2)/(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           unphased_color_1of1 * binCounts_unphased /(sum(binCounts_phased)+sum(binCounts_unphased));
					elseif (localCopyEstimate == 2)
						binCounts_phased   = zeros(1,3);
						binCounts_unphased = zeros(1,2);
						if (length(ratioData_phased) > 0)
							for i = 1:length(ratioData_phased)
								    if (ratioData_phased(i) < 1/4);   binCounts_phased(1) = binCounts_phased(1)+1;
								elseif (ratioData_phased(i) > 3/4);   binCounts_phased(3) = binCounts_phased(3)+1;
								else                                  binCounts_phased(2) = binCounts_phased(2)+1;
								end;
							end;
						end;
						if (length(ratioData_unphased) > 0)
							for i = 1:length(ratioData_unphased)
								    if (ratioData_unphased(i) < 1/4);   binCounts_unphased(1) = binCounts_unphased(1)+1;
								elseif (ratioData_unphased(i) > 3/4);   binCounts_unphased(1) = binCounts_unphased(1)+1;
								else                                    binCounts_unphased(2) = binCounts_unphased(2)+1;
								end;
							end;
						end;
						colorMix = colorAA             * binCounts_phased(1)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorAB             * binCounts_phased(2)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorBB             * binCounts_phased(3)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           unphased_color_2of2 * binCounts_unphased(1)/(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           unphased_color_1of2 * binCounts_unphased(2)/(sum(binCounts_phased)+sum(binCounts_unphased));
					elseif (localCopyEstimate == 3)
						binCounts_phased   = zeros(1,4);
						binCounts_unphased = zeros(1,2);
						if (length(ratioData_phased) > 0)
							for i = 1:length(ratioData_phased)
								    if (ratioData_phased(i) < 1/6);  binCounts_phased(1) = binCounts_phased(1)+1;
								elseif (ratioData_phased(i) < 1/2);  binCounts_phased(2) = binCounts_phased(2)+1;
								elseif (ratioData_phased(i) > 5/6);  binCounts_phased(4) = binCounts_phased(4)+1;
								elseif (ratioData_phased(i) > 1/2);  binCounts_phased(3) = binCounts_phased(3)+1;
								else
									binCounts_phased(2) = binCounts_phased(2)+0.5;
									binCounts_phased(3) = binCounts_phased(3)+0.5;
								end;
							end;
						end;
						if (length(ratioData_unphased) > 0)
							for i = 1:length(ratioData_unphased)
								    if (ratioData_unphased(i) < 1/6);  binCounts_unphased(1) = binCounts_unphased(1)+1;
								elseif (ratioData_unphased(i) > 5/6);  binCounts_unphased(1) = binCounts_unphased(1)+1;
								else                                   binCounts_unphased(2) = binCounts_unphased(2)+1;
								end;
							end;
						end;
						colorMix = colorAAA            * binCounts_phased(1)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorAAB            * binCounts_phased(2)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorABB            * binCounts_phased(3)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorBBB            * binCounts_phased(4)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           unphased_color_3of3 * binCounts_unphased(1)/(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           unphased_color_2of3 * binCounts_unphased(2)/(sum(binCounts_phased)+sum(binCounts_unphased));
					elseif (localCopyEstimate == 4)
						binCounts_phased   = zeros(1,5);
						binCounts_unphased = zeros(1,3);
						if (length(ratioData_phased) > 0)
							for i = 1:length(ratioData_phased)
								    if (ratioData_phased(i) < 1/8);  binCounts_phased(1) = binCounts_phased(1)+1;
								elseif (ratioData_phased(i) < 3/8);  binCounts_phased(2) = binCounts_phased(2)+1;
								elseif (ratioData_phased(i) > 7/8);  binCounts_phased(5) = binCounts_phased(5)+1;
								elseif (ratioData_phased(i) > 5/8);  binCounts_phased(4) = binCounts_phased(4)+1;
								else                                 binCounts_phased(3) = binCounts_phased(3)+1;
								end;
							end;
						end;
						if (length(ratioData_unphased) > 0)
							for i = 1:length(ratioData_unphased)
								    if (ratioData_unphased(i) < 1/8);  binCounts_unphased(1) = binCounts_unphased(1)+1;
								elseif (ratioData_unphased(i) < 3/8);  binCounts_unphased(2) = binCounts_unphased(2)+1;
								elseif (ratioData_unphased(i) > 7/8);  binCounts_unphased(1) = binCounts_unphased(1)+1;
								elseif (ratioData_unphased(i) > 8/8);  binCounts_unphased(2) = binCounts_unphased(2)+1;
								else                                   binCounts_unphased(3) = binCounts_unphased(3)+1;
								end;
							end;
						end;
						colorMix = colorAAAA           * binCounts_phased(1)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorAAAB           * binCounts_phased(2)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorAABB           * binCounts_phased(3)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorABBB           * binCounts_phased(4)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorBBBB           * binCounts_phased(5)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           unphased_color_4of4 * binCounts_unphased(1)/(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           unphased_color_3of4 * binCounts_unphased(2)/(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           unphased_color_2of4 * binCounts_unphased(3)/(sum(binCounts_phased)+sum(binCounts_unphased));
					elseif (localCopyEstimate == 5)
						binCounts_phased   = zeros(1,6);
						binCounts_unphased = zeros(1,3);
						if (length(ratioData_phased) > 0)
							for i = 1:length(ratioData_phased)
								    if (ratioData_phased(i) < 1/10);  binCounts_phased(1) = binCounts_phased(1)+1;
								elseif (ratioData_phased(i) < 3/10);  binCounts_phased(2) = binCounts_phased(2)+1;
								elseif (ratioData_phased(i) < 1/2);   binCounts_phased(3) = binCounts_phased(3)+1;
								elseif (ratioData_phased(i) > 9/10);  binCounts_phased(6) = binCounts_phased(6)+1;
								elseif (ratioData_phased(i) > 7/10);  binCounts_phased(5) = binCounts_phased(5)+1;
								elseif (ratioData_phased(i) > 1/2);   binCounts_phased(4) = binCounts_phased(4)+1;
								else
									binCounts_phased(3) = binCounts_phased(3)+0.5;
									binCounts_phased(4) = binCounts_phased(4)+0.5;
								end;
							end;
						end;
						if (length(ratioData_unphased) > 0)
							for i = 1:length(ratioData_unphased)
								    if (ratioData_unphased(i) < 1/10);  binCounts_unphased(1) = binCounts_unphased(1)+1;
								elseif (ratioData_unphased(i) < 3/10);  binCounts_unphased(2) = binCounts_unphased(2)+1;
								elseif (ratioData_unphased(i) > 9/10);  binCounts_unphased(1) = binCounts_unphased(1)+1;
								elseif (ratioData_unphased(i) > 7/10);  binCounts_unphased(2) = binCounts_unphased(2)+1;
								else                                    binCounts_unphased(3) = binCounts_unphased(3)+1;
								end;
							end;
						end;
						colorMix = colorAAAAA          * binCounts_phased(1)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorAAAAB          * binCounts_phased(2)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorAAABB          * binCounts_phased(3)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorAABBB          * binCounts_phased(4)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorABBBB          * binCounts_phased(5)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorBBBBB          * binCounts_phased(6)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           unphased_color_5of5 * binCounts_unphased(1)/(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           unphased_color_4of5 * binCounts_unphased(2)/(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           unphased_color_3of5 * binCounts_unphased(3)/(sum(binCounts_phased)+sum(binCounts_unphased));
					elseif (localCopyEstimate == 6)
						binCounts_phased   = zeros(1,7);
						binCounts_unphased = zeros(1,4);
						if (length(ratioData_phased) > 0)
							for i = 1:length(ratioData_phased)
								    if (ratioData_phased(i) < 1/12);  binCounts_phased(1) = binCounts_phased(1)+1;
								elseif (ratioData_phased(i) < 3/12);  binCounts_phased(2) = binCounts_phased(2)+1;
								elseif (ratioData_phased(i) < 5/12);  binCounts_phased(3) = binCounts_phased(3)+1;
								elseif (ratioData_phased(i) > 11/12); binCounts_phased(7) = binCounts_phased(7)+1;
								elseif (ratioData_phased(i) > 9/12);  binCounts_phased(6) = binCounts_phased(6)+1;
								elseif (ratioData_phased(i) > 7/12);  binCounts_phased(5) = binCounts_phased(5)+1;
								else                                  binCounts_phased(4) = binCounts_phased(4)+1;
								end;
							end;
						end;
						if (length(ratioData_unphased) > 0)
							for i = 1:length(ratioData_unphased)
								    if (ratioData_unphased(i) < 1/12);  binCounts_unphased(1) = binCounts_unphased(1)+1;
								elseif (ratioData_unphased(i) < 3/12);  binCounts_unphased(2) = binCounts_unphased(2)+1;
								elseif (ratioData_unphased(i) < 5/12);  binCounts_unphased(3) = binCounts_unphased(3)+1;
								elseif (ratioData_unphased(i) > 11/12); binCounts_unphased(1) = binCounts_unphased(1)+1;
								elseif (ratioData_unphased(i) > 9/12);  binCounts_unphased(2) = binCounts_unphased(2)+1;
								elseif (ratioData_unphased(i) > 7/12);  binCounts_unphased(3) = binCounts_unphased(3)+1;
								else                                    binCounts_unphased(4) = binCounts_unphased(4)+1;
								end;
							end;
						end;
						colorMix = colorAAAAAA         * binCounts_phased(1)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorAAAAAB         * binCounts_phased(2)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorAAAABB         * binCounts_phased(3)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorAAABBB         * binCounts_phased(4)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorAABBBB         * binCounts_phased(5)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorABBBBB         * binCounts_phased(6)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorBBBBBB         * binCounts_phased(7)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           unphased_color_6of6 * binCounts_unphased(1)/(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           unphased_color_5of6 * binCounts_unphased(2)/(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           unphased_color_4of6 * binCounts_unphased(3)/(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           unphased_color_3of6 * binCounts_unphased(4)/(sum(binCounts_phased)+sum(binCounts_unphased));
					elseif (localCopyEstimate == 7)
						binCounts_phased   = zeros(1,8);
						binCounts_unphased = zeros(1,4);
						if (length(ratioData_phased) > 0)
							for i = 1:length(ratioData_phased)
								    if (ratioData_phased(i) < 1/14);  binCounts_phased(1) = binCounts_phased(1)+1;
								elseif (ratioData_phased(i) < 3/14);  binCounts_phased(2) = binCounts_phased(2)+1;
								elseif (ratioData_phased(i) < 5/14);  binCounts_phased(3) = binCounts_phased(3)+1;
								elseif (ratioData_phased(i) < 1/2);   binCounts_phased(4) = binCounts_phased(4)+1;
								elseif (ratioData_phased(i) > 13/14); binCounts_phased(8) = binCounts_phased(8)+1;
								elseif (ratioData_phased(i) > 11/14); binCounts_phased(7) = binCounts_phased(7)+1;
								elseif (ratioData_phased(i) > 9/14);  binCounts_phased(6) = binCounts_phased(6)+1;
								elseif (ratioData_phased(i) > 1/2);   binCounts_phased(5) = binCounts_phased(5)+1;
								else
									binCounts_phased(4) = binCounts_phased(4)+0.5;
									binCounts_phased(5) = binCounts_phased(5)+0.5;
								end;
							end;
						end;
						if (length(ratioData_unphased) > 0)
							for i = 1:length(ratioData_unphased)
								    if (ratioData_unphased(i) < 1/14);  binCounts_unphased(1) = binCounts_unphased(1)+1;
								elseif (ratioData_unphased(i) < 3/14);  binCounts_unphased(2) = binCounts_unphased(2)+1;
								elseif (ratioData_unphased(i) < 5/14);  binCounts_unphased(3) = binCounts_unphased(3)+1;
								elseif (ratioData_unphased(i) > 13/14); binCounts_unphased(1) = binCounts_unphased(1)+1;
								elseif (ratioData_unphased(i) > 11/14); binCounts_unphased(2) = binCounts_unphased(2)+1;
								elseif (ratioData_unphased(i) > 9/14);  binCounts_unphased(3) = binCounts_unphased(3)+1;
								else                                    binCounts_unphased(4) = binCounts_unphased(4)+1;
								end;
							end;
						end;
						colorMix = colorAAAAAAA        * binCounts_phased(1)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorAAAAAAB        * binCounts_phased(2)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorAAAAABB        * binCounts_phased(3)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorAAAABBB        * binCounts_phased(4)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorAAABBBB        * binCounts_phased(5)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorAABBBBB        * binCounts_phased(6)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorABBBBBB        * binCounts_phased(7)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorBBBBBBB        * binCounts_phased(8)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           unphased_color_7of7 * binCounts_unphased(1)/(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           unphased_color_6of7 * binCounts_unphased(2)/(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           unphased_color_5of7 * binCounts_unphased(3)/(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           unphased_color_4of7 * binCounts_unphased(4)/(sum(binCounts_phased)+sum(binCounts_unphased));
					elseif (localCopyEstimate == 8)
						binCounts_phased   = zeros(1,9);
						binCounts_unphased = zeros(1,5);
						if (length(ratioData_phased) > 0)
							for i = 1:length(ratioData_phased)
								    if (ratioData_phased(i) < 1/16);  binCounts_phased(1) = binCounts_phased(1)+1;
								elseif (ratioData_phased(i) < 3/16);  binCounts_phased(2) = binCounts_phased(2)+1;
								elseif (ratioData_phased(i) < 5/16);  binCounts_phased(3) = binCounts_phased(3)+1;
								elseif (ratioData_phased(i) < 7/16);  binCounts_phased(4) = binCounts_phased(4)+1;
								elseif (ratioData_phased(i) > 15/16); binCounts_phased(9) = binCounts_phased(9)+1;
								elseif (ratioData_phased(i) > 13/16); binCounts_phased(8) = binCounts_phased(8)+1;
								elseif (ratioData_phased(i) > 11/16); binCounts_phased(7) = binCounts_phased(7)+1;
								elseif (ratioData_phased(i) > 9/16);  binCounts_phased(6) = binCounts_phased(6)+1;
								else                                  binCounts_phased(5) = binCounts_phased(5)+1;
								end;
							end;
						end;
						if (length(ratioData_unphased) > 0)
							for i = 1:length(ratioData_unphased)
								    if (ratioData_unphased(i) < 1/16);  binCounts_unphased(1) = binCounts_unphased(1)+1;
								elseif (ratioData_unphased(i) < 3/16);  binCounts_unphased(2) = binCounts_unphased(2)+1;
								elseif (ratioData_unphased(i) < 5/16);  binCounts_unphased(3) = binCounts_unphased(3)+1;
								elseif (ratioData_unphased(i) < 6/16);  binCounts_unphased(4) = binCounts_unphased(4)+1;
								elseif (ratioData_unphased(i) > 15/16); binCounts_unphased(1) = binCounts_unphased(1)+1;
								elseif (ratioData_unphased(i) > 13/16); binCounts_unphased(2) = binCounts_unphased(2)+1;
								elseif (ratioData_unphased(i) > 11/16); binCounts_unphased(3) = binCounts_unphased(3)+1;
								elseif (ratioData_unphased(i) > 9/16);  binCounts_unphased(4) = binCounts_unphased(4)+1;
								else                                    binCounts_unphased(5) = binCounts_unphased(5)+1;
								end;
							end;
						end;
						colorMix = colorAAAAAAAA       * binCounts_phased(1)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorAAAAAAAB       * binCounts_phased(2)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorAAAAAABB       * binCounts_phased(3)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorAAAAABBB       * binCounts_phased(4)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorAAAABBBB       * binCounts_phased(5)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorAAABBBBB       * binCounts_phased(6)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorAABBBBBB       * binCounts_phased(7)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorABBBBBBB       * binCounts_phased(8)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorBBBBBBBB       * binCounts_phased(9)  /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           unphased_color_8of8 * binCounts_unphased(1)/(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           unphased_color_7of8 * binCounts_unphased(2)/(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           unphased_color_6of8 * binCounts_unphased(3)/(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           unphased_color_5of8 * binCounts_unphased(4)/(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           unphased_color_4of8 * binCounts_unphased(5)/(sum(binCounts_phased)+sum(binCounts_unphased));
					elseif (localCopyEstimate >= 9)
						binCounts_phased   = zeros(1,10);
						binCounts_unphased = zeros(1,5);
						if (length(ratioData_phased) > 0)
							for i = 1:length(ratioData_phased)
								    if (ratioData_phased(i) < 1/18);  binCounts_phased(1)  = binCounts_phased(1) +1;
								elseif (ratioData_phased(i) < 3/18);  binCounts_phased(2)  = binCounts_phased(2) +1;
								elseif (ratioData_phased(i) < 5/18);  binCounts_phased(3)  = binCounts_phased(3) +1;
								elseif (ratioData_phased(i) < 7/18);  binCounts_phased(4)  = binCounts_phased(4) +1;
								elseif (ratioData_phased(i) < 1/2);   binCounts_phased(5)  = binCounts_phased(5) +1;
								elseif (ratioData_phased(i) > 17/18); binCounts_phased(10) = binCounts_phased(10)+1;
								elseif (ratioData_phased(i) > 15/18); binCounts_phased(9)  = binCounts_phased(9) +1;
								elseif (ratioData_phased(i) > 13/18); binCounts_phased(8)  = binCounts_phased(8) +1;
								elseif (ratioData_phased(i) > 11/18); binCounts_phased(7)  = binCounts_phased(7) +1;
								elseif (ratioData_phased(i) > 1/2);   binCounts_phased(6)  = binCounts_phased(6) +1;
								else
									binCounts_phased(5) = binCounts_phased(5)+0.5;
									binCounts_phased(6) = binCounts_phased(6)+0.5;
								end;
							end;
						end;
						if (length(ratioData_unphased) > 0)
							for i = 1:length(ratioData_unphased)
								    if (ratioData_unphased(i) < 1/18);  binCounts_unphased(1)  = binCounts_unphased(1) +1;
								elseif (ratioData_unphased(i) < 3/18);  binCounts_unphased(2)  = binCounts_unphased(2) +1;
								elseif (ratioData_unphased(i) < 5/18);  binCounts_unphased(3)  = binCounts_unphased(3) +1;
								elseif (ratioData_unphased(i) < 7/18);  binCounts_unphased(4)  = binCounts_unphased(4) +1;
								elseif (ratioData_unphased(i) > 17/18); binCounts_unphased(1)  = binCounts_unphased(1) +1;
								elseif (ratioData_unphased(i) > 15/18); binCounts_unphased(2)  = binCounts_unphased(2) +1;
								elseif (ratioData_unphased(i) > 13/18); binCounts_unphased(3)  = binCounts_unphased(3) +1;
								elseif (ratioData_unphased(i) > 11/18); binCounts_unphased(4)  = binCounts_unphased(4) +1;
								else                                    binCounts_unphased(5) = binCounts_unphased(5)+1;
								end;
							end;
						end;
						colorMix = colorAAAAAAAAA      * binCounts_phased(1 ) /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorAAAAAAAAB      * binCounts_phased(2 ) /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorAAAAAAABB      * binCounts_phased(3 ) /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorAAAAAABBB      * binCounts_phased(4 ) /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorAAAAABBBB      * binCounts_phased(5 ) /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorAAAABBBBB      * binCounts_phased(6 ) /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorAAABBBBBB      * binCounts_phased(7 ) /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorAABBBBBBB      * binCounts_phased(8 ) /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorABBBBBBBB      * binCounts_phased(9 ) /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           colorBBBBBBBBB      * binCounts_phased(10) /(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           unphased_color_9of9 * binCounts_unphased(1)/(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           unphased_color_8of9 * binCounts_unphased(2)/(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           unphased_color_7of9 * binCounts_unphased(3)/(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           unphased_color_6of9 * binCounts_unphased(4)/(sum(binCounts_phased)+sum(binCounts_unphased)) + ...
						           unphased_color_5of9 * binCounts_unphased(5)/(sum(binCounts_phased)+sum(binCounts_unphased));
					end;

					c_post =   colorMix   *   min(1,SNPs_to_fullData_ratio{chr}(chr_bin)) + ...
					           colorNoData*(1-min(1,SNPs_to_fullData_ratio{chr}(chr_bin)));
	            end;
	        else
	            c_post = colorInit;
	        end;
	        colors(chr_bin,1) = c_post(1);
	        colors(chr_bin,2) = c_post(2);
	        colors(chr_bin,3) = c_post(3);
	    end;
		
	    % standard : draw colorbars.
	    for chr_bin = 1:length(HETplot2{chr})+1;
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

		%% standard : cgh plot section.
		c_ = [0 0 0];
		fprintf(['main-plot : chr' num2str(chr) ':' num2str(length(CNVplot2{chr})) '\n']);
		for chr_bin = 1:length(CNVplot2{chr});
			x_ = [chr_bin chr_bin chr_bin-1 chr_bin-1];
			if (CNVplot2{chr}(chr_bin) == 0)
				CNVhistValue = 1;
			else
				CNVhistValue = CNVplot2{chr}(chr_bin);
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

		% standard : draw lines across plots for easier interpretation of CNV regions.
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
    
	    % standard : modify y axis limits to show annotation locations if any are provided.
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
			title([ project ' vs. (hapmap)' hapmap ' SNP/LOH map'],'Interpreter','none','FontSize',12);
	    end;
	    hold on;
	    % standard : end axes labels etc.
    
		% standard : show centromere outlines and horizontal marks.
		x1 = cen_start(chr)/bases_per_bin;
		x2 = cen_end(chr)/bases_per_bin;
		leftEnd  = 0.5*5000/bases_per_bin;
		rightEnd = (chr_size(chr) - 0.5*5000)/bases_per_bin;
		if (Centromere_format == 0)
			% standard chromosome cartoons in a way which will not cause segfaults when running via commandline.
			dx = cen_tel_Xindent; %5*5000/bases_per_bin;
			dy = cen_tel_Yindent; %maxY/10;
			% draw white triangles at corners and centromere locations.
			fill([leftEnd   leftEnd   leftEnd+dx ],       [maxY-dy   maxY      maxY],         [1.0 1.0 1.0], 'LineStyle', 'none');    % top left corner.
			fill([leftEnd   leftEnd   leftEnd+dx ],       [dy        0         0   ],         [1.0 1.0 1.0], 'LineStyle', 'none');    % bottom left corner.
			fill([rightEnd  rightEnd  rightEnd-dx],       [maxY-dy   maxY      maxY],         [1.0 1.0 1.0], 'LineStyle', 'none');    % top right corner.
			fill([rightEnd  rightEnd  rightEnd-dx],       [dy        0         0   ],         [1.0 1.0 1.0], 'LineStyle', 'none');    % bottom right corner.
			fill([x1-dx     x1        x2           x2+dx],[maxY      maxY-dy   maxY-dy  maxY],[1.0 1.0 1.0], 'LineStyle', 'none');    % top centromere.
			fill([x1-dx     x1        x2           x2+dx],[0         dy        dy       0   ],[1.0 1.0 1.0], 'LineStyle', 'none');    % bottom centromere.
			% draw outlines of chromosome cartoon.   (drawn after horizontal lines to that cartoon edges are not interrupted by horiz lines.
			plot([leftEnd   leftEnd   leftEnd+dx   x1-dx   x1        x2        x2+dx    rightEnd-dx   rightEnd   rightEnd   rightEnd-dx   x2+dx   x2   x1   x1-dx   leftEnd+dx   leftEnd],...
			     [dy        maxY-dy   maxY         maxY    maxY-dy   maxY-dy   maxY     maxY          maxY-dy    dy         0             0       dy   dy   0       0            dy     ],...
			      'Color',[0 0 0]);        
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


	    %% Linear figure draw section
	    if (Linear_display == true)
	        figure(Linear_fig);
	        Linear_width = Linear_Chr_max_width*chr_size(chr)/Linear_genome_size;
	        subplot('Position',[Linear_left Linear_base Linear_width Linear_height]);
	        Linear_left = Linear_left + Linear_width + Linear_left_chr_gap;
	        hold on;
	        title(chr_label{chr},'Interpreter','none','FontSize',10);

	        % linear : draw colorbars.
	        for chr_bin = 1:length(HETplot2{chr})+1;
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

			% linear : cgh plot section.
			c_ = [0 0 0];
			fprintf(['linear-plot : chr' num2str(chr) ':' num2str(length(CNVplot2{chr})) '\n']);
			for chr_bin = 1:length(CNVplot2{chr});
				x_ = [chr_bin chr_bin chr_bin-1 chr_bin-1];
				if (CNVplot2{chr}(chr_bin) == 0)
					CNVhistValue = 1;
				else
					CNVhistValue = CNVplot2{chr}(chr_bin);
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

			% linear : draw lines across plots for easier interpretation of CNV regions.
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
			% linear : end cgh plot section.

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
			x1 = cen_start(chr)/bases_per_bin;
			x2 = cen_end(chr)/bases_per_bin;
			leftEnd  = 0.5*5000/bases_per_bin;
			rightEnd = (chr_size(chr) - 0.5*5000)/bases_per_bin;
			if (Centromere_format == 0)
				% standard chromosome cartoons in a way which will not cause segfaults when running via commandline.
				dx = cen_tel_Xindent; %5*5000/bases_per_bin;
				dy = cen_tel_Yindent; %maxY/10;
				% draw white triangles at corners and centromere locations.
				fill([leftEnd   leftEnd   leftEnd+dx ],       [maxY-dy   maxY      maxY],         [1.0 1.0 1.0], 'linestyle', 'none');  % top left corner.
				fill([leftEnd   leftEnd   leftEnd+dx ],       [dy        0         0   ],         [1.0 1.0 1.0], 'linestyle', 'none');  % bottom left corner.
				fill([rightEnd  rightEnd  rightEnd-dx],       [maxY-dy   maxY      maxY],         [1.0 1.0 1.0], 'linestyle', 'none');  % top right corner.
				fill([rightEnd  rightEnd  rightEnd-dx],       [dy        0         0   ],         [1.0 1.0 1.0], 'linestyle', 'none');  % bottom right corner.
				fill([x1-dx     x1        x2           x2+dx],[maxY      maxY-dy   maxY-dy  maxY],[1.0 1.0 1.0], 'linestyle', 'none');  % top centromere.
				fill([x1-dx     x1        x2           x2+dx],[0         dy        dy       0   ],[1.0 1.0 1.0], 'linestyle', 'none');  % bottom centromere.
				% draw outlines of chromosome cartoon.   (drawn after horizontal lines to that cartoon edges are not interrupted by horiz lines.
				plot([leftEnd   leftEnd   leftEnd+dx   x1-dx   x1        x2        x2+dx   rightEnd-dx   rightEnd   rightEnd   rightEnd-dx   x2+dx   x2   x1   x1-dx   leftEnd+dx   leftEnd],...
				      [dy        maxY-dy   maxY         maxY    maxY-dy   maxY-dy   maxY    maxY          maxY-dy    dy         0             0       dy   dy   0       0            dy],...
				      'Color',[0 0 0]);
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

			% linear : Final formatting stuff.
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
				ylabel({project;'vs. (hapmap)';hamap}, 'Rotation', 0, 'HorizontalAlign', 'right', 'VerticalAlign', 'bottom','Interpreter','none','FontSize',5);
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
			%end final reformatting.
	        
			% shift back to main figure generation.
			figure(fig);
			hold on;
		end;
	end;
end;


%% ========================================================================
% end stuff
%==========================================================================

%% Save figures.
set(fig,'PaperPosition',[0 0 8 6]*2);
saveas(fig,        [projectDir 'fig.CNV-SNP-map.1.eps'], 'epsc');
saveas(fig,        [projectDir 'fig.CNV-SNP-map.1.png'], 'png');
set(Linear_fig,'PaperPosition',[0 0 8 0.62222222]*2);
saveas(Linear_fig, [projectDir 'fig.CNV-SNP-map.2.eps'], 'epsc');
saveas(Linear_fig, [projectDir 'fig.CNV-SNP-map.2.png'], 'png');

%% Delete figures from memory.
delete(fig);
delete(Linear_fig);

end
