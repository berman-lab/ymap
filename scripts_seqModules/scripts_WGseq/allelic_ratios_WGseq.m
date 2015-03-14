function [] = allelic_ratios_WGseq(main_dir,user,genomeUser,project,parent,hapmap,genome,ploidyEstimateString,ploidyBaseString, ...
                                   SNP_verString,LOH_verString,CNV_verString,displayBREAKS);
addpath('../');

fprintf('Generating FirePlot of SNP allelic ratio data across genome.');
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
Linear_display              = true;
Linear_displayBREAKS        = false;

projectDir = [main_dir 'users/' user '/projects/' project '/'];
genomeDir  = [main_dir 'users/' genomeUser '/genomes/' genome '/'];


[centromeres, chr_sizes, figure_details, annotations, ploidy_default] = Load_genome_information(genomeDir);
[Aneuploidy]                                                          = Load_dataset_information(projectDir);

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
end
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
		fprintf(['\t[' num2str(annotations(i).chr) ':' annotations(i).type ':' num2str(annotations(i).start) ':' ...
		         num2str(annotations(i).end) ':' annotations(i).fillcolor ':' annotations(i).edgecolor ':' num2str(annotations(i).size) ']\n']);
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


%%============================================================================================================
% Define files for processing.
%-------------------------------------------------------------------------------------------------------------
datafile = [projectDir 'putative_SNPs_v4.txt'];


%%============================================================================================================
% Preallocate data vectors the length of each chromosome.
%-------------------------------------------------------------------------------------------------------------
chr_SNP_data_positions = cell(length(chr_size),1);
chr_SNP_data_ratios    = cell(length(chr_size),1);
chr_count              = cell(length(chr_size),1);
for chrID = 1:length(chr_size)
	if (chr_in_use(chrID) == 1)
		chr_SNP_data_positions{chrID} = zeros(chr_size(chrID),1);
		chr_SNP_data_ratios{   chrID} = zeros(chr_size(chrID),1);
		chr_count{             chrID} = zeros(chr_size(chrID),1);
		chr_lines_analyzed(    chrID) = 0;
	end;
end;


%%============================================================================================================
% Process project 1 dataset.
%-------------------------------------------------------------------------------------------------------------
data = fopen(datafile, 'r');
while not (feof(data))
	dataLine = fgetl(data);
	if (length(dataLine) > 0)
		% process the loaded line into data channels.
		SNP_chr_name   = sscanf(dataLine, '%s',1);
		SNP_coordinate = sscanf(dataLine, '%s',2);   for i = 1:size(sscanf(dataLine,'%s',1),2);   SNP_coordinate(1) = [];   end;
		SNP_reference  = sscanf(dataLine, '%s',3);   for i = 1:size(sscanf(dataLine,'%s',2),2);   SNP_reference(1)  = [];   end;
		SNP_countA     = sscanf(dataLine, '%s',4);   for i = 1:size(sscanf(dataLine,'%s',3),2);   SNP_countA(1)     = [];   end;
		SNP_countT     = sscanf(dataLine, '%s',5);   for i = 1:size(sscanf(dataLine,'%s',4),2);   SNP_countT(1)     = [];   end;
		SNP_countG     = sscanf(dataLine, '%s',6);   for i = 1:size(sscanf(dataLine,'%s',5),2);   SNP_countG(1)     = [];   end;
		SNP_countC     = sscanf(dataLine, '%s',7);   for i = 1:size(sscanf(dataLine,'%s',6),2);   SNP_countC(1)     = [];   end;
		chr_num        = find(strcmp(SNP_chr_name, chr_name));
		if (length(chr_num) > 0)
			SNP_countA       = str2num(SNP_countA);
			SNP_countT       = str2num(SNP_countT);
			SNP_countG       = str2num(SNP_countG);
			SNP_countC       = str2num(SNP_countC);
			count_vector     = [SNP_countA SNP_countT SNP_countG SNP_countC];
			SNP_coordinate   = str2num(SNP_coordinate);
			chr_lines_analyzed(chr_num) = chr_lines_analyzed(chr_num)+1;
			chr_SNP_data_positions{chr_num}(chr_lines_analyzed(chr_num)) = SNP_coordinate;
			chr_SNP_data_ratios   {chr_num}(chr_lines_analyzed(chr_num)) = max(count_vector)/sum(count_vector);
			chr_count             {chr_num}(chr_lines_analyzed(chr_num)) = sum(count_vector);
		end;
	end;
end;
fclose(data);


%%============================================================================================================
% Clean up data vectors.
%-------------------------------------------------------------------------------------------------------------
for chrID = 1:length(chr_size)
	if (chr_in_use(chrID) == 1)
		chr_SNP_data_ratios{   chrID}(chr_SNP_data_positions{chrID} == 0) = [];
		chr_count{             chrID}(chr_SNP_data_positions{chrID} == 0) = [];
		chr_SNP_data_positions{chrID}(chr_SNP_data_positions{chrID} == 0) = [];
		chr_SNP_data_ratios{   chrID}(chr_count{chrID} <= 20)             = [];
		chr_SNP_data_positions{chrID}(chr_count{chrID} <= 20)             = [];
		chr_count{             chrID}(chr_count{chrID} <= 20)             = [];
	end;
end;


%%============================================================================================================
% Save processed data file.
%-------------------------------------------------------------------------------------------------------------
save([projectDir 'SNP_' SNP_verString '.all1.mat'],'chr_SNP_data_ratios','chr_SNP_data_positions','chr_count');



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
% Control variables for Candida albicans SC5314.
%--------------------------------------------------------------------------
projectDir = [main_dir 'users/' user '/projects/' project '/'];
genomeDir  = [main_dir 'users/' genomeUser '/genomes/' genome '/'];
if (strcmp(hapmap,genome) == 0)
	useHapmap = true;
	if (exist([main_dir 'users/default/hapmaps/' hapmap '/'], 'dir') == 7)
		hapmapDir = [main_dir 'users/default/hapmaps/' hapmap '/'];   % system hapmap.
	else
		hapmapDir = [main_dir 'users/' user '/hapmaps/' hapmap '/'];  % user hapmap.
	end;
else
	useHapmap = false;
end;
if (strcmp(project,parent) == 0)
	useParent = true;
	if (exist([main_dir 'users/default/projects/' parent '/'], 'dir') == 7)
		parentDir = [main_dir 'users/default/projects/' parent '/'];   % system parent.
	else
		parentDir = [main_dir 'users/' user '/projects/' parent '/'];  % user parent.
	end;
else
	useParent = false
	parentDir = projectDir;
end;


[centromeres, chr_sizes, figure_details, annotations, ploidy_default] = Load_genome_information(genomeDir);
[Aneuploidy]                                                          = Load_dataset_information(projectDir);

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
		fprintf(['\t[' num2str(annotations(i).chr) ':' annotations(i).type ':' num2str(annotations(i).start) ':' num2str(annotations(i).end) ':' ...
		         annotations(i).fillcolor ':' annotations(i).edgecolor ':' num2str(annotations(i).size) ']\n']);
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


%%=========================================================================
%%= No further control variables below. ===================================
%%=========================================================================

% Sanitize user input of euploid state.
ploidyBase = 2;
fprintf(['\nEuploid base = "' num2str(ploidyBase) '"\n']);

% basic plot parameters not defined per genome.
TickSize         = 0; % -0.005;  %negative for outside, percentage of longest chr figure.
bases_per_bin    = max(chr_size)/700;
maxY             = 50;   % number of Y-bins in 2D smoothed histogram.
cen_tel_Xindent  = 5;
cen_tel_Yindent  = maxY/5;

%define colors for colorBars plot
colorNoData = [1.0   1.0   1.0  ]; %used when no data is available for the bin.
colorInit   = [0.5   0.5   0.5  ]; %external; used in blending at ends of chr.
colorHET    = [0.0   0.0   0.0  ]; % near 1:1 ratio SNPs
colorOddHET = [0.0   1.0   0.0  ]; % Het, but not near 1:1 ratio SNPs.
colorHOM    = [1.0   0.0   0.0  ]; % Hom SNPs;

colorAB     = [0.667 0.667 0.667]; % heterozygous.
colorA      = [1.0   0.0   1.0  ]; % homozygous a:magenta.
colorB      = [0.0   1.0   1.0  ]; % homozygous b:cyan.

fprintf(['\nGenerating LOH-map figure from ''' project ''' vs. (hapmap)''' hapmap ''' data.\n']);

% Initializes vectors used to hold number of SNPs in each interpretation catagory for each chromosome region.
for chr = 1:length(chr_sizes)
	% 4 SNP interpretation catagories tracked.
	%	1 : phased ratio data.
	%	2 : unphased ratio data.
	%   3 : phased coordinate data.
	%   4 : unphased coordinate data.
	chr_length = ceil(chr_size(chr)/bases_per_bin);
	for j = 1:4
		chr_SNPdata{chr,j} = cell(1,chr_length);
	end;
	% fprintf(['0|' num2str(chr) ':' num2str(length(chr_SNPdata{chr,1})) '\n']);
end;


%%================================================================================================
% Load SNP/LOH data.
%-------------------------------------------------------------------------------------------------
%    chr_SNP_data_ratios    : allelic ratios of SNP data.
%    chr_SNP_data_positions : coordinates of SNP data.
%    chr_count              : number of chromosomes in dataset.
load([projectDir 'SNP_' SNP_verString '.all1.mat']);


%-------------------------------------------------------------------------------------------------
% Load CNV data.
%-------------------------------------------------------------------------------------------------
%    genome_CNV : string defining name of genome in use.
%    CNVplot2   : 
load ([projectDir 'Common_CNV.mat']);

%	bin_position = ceil(position/bases_per_bin);
%	for chr = 1:num_chrs                   % chr : chromosome. 
%		for i = 1:length(CNVplot2{chr});   % i   : bin along chromosome.		
%		end;
%	end;


%% -----------------------------------------------------------------------------------------
% Median normalize CNV data before figure generation.
%-------------------------------------------------------------------------------------------
CNVdata_all = [];
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		CNVdata_all = [CNVdata_all   CNVplot2{chr}];
	end;
end;
medianCNV = median(CNVdata_all)
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		CNVplot2{chr} = CNVplot2{chr}/medianCNV;
	end;
end;


%% -----------------------------------------------------------------------------------------
% Setup for main figure generation.
%------------------------------------------------------------------------------------------
% threshold for full color saturation in SNP/LOH figure.
% synced to bases_per_bin as below, or defaulted to 50.
full_data_threshold = floor(bases_per_bin/100);

fig = figure(1);
set(gcf, 'Position', [0 70 1024 600]);
largestChr = find(chr_width == max(chr_width));
largestChr = largestChr(1);


%% -----------------------------------------------------------------------------------------
% Setup for linear-view figure generation.
%-------------------------------------------------------------------------------------------
if (Linear_display == true)
	Linear_fig = figure(2);
	Linear_genome_size   = sum(chr_size);
	Linear_Chr_max_width = 0.91;               % width for all chromosomes across figure.  1.00 - leftMargin - rightMargin - subfigure gaps.
	Linear_left_start    = 0.02;               % left margin (also right margin).   (formerly 0.01)
	Linear_left_chr_gap  = 0.07/(num_chrs-1);  % gaps between chr subfigures.
	Linear_height        = 0.6;
	Linear_base          = 0.1;
	Linear_TickSize      = -0.01;  % negative for outside, percentage of longest chr figure.
	maxY                 = 50;     % number of Y-bins in 2D smoothed histogram.
	Linear_left          = Linear_left_start;
	axisLabelPosition_horiz = 0.01125;
end;
axisLabelPosition_vert = 0.01125;

%% -----------------------------------------------------------------------------------------
% Make figures
%-------------------------------------------------------------------------------------------
first_chr = true;

%% Determine statistics of data density across entire genome.
all_data        = [];
chr_mean        = zeros(1,num_chrs);
chr_mean_scaler = zeros(1,num_chrs);
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		chr_length                                  = ceil(chr_size(chr)/bases_per_bin);
		dataX                                       = ceil(chr_SNP_data_positions{chr}/bases_per_bin)';
		dataY1                                      = chr_SNP_data_ratios{chr};
		dataY2                                      = (dataY1*maxY)';
		CNVplot2_temp{chr}                          = round(CNVplot2{chr});
		CNVplot2_temp{chr}(CNVplot2_temp{chr} == 0) = 1;
		dataX_CNVcorrection                         = 1./CNVplot2_temp{chr};
		if (length(dataX) > 0)
			% 2D smoothed hisogram with correction term.
			[imageX{chr},imageY{chr},imageC{chr}, imageD{chr}] = smoothhist2D_4_Xcorrected([dataX dataX], [dataY2 (maxY-dataY2)], 0.5,[chr_length maxY],[chr_length maxY], dataX_CNVcorrection, 1,1);
			all_data = [all_data imageD{chr}];
		end;
		chr_mean(chr) = mean(imageD{chr}(:));
		fprintf(['Chromosome ' num2str(chr) ' smoothhist2D average value = ' num2str(chr_mean) '\n']);
	end;
end;
max_mean            = max(chr_mean);
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		chr_mean_scaler(chr) = max_mean/chr_mean(chr);
	end;
end;

% value statistics.
median_val    = median(all_data(:));
mean_val      = mean(all_data(:));
mode_val      = mode(all_data(:));
min_val       = min(all_data(:));
max_val       = max(all_data(:));


%% Generate chromosome figures.
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
	    figure(fig);
	    % make standard chr cartoons.
	    left   = chr_posX(chr);
	    bottom = chr_posY(chr);
	    width  = chr_width(chr);
	    height = chr_height(chr);
	    subplot('Position',[left bottom width height]);
		hold on;
		fprintf(['\tfigposition = [' num2str(left) ' | ' num2str(bottom) ' | ' num2str(width) ' | ' num2str(height) ']\n']);


		% standard : axes labels etc.
	    xlim([0,chr_size(chr)/bases_per_bin]);
    
		%% standard : modify y axis limits to show annotation locations if any are provided.
		if (length(annotations) > 0)
			ylim([-maxY/10*1.5,maxY]);
		else
			ylim([0,maxY]);
		end;
		set(gca,'YTick',[]);
		set(gca,'YTickLabel',[]);
		set(gca,'TickLength',[(TickSize*chr_size(largestChr)/chr_size(chr)) 0]); %ensures same tick size on all subfigs.

		text(-50000/5000/2*3, maxY/2,     chr_label{chr}, 'Rotation',90, 'HorizontalAlignment','center', 'VerticalAlign','bottom', 'Fontsize',20);

		set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
		set(gca,'XTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8','2.0','2.2','2.4','2.6','2.8','3.0','3.2'});

		text(axisLabelPosition_vert, maxY/4*0, '0'  ,'HorizontalAlignment','right','Fontsize',10);
		text(axisLabelPosition_vert, maxY/4*1, '1/4','HorizontalAlignment','right','Fontsize',10);
		text(axisLabelPosition_vert, maxY/4*2, '1/2','HorizontalAlignment','right','Fontsize',10);
		text(axisLabelPosition_vert, maxY/4*3, '3/4','HorizontalAlignment','right','Fontsize',10);
		text(axisLabelPosition_vert, maxY/4*4, '1'  ,'HorizontalAlignment','right','Fontsize',10);

		set(gca,'FontSize',12);
		if (chr == find(chr_posY == max(chr_posY)))
			title([ project ' allelic fraction map'],'Interpreter','none','FontSize',24);
		end;
		% standard : end axes labels etc.

		% standard : show allelic ratio data.
		chr_length                   = ceil(chr_size(chr)/bases_per_bin);
		dataX                        = ceil(chr_SNP_data_positions{chr}/bases_per_bin)';
		dataY1                       = chr_SNP_data_ratios{chr};
		dataY2                       = (dataY1*maxY)';

		CNVplot2_temp{chr}                          = round(CNVplot2{chr});
		CNVplot2_temp{chr}(CNVplot2_temp{chr} == 0) = 1;
		dataX_CNVcorrection                         = 1./CNVplot2_temp{chr};

		% removing data points in bins with less than 20 is needed for ddRADseq visualization.
		% dataX(chr_count{chr}  <= 20) = [];
		% dataY2(chr_count{chr} <= 20) = [];
		if (length(dataX) > 0)
			% 2D smoothed hisogram with correction term.
			[imageX{chr},imageY{chr},imageC{chr}, discard] = smoothhist2D_4_Xcorrected([dataX dataX], [dataY2 (maxY-dataY2)], 0.5,[chr_length maxY],[chr_length maxY], dataX_CNVcorrection, mean_val, chr_mean_scaler(chr));

			% Image correction method to de-emphasize the near homozygous data points.
			%    The square factor correction was determined empirically, from the relative amounts of data near homozygous and heterozygous.
			%    Improvements in sequencing technology that reduce sequencing error and reduce near-homozygous data will require adjusting this.
			imageC_correction          = imageC{chr}*0;
			for y = 1:maxY
				imageC_correction(y,:) = 1-abs(y-maxY/2)/(maxY/2);
			end;
			imageC{chr} = imageC{chr}.*(1+imageC_correction.^2*16);

			image(imageX{chr}, imageY{chr}, imageC{chr});
		end;
		% standard : end show allelic ratio data.

		% standard : show ChARM breakpoints.
		if (displayBREAKS == true) && (show_annotations == true)
			chr_length = ceil(chr_size(chr)/bases_per_bin);
			for segment = 2:length(chr_breaks{chr})-1
				bP = chr_breaks{chr}(segment)*chr_length;
				plot([bP bP], [(-maxY/10*2.5) 0],  'Color',[1 0 0],'LineWidth',2);
			end;
		end;

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
	        % top left corner.
	        c_ = [1.0 1.0 1.0];
	        x_ = [leftEnd   leftEnd   leftEnd+dx];
	        y_ = [maxY-dy   maxY      maxY      ];
	        f = fill(x_,y_,c_);
	        set(f,'linestyle','none');
	        % bottom left corner.
	        x_ = [leftEnd   leftEnd   leftEnd+dx];
	        y_ = [dy        0         0         ];
	        f = fill(x_,y_,c_);
	        set(f,'linestyle','none');
	        % top right corner.
	        x_ = [rightEnd   rightEnd   rightEnd-dx];
	        y_ = [maxY-dy    maxY       maxY      ];
	        f = fill(x_,y_,c_);
	        set(f,'linestyle','none');
	        % bottom right corner.
	        x_ = [rightEnd   rightEnd   rightEnd-dx];
	        y_ = [dy         0          0         ];
	        f = fill(x_,y_,c_);
	        set(f,'linestyle','none');
	        % top centromere.
	        x_ = [x1-dx   x1        x2        x2+dx];
	        y_ = [maxY    maxY-dy   maxY-dy   maxY];
	        f = fill(x_,y_,c_);
	        set(f,'linestyle','none');
	        % bottom centromere.
	        x_ = [x1-dx   x1   x2   x2+dx];
	        y_ = [0       dy   dy   0    ];
	        f = fill(x_,y_,c_);
	        set(f,'linestyle','none');
	        % draw outlines of chromosome cartoon.   (drawn after horizontal lines to that cartoon edges are not interrupted by horiz lines.
	        plot([leftEnd   leftEnd   leftEnd+dx   x1-dx   x1        x2        x2+dx   rightEnd-dx   rightEnd   rightEnd   rightEnd-dx   x2+dx   x2   x1   x1-dx   leftEnd+dx   leftEnd],...
	             [dy        maxY-dy   maxY         maxY    maxY-dy   maxY-dy   maxY    maxY          maxY-dy    dy         0             0       dy   dy   0       0            dy     ],...
	            'Color',[0 0 0]);
	    end;
	    % standard : end show centromere.
    
	    % standard : show annotation locations
	    if (show_annotations) && (length(annotations) > 0)
	        plot([leftEnd rightEnd], [-maxY/10*1.5 -maxY/10*1.5],'color',[0 0 0]);
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
	    end;
	    % standard : end show annotation locations.
		hold off;

	    %% Linear figure draw section
	    if (Linear_display == true)
	        figure(Linear_fig);
	        Linear_width = Linear_Chr_max_width*chr_size(chr)/Linear_genome_size;
	        subplot('Position',[Linear_left Linear_base Linear_width Linear_height]);
			hold on;
	        Linear_left = Linear_left + Linear_width + Linear_left_chr_gap;
	        title(chr_label{chr},'Interpreter','none','FontSize',20);

			% linear : show segmental anueploidy breakpoints.
			if (Linear_displayBREAKS == true) && (show_annotations == true)
				chr_length = ceil(chr_size(chr)/bases_per_bin);
				for segment = 2:length(chr_breaks{chr})-1
					bP = chr_breaks{chr}(segment)*chr_length;
					plot([bP bP], [(-maxY/10*2.5) 0],  'Color',[1 0 0],'LineWidth',2);
				end;
			end;

			% linear : show allelic ratio data as 2D-smoothed scatter-plot.
			image(imageX{chr}, imageY{chr}, imageC{chr});
			% linear : end show allelic ratio data.


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
				c_ = [1.0 1.0 1.0];
				% top left corner.
				x_ = [leftEnd   leftEnd   leftEnd+dx];        y_ = [maxY-dy   maxY      maxY        ];    f = fill(x_,y_,c_);    set(f,'linestyle','none');
				% bottom left corner.     
				x_ = [leftEnd   leftEnd   leftEnd+dx];        y_ = [dy        0         0           ];    f = fill(x_,y_,c_);    set(f,'linestyle','none');
				% top right corner.
				x_ = [rightEnd   rightEnd   rightEnd-dx];     y_ = [maxY-dy    maxY       maxY      ];    f = fill(x_,y_,c_);    set(f,'linestyle','none');
				% bottom right corner.
				x_ = [rightEnd   rightEnd   rightEnd-dx];     y_ = [dy         0          0         ];    f = fill(x_,y_,c_);    set(f,'linestyle','none');
				% top centromere.
				x_ = [x1-dx   x1        x2        x2+dx];     y_ = [maxY    maxY-dy   maxY-dy   maxY];    f = fill(x_,y_,c_);    set(f,'linestyle','none');
				% bottom centromere.
				x_ = [x1-dx   x1   x2   x2+dx];               y_ = [0       dy   dy   0    ];             f = fill(x_,y_,c_);    set(f,'linestyle','none');
				% draw outlines of chromosome cartoon.   (drawn after horizontal lines to that cartoon edges are not interrupted by horiz lines.
				plot([leftEnd   leftEnd   leftEnd+dx   x1-dx   x1        x2        x2+dx   rightEnd-dx   rightEnd   rightEnd   rightEnd-dx   x2+dx   x2   x1   x1-dx   leftEnd+dx   leftEnd],...
				     [dy        maxY-dy   maxY         maxY    maxY-dy   maxY-dy   maxY    maxY          maxY-dy    dy         0             0       dy   dy   0       0            dy],...
				      'Color',[0 0 0]);
			end;
			% linear : end show centromere.

			% linear : show annotation locations
			if (show_annotations) && (length(annotations) > 0)
				plot([leftEnd rightEnd], [-maxY/10*1.5 -maxY/10*1.5],'color',[0 0 0]);
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
		set(gca,'YTick',[]);
		set(gca,'YTickLabel',[]);
	        set(gca,'TickLength',[(Linear_TickSize*chr_size(largestChr)/chr_size(chr)) 0]); %ensures same tick size on all subfigs.
	        set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
		set(gca,'XTickLabel',[]);
	        if (first_chr == true)
			% This section sets the Y-axis labelling.
			text(axisLabelPosition_horiz, maxY/4*0, '0'  ,'HorizontalAlignment','right','Fontsize',10);
			text(axisLabelPosition_horiz, maxY/4*1, '1/4','HorizontalAlignment','right','Fontsize',10);
			text(axisLabelPosition_horiz, maxY/4*2, '1/2','HorizontalAlignment','right','Fontsize',10);
			text(axisLabelPosition_horiz, maxY/4*3, '3/4','HorizontalAlignment','right','Fontsize',10);
			text(axisLabelPosition_horiz, maxY/4*4, '1'  ,'HorizontalAlignment','right','Fontsize',10);
		end;
		set(gca,'FontSize',12);
		% linear : end final reformatting.

		hold off;
	        
	        % shift back to main figure generation.
	        figure(fig);
			first_chr = false;
	    end;
	end;
end;

%% Save figures.
set(fig,'PaperPosition',[0 0 8 6]*2);
saveas(fig,        [projectDir 'fig.allelic_ratio-map.b1.eps'], 'epsc');
saveas(fig,        [projectDir 'fig.allelic_ratio-map.b1.png'], 'png');
set(Linear_fig,'PaperPosition',[0 0 8 0.62222222]*2);
saveas(Linear_fig, [projectDir 'fig.allelic_ratio-map.b2.eps'], 'epsc');
saveas(Linear_fig, [projectDir 'fig.allelic_ratio-map.b2.png'], 'png');

%% Delete figures from memory.
delete(fig);
delete(Linear_fig);

%% ========================================================================
% end stuff
%==========================================================================
end
