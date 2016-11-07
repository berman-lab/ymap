function [] = allelic_ratios_WGseq(main_dir,user,genomeUser,project,parent,hapmap,genome,ploidyEstimateString,ploidyBaseString,SNP_verString,LOH_verString,CNV_verString,displayBREAKS);
addpath('../');

fprintf(['main_dir             = "' main_dir             '"\n']);
fprintf(['user                 = "' user                 '"\n']);
fprintf(['genomeUser           = "' genomeUser           '"\n']);
fprintf(['project              = "' project              '"\n']);
fprintf(['parent               = "' parent               '"\n']);
fprintf(['hapmap               = "' hapmap               '"\n']);
fprintf(['genome               = "' genome               '"\n']);
fprintf(['ploidyEstimateString = "' ploidyEstimateString '"\n']);
fprintf(['ploidyBaseString     = "' ploidyBaseString     '"\n']);

fprintf('\n\n\t*===============================================================*\n');
fprintf(    '\t| Fireplot generation in script "allelic_ratios_WGseq.m".       |\n');
fprintf(    '\t*---------------------------------------------------------------*\n');
tic;
fprintf('\t|\tGenerating FirePlot of SNP allelic ratio data across genome.\n');
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


fprintf('\t|\tDetermine if hapmap is in use.\n');
%
% For right now, ('parent' == 'hapmap') always because of earlier mixed use of variables.
% Determine if 'hapmap' is in use by checking user and system hapmap directories.
%
% Possible error case where 'parent' and 'hapmap' have same name string.
% Will be resolved with later disambiguation of parent/hapmap variable earlier in module.
%
if (exist([main_dir 'users/default/hapmaps/' hapmap '/'], 'dir') == 7)
	hapmapDir = [main_dir 'users/default/hapmaps/' hapmap '/'];   % system hapmap.
	useHapmap = true;
elseif (exist([main_dir 'users/' user '/hapmaps/' hapmap '/'], 'dir') == 7)
	hapmapDir = [main_dir 'users/' user '/hapmaps/' hapmap '/'];  % user hapmap.
	useHapmap = true;
else
	useHapmap = false;
end;


fprintf('\t|\tDetermine if parent project is in use.\n');
%
% The 'parent' will == the 'project' when no 'parent' is selected in setup.
%
if (strcmp(project,parent) == 0)
	useParent = true;
	if (exist([main_dir 'users/default/projects/' parent '/'], 'dir') == 7)
		parentDir = [main_dir 'users/default/projects/' parent '/'];   % system parent.
	else
		parentDir = [main_dir 'users/' user '/projects/' parent '/'];  % user parent.
	end;
else
	useParent = false;
	parentDir = projectDir;
end;


fprintf('\t|\tLoading dataset information.\n');
[centromeres, chr_sizes, figure_details, annotations, ploidy_default] = Load_genome_information(genomeDir);
[Aneuploidy]                                                          = Load_dataset_information(projectDir);

num_chrs = length(chr_sizes);
for chrID = 1:length(chr_sizes)
	chr_size( chrID) = 0;
	cen_start(chrID) = 0;
	cen_end(  chrID) = 0;
end;
for chrID = 1:length(chr_sizes)
	chr_size(chr_sizes(   chrID).chr) = chr_sizes(  chrID).size;
	cen_start(centromeres(chrID).chr) = centromeres(chrID).start;
	cen_end(centromeres(  chrID).chr) = centromeres(chrID).end;
end
if (length(annotations) > 0)
	fprintf(['\nAnnotations for ' genome '.\n']);
	for annoteID = 1:length(annotations)
		annotation_chr(      annoteID) = annotations(annoteID).chr;
		annotation_type{     annoteID} = annotations(annoteID).type;
		annotation_start(    annoteID) = annotations(annoteID).start;
		annotation_end(      annoteID) = annotations(annoteID).end;
		annotation_fillcolor{annoteID} = annotations(annoteID).fillcolor;
		annotation_edgecolor{annoteID} = annotations(annoteID).edgecolor;
		annotation_size(     annoteID) = annotations(annoteID).size;
		fprintf(['\t[' num2str(annotations(annoteID).chr) ':' annotations(annoteID).type ':' num2str(annotations(annoteID).start) ':' ...
		               num2str(annotations(annoteID).end) ':' annotations(annoteID).fillcolor ':' annotations(annoteID).edgecolor ':' num2str(annotations(annoteID).size) ']\n']);
	end;
end;
for figureDetailID = 1:length(figure_details)
	if (figure_details(figureDetailID).chr == 0)
		if (strcmp(figure_details(figureDetailID).label,'Key') == 1)
			key_posX   = figure_details(figureDetailID).posX;
			key_posY   = figure_details(figureDetailID).posY;
			key_width  = figure_details(figureDetailID).width;
			key_height = figure_details(figureDetailID).height;
		end;
	else
		chr_id    (figure_details(figureDetailID).chr) = figure_details(figureDetailID).chr;
		chr_label {figure_details(figureDetailID).chr} = figure_details(figureDetailID).label;
		chr_name  {figure_details(figureDetailID).chr} = figure_details(figureDetailID).name;
		chr_posX  (figure_details(figureDetailID).chr) = figure_details(figureDetailID).posX;
		chr_posY  (figure_details(figureDetailID).chr) = figure_details(figureDetailID).posY;
		chr_width (figure_details(figureDetailID).chr) = figure_details(figureDetailID).width;
		chr_height(figure_details(figureDetailID).chr) = figure_details(figureDetailID).height;
		chr_in_use(figure_details(figureDetailID).chr) = str2num(figure_details(figureDetailID).useChr);
	end;
end;

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


%%================================================================================================
% Load FASTA file name from 'reference.txt' file for project.
%-------------------------------------------------------------------------------------------------
fprintf('\t|\tLoad FASTA file name for project.\n');
userReference    = [main_dir 'users/' user '/genomes/' genome '/reference.txt'];
defaultReference = [main_dir 'users/default/genomes/' genome '/reference.txt'];
if (exist(userReference,'file') == 0)
	FASTA_string = strtrim(fileread(defaultReference));
else
	FASTA_string = strtrim(fileread(userReference));
end;
[FastaPath,FastaName,FastaExt] = fileparts(FASTA_string);


%%================================================================================================
% Preallocate data vectors the length of each chromosome.
%-------------------------------------------------------------------------------------------------
fprintf('\t|\tPreallocating data vectors the length of each chromosome.\n');
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


%%================================================================================================
% Process project 1 dataset.
%-------------------------------------------------------------------------------------------------
if (useHapmap)
	% Load only putative SNP data corresponding to hapmap loci.
	fprintf('\t|\tLoad SNP information from "trimmed_SNPs_v5.txt" file for project.\n');
	fprintf('\t|\t\t');
	datafile   = [projectDir 'trimmed_SNPs_v5.txt'];
else
	% Load all putative SNP data.
	fprintf('\t|\tLoad SNP information from "putative_SNPs_v4.txt" file for project.\n');
	fprintf('\t|\t\t');
	datafile   = [projectDir 'putative_SNPs_v4.txt'];
end;

data       = fopen(datafile, 'r');
count      = 0;
old_chr    = 0;
gap_string = '';
% reading the line before checking for end of file to avoid reading empty
% file
dataLine = fgetl(data);
while not (feof(data))
	if (length(dataLine) > 0)
		% process the loaded line into data channels.
		% if using hapmap no SNP reference is in the file so avoid reading
		% it
		if (useHapmap)
		    lineVariables = textscan(dataLine, '%s %d %d %d %d %d');
		    SNP_chr_name   = lineVariables{1}{1};
		    SNP_coordinate = lineVariables{2};
		    SNP_countA     = lineVariables{3};
		    SNP_countT     = lineVariables{4};
		    SNP_countG     = lineVariables{5};
		    SNP_countC     = lineVariables{6};
		else
		    lineVariables = textscan(dataLine, '%s %d %s %d %d %d %d');
		    SNP_chr_name   = lineVariables{1}{1};
		    SNP_coordinate = lineVariables{2};
		    SNP_reference  = lineVariables{3}{1};
		    SNP_countA     = lineVariables{4};
		    SNP_countT     = lineVariables{5};
		    SNP_countG     = lineVariables{6};
		    SNP_countC     = lineVariables{7};
		end;
		chr_num = find(strcmp(SNP_chr_name, chr_name));
		% running if it's not a comment line and the chromosome is found
		if (~strcmp(SNP_chr_name,'###') && length(find(strcmp(SNP_chr_name, chr_name))) > 0)
			count = count+1;
			if (old_chr ~= chr_num)
				fprintf(['\n\t|\t' SNP_chr_name '\n\t|\t' gap_string]);
			end;
			if (mod(count,300) == 0)
				fprintf('.');
				gap_string = [gap_string ' '];
			end;
			if (count == 24000)
				fprintf('\n\t|\t');
				count = 0;
				gap_string = '';
			end;
			count_vector     = [SNP_countA SNP_countT SNP_countG SNP_countC];
			chr_lines_analyzed(chr_num) = chr_lines_analyzed(chr_num)+1;
			chr_SNP_data_positions{chr_num}(chr_lines_analyzed(chr_num)) = SNP_coordinate;
			chr_SNP_data_ratios   {chr_num}(chr_lines_analyzed(chr_num)) = max(count_vector)/sum(count_vector);
			chr_count             {chr_num}(chr_lines_analyzed(chr_num)) = sum(count_vector);
			old_chr          = chr_num;
		end;
	end;
    % read next line
	dataLine = fgetl(data);
end;
fclose(data);

%%================================================================================================
% Clean up data vectors.
%-------------------------------------------------------------------------------------------------
fprintf('\n\t|\tClean up data vectors.\n');
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


%%================================================================================================
% Save processed SNP/LOH data file.
%-------------------------------------------------------------------------------------------------
%    chr_SNP_data_ratios    : allelic ratios of SNP data.
%    chr_SNP_data_positions : coordinates of SNP data.
%    chr_count              : number of chromosomes in dataset.
%
fprintf('\t|\tSave processed SNP/LOH data to file "SNP_v4.all1.mat" for project.\n');
save([projectDir 'SNP_' SNP_verString '.all1.mat'],'chr_SNP_data_ratios','chr_SNP_data_positions','chr_count');


%%================================================================================================
% Setup basic figure parameters.
%-------------------------------------------------------------------------------------------------
fprintf('\t|\tDefine basic figure parameters, not specific to genome.\n');
% basic plot parameters not defined per genome.
TickSize         = 0; % -0.005;  %negative for outside, percentage of longest chr figure.
bases_per_bin    = max(chr_size)/700;
maxY             = 50;   % number of Y-bins in 2D smoothed histogram.
cen_tel_Xindent  = 5;
cen_tel_Yindent  = maxY/5;
largestChr       = find(chr_width == max(chr_width));
largestChr       = largestChr(1);


%%================================================================================================
% Setup for main-view figure generation.
%-------------------------------------------------------------------------------------------------
% load size definitions
[linear_fig_height,linear_fig_width,Linear_left_start,Linear_chr_gap,Linear_Chr_max_width,Linear_height...
    ,Linear_base,rotate,linear_chr_font_size,linear_axis_font_size,linear_gca_font_size,stacked_fig_height,...
    stacked_fig_width,stacked_chr_font_size,stacked_title_size,stacked_axis_font_size,...
    gca_stacked_font_size,stacked_copy_font_size,max_chrom_label_size] = Load_size_info(chr_in_use,num_chrs,chr_label,chr_size);

fprintf('\t|\tInitialize main figure.\n');
fig = figure(1);

%%================================================================================================
% Setup for linear-view figure generation.
%-------------------------------------------------------------------------------------------------
if (Linear_display == true)
	fprintf('\t|\tInitialize linear figure.\n');
	Linear_fig = figure(2);
	Linear_genome_size   = sum(chr_size);
	Linear_TickSize      = -0.01;  % negative for outside, percentage of longest chr figure.
	maxY                 = 50;     % number of Y-bins in 2D smoothed histogram.
	Linear_left          = Linear_left_start;
	axisLabelPosition_horiz = 0.01125;
end;
axisLabelPosition_vert = 0.01125;


%%================================================================================================
% Make figures
%-------------------------------------------------------------------------------------------------
first_chr = true;

%% Determine statistics of data density across entire genome.
fprintf('\t|\tDetermine statistics of data density for chromosomes.\n');
all_data        = [];
chr_mean        = zeros(1,num_chrs);
chr_mean_scaler = zeros(1,num_chrs);
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		chr_length          = ceil(chr_size(chr)/bases_per_bin);
		dataX               = ceil(chr_SNP_data_positions{chr}/bases_per_bin)';
		dataY1              = chr_SNP_data_ratios{chr};
		dataY2              = (dataY1*maxY)';
		dataX_CNVcorrection = ones(1,chr_length);
		if (length(dataX) > 0)
			% 2D smoothed hisogram with correction term.
			[imageX{chr},imageY{chr},imageC{chr}, imageD{chr}] = smoothhist2D_4_Xcorrected([dataX dataX 0 chr_length], [dataY2 (maxY-dataY2) 0 0], 0.5,[chr_length maxY],[chr_length maxY], dataX_CNVcorrection, 1,1);
			all_data = [all_data imageD{chr}];
		end;
		if (length(dataX) > 0)
		    chr_mean(chr) = mean(imageD{chr}(:));
		else
		    chr_mean(chr) = 0;
		end;
		fprintf(['\t|\t\tChr' num2str(chr) ' smoothhist2D average value = ' num2str(chr_mean) '\n']);
	end;
end;
max_mean = max(chr_mean);
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1 && chr_mean(chr) ~= 0)
		chr_mean_scaler(chr) = max_mean/chr_mean(chr);
	else
		chr_mean_scaler(chr) = 0;
	end;
end;
median_val = median(all_data(:));
mean_val   = mean(all_data(:));
mode_val   = mode(all_data(:));
min_val    = min(all_data(:));
max_val    = max(all_data(:));


%% Generate chromosome figures.
fprintf('\t|\tGenerate final chromosome figures.\n');
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
		text(-50000/5000/2*3, maxY/2,     chr_label{chr}, 'Rotation',90, 'HorizontalAlignment','center', 'VerticalAlign','bottom', 'Fontsize',stacked_chr_font_size);
		set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
		set(gca,'XTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8','2.0','2.2','2.4','2.6','2.8','3.0','3.2'});
		text(axisLabelPosition_vert, maxY/4*0, '0'  ,'HorizontalAlignment','right','Fontsize',stacked_axis_font_size);
		text(axisLabelPosition_vert, maxY/4*1, '1/4','HorizontalAlignment','right','Fontsize',stacked_axis_font_size);
		text(axisLabelPosition_vert, maxY/4*2, '1/2','HorizontalAlignment','right','Fontsize',stacked_axis_font_size);
		text(axisLabelPosition_vert, maxY/4*3, '3/4','HorizontalAlignment','right','Fontsize',stacked_axis_font_size);
		text(axisLabelPosition_vert, maxY/4*4, '1'  ,'HorizontalAlignment','right','Fontsize',stacked_axis_font_size);

		set(gca,'FontSize',gca_stacked_font_size);
		if (chr == find(chr_posY == max(chr_posY)))
			title([ project ' allelic fraction map'],'Interpreter','none','FontSize',stacked_title_size);
		end;
		% standard : end axes labels etc.

		% standard : show allelic ratio data.
		chr_length                   = ceil(chr_size(chr)/bases_per_bin);
		dataX                        = ceil(chr_SNP_data_positions{chr}/bases_per_bin)';
		dataY1                       = chr_SNP_data_ratios{chr};
		dataY2                       = (dataY1*maxY)';
		dataX_CNVcorrection          = ones(1,chr_length);;

		if (length(dataX) > 0)
			% 2D smoothed hisogram with correction term.
			fprintf(['\t|\t\tGenerating chr' num2str(chr) ' final smoothed 2D histogram.\n']);
			[imageX{chr},imageY{chr},imageC{chr}, discard] = smoothhist2D_4_Xcorrected([dataX dataX 0 chr_length], [dataY2 (maxY-dataY2) 0 0], 0.5,[chr_length maxY],[chr_length maxY], dataX_CNVcorrection, mean_val, chr_mean_scaler(chr));

			fprintf('\t|\t\t\tDe-emphasizing near-homozygous data.\n');
			% Image correction method to de-emphasize the near homozygous data points.
			%    The square factor correction was determined empirically, from the relative amounts of data near homozygous and heterozygous.
			%    Improvements in sequencing technology that reduce sequencing error and reduce near-homozygous data will require adjusting this.
			imageC_correction          = imageC{chr}*0;
			for y = 1:maxY
				imageC_correction(y,:) = 1-abs(y-maxY/2)/(maxY/2);
			end;
			imageC{chr} = imageC{chr}.*(1+imageC_correction.^2*16);

			fprintf('\t|\t\t\tDrawing histogram to figure.\n');
			image(imageX{chr}, imageY{chr}, imageC{chr});
		end;
		% standard : end show allelic ratio data.

		% standard : show ChARM breakpoints.
		if (displayBREAKS == true) && (show_annotations == true)
			fprintf('\t|\t\t\tShow ChARM breakpoints.\n');
			chr_length = ceil(chr_size(chr)/bases_per_bin);
			for segment = 2:length(chr_breaks{chr})-1
				bP = chr_breaks{chr}(segment)*chr_length;
				plot([bP bP], [(-maxY/10*2.5) 0],  'Color',[1 0 0],'LineWidth',2);
			end;
		end;

		% standard : show centromere outlines and horizontal marks.
		fprintf('\t|\t\t\tDraw centromere and horizontal lines.\n');
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
    
		%% standard : show annotation locations
		fprintf('\t|\t\t\tShow annotation locations.\n');
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
			Linear_left = Linear_left + Linear_width + Linear_chr_gap;
            
			% linear : show segmental anueploidy breakpoints.
			if (Linear_displayBREAKS == true) && (show_annotations == true)
				fprintf('\t|\t\t\tShow ChARM breakpoints on linear figure.\n');
				chr_length = ceil(chr_size(chr)/bases_per_bin);
				for segment = 2:length(chr_breaks{chr})-1
					bP = chr_breaks{chr}(segment)*chr_length;
					plot([bP bP], [(-maxY/10*2.5) 0],  'Color',[1 0 0],'LineWidth',2);
				end;
			end;

			% linear : show allelic ratio data as 2D-smoothed scatter-plot.
			fprintf('\t|\t\t\tDraw 2D smoothed histogram of allelic ratio data in linear figure.\n');
			% display only if processing succeeded (variables will no be
			% present if the data is zero)
			if (exist('imageX') && exist('imageY') && exist('imageC'))
				image(imageX{chr}, imageY{chr}, imageC{chr});
			end;
			% linear : end show allelic ratio data.

			% linear : show centromere.
			fprintf('\t|\t\t\tDraw centromere in linear figure.\n');
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

			% linear : show annotation locations.
			if (show_annotations) && (length(annotations) > 0)
				fprintf('\t|\t\t\tShow annotation locations in linear figure.\n');
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
				text(axisLabelPosition_horiz, maxY/4*0, '0'  ,'HorizontalAlignment','right','Fontsize',linear_axis_font_size);
				text(axisLabelPosition_horiz, maxY/4*1, '1/4','HorizontalAlignment','right','Fontsize',linear_axis_font_size);
				text(axisLabelPosition_horiz, maxY/4*2, '1/2','HorizontalAlignment','right','Fontsize',linear_axis_font_size);
				text(axisLabelPosition_horiz, maxY/4*3, '3/4','HorizontalAlignment','right','Fontsize',linear_axis_font_size);
				text(axisLabelPosition_horiz, maxY/4*4, '1'  ,'HorizontalAlignment','right','Fontsize',linear_axis_font_size);
			end;
			set(gca,'FontSize',linear_gca_font_size);
			% linear : end final reformatting.
			% adding title in the middle of the cartoon
			% note: adding title is done in the end since if placed upper
			% in the code somehow the plot function changes the title position			
			if (rotate == 0 && chr_size(chr) ~= 0 )
				title(chr_label{chr},'Interpreter','none','FontSize',linear_chr_font_size,'Rotation',rotate);
			else
				text((chr_size(chr)/bases_per_bin)/2,maxY+0.5,chr_label{chr},'Interpreter','none','FontSize',linear_chr_font_size,'Rotation',rotate);
			end;

			hold off;
	        
			% shift back to main figure generation.
			figure(fig);
			first_chr = false;
		end;
	end;
end;

%% Save figures.
% commenting out stacked figure since it's not diplayed, left for debugging
%{
set(fig,'PaperPosition',[0 0 stacked_fig_width stacked_fig_height]);
fprintf('\t|\tSaving standard figure in EPS format.\n');
saveas(fig,        [projectDir 'fig.allelic_ratio-map.b1.eps'], 'epsc');
fprintf('\t|\tSaving standard figure in PNG format.\n');
saveas(fig,        [projectDir 'fig.allelic_ratio-map.b1.png'], 'png');
%}
set(Linear_fig,'PaperPosition',[0 0 linear_fig_width linear_fig_height]);
fprintf('\t|\tSaving linear figure in EPS format.\n');
saveas(Linear_fig, [projectDir 'fig.allelic_ratio-map.b2.eps'], 'epsc');
fprintf('\t|\tSaving linear figure in PNG format.\n');
saveas(Linear_fig, [projectDir 'fig.allelic_ratio-map.b2.png'], 'png');

%% Delete figures from memory.
delete(fig);
delete(Linear_fig);

time_end = toc;
fprintf('\t|\t%d min, %f sec.\n',floor(time_end/60),rem(time_end,60));
fprintf('\t*---------------------------------------------------------------*\n');
fprintf('\t| Fireplot generation in "allelic_ratios_WGseq.m" completed.    |\n');
fprintf('\t*===============================================================*\n');
end
