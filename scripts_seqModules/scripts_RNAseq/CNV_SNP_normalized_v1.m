function [] = CNV_SNP_normalized_v1(projectName1_parent,projectName2_child,genome_data,ploidyString, ...
                  ave_copy_num,CNV_verString,SNP_verString,LOH_verString,workingDir,figureDir,displayBREAKS);
addpath('../);

%% ========================================================================
% Generate CGH-type figures from RNAseq data, using a reference dataset to correct for genome position-dependant biases.
%==========================================================================
%    Centromere_format          : Controls how centromeres are depicted.   [0..2]   '2' is pinched cartoon default.
%    bases_per_bin              : Controls bin sizes for CGH fractions of plot.
%    scale_type                 : 'Ratio' or 'Log2Ratio' y-axis scaling of copy number.
%                                 'Log2Ratio' does not properly scale CGH data by ploidy.
Centromere_format           = 0;
Yscale_nearest_even_ploidy  = true;
CNVhistPlot                 = true;
SNPhistPlot                 = true;		if CNVhistPlot; SNPhistPlot = false; end;
ChrNum                      = true;
show_annotations            = true;
colorBars                   = true;
blendColorBars              = false;
Linear_display              = true;

%%=========================================================================
% Control variables.
%--------------------------------------------------------------------------
% Defines chr sizes in bp. (diploid total=28,567,7888)
% Defines centromere locations in bp.
% Defines annotation locations in bp.

[centromeres, chr_sizes, figure_details, annotations, ploidy_default] = Load_genome_information(genomeDir);
[Aneuploidy]                                                          = Load_dataset_information(projectName2_child);

for i = 1:length(chr_sizes)
    chr_size(chr_sizes(i).chr)    = chr_sizes(i).size;
end;
for i = 1:length(centromeres)
    cen_start(centromeres(i).chr) = centromeres(i).start;
    cen_end(centromeres(i).chr)   = centromeres(i).end;
end;
if (length(annotations) > 0)
    fprintf(['\nAnnotations for ' genome_data '.\n']);
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
% bases_per_bin			= 5000;
bases_per_bin			= max(chr_size)/700;
chr_length_scale_multiplier	= 1/bases_per_bin;

%%=========================================================================
%%= No further control variables below. ===================================
%%=========================================================================
if (strcmp(projectName1_parent,projectName2_child) == 1)
    fprintf(['\nGenerating CNV-SNP map figure from ' projectName1_parent ' genome data (CNV normalized by reference).\n']);
else
    fprintf(['\nGenerating CNV-LOH map figure from ' projectName1_parent '(parent) and ' projectName2_child '(child) genome data (CNV normalized by reference).\n']);
end;

%% Load CNV figure data from earlier analysis.
if (exist([workingDir 'matlab_dir/' projectName2_child '.normalized_CNV_' CNV_verString '.ploidy_' ploidyString '.mat'],'file') == 0)
    fprintf('\nMAT(CNV_normalized) file not found, exiting.\n');
    exit;
else
    fprintf('\nMAT(CNV_normalized) file found, loading.\n');
    load([workingDir 'matlab_dir/' projectName2_child '.normalized_CNV_' CNV_verString '.ploidy_' ploidyString '.mat']);
    % 'chr_CNVdata','CGD_chr_CNVdata'.
end;

%% Load SNP/LOH figure data from earlier analysis.
% file names differ if the SNPs of one or two strains are being examined.
if (strcmp(projectName1_parent,projectName2_child) == 1)
    fileName_SNP = [projectName2_child '.SNP_' SNP_verString '.mat'];
else
    fileName_SNP = [projectName1_parent '->' projectName2_child '.LOH_' LOH_verString '.mat'];
end;

if (exist([workingDir 'matlab_dir/' fileName_SNP],'file') == 0)
    fprintf('\nMAT(SNP) file not found, exiting.\n');
    exit;
else
    fprintf('\nMAT(SNP) file found, loading.\n');
    load([workingDir 'matlab_dir/' fileName_SNP]);
    % 'chr_SNPdata','ave_copy_number'.
end;

% basic plot parameters not defined per genome.
TickSize         = -0.005;  %negative for outside, percentage of longest chr figure.
bases_per_bin    = max(chr_size)/700;
maxY             = ploidyBase*2;
cen_tel_Xindent  = 5;
cen_tel_Yindent  = maxY/5;

%define colors for colorBars plot  
colorNoData = [1.0   1.0   1.0  ]; %used when no data is available for the bin.
colorInit   = [0.5   0.5   0.5  ]; %external; used in blending at ends of chr.
colorHET    = [0.0   0.0   1.0  ]; % near 1:1 ratio SNPs
colorOddHET = [0.0   1.0   0.0  ]; % Het, but not near 1:1 ratio SNPs.
colorHOM    = [1.0   0.0   0.0  ]; % Hom SNPs;

%% -----------------------------------------------------------------------------------------
% Setup for main figure generation.
%-------------------------------------------------------------------------------------------
fig = figure(1);
set(gcf, 'Position', [0 70 1024 600]);

%% CNV pre-figure calculations.

all_CNV_data = [];
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		for pos = 1:length(chr_CNVdata{chr,1})
			% Plot the sum of the data in each region, divided by the number of data points in each region.
			% Then divided by this value calculated for SC5314 reference data.
			if (chr_CNVdata{chr,2}(pos) == 0) || (chr_refCNVdata{chr,2}(pos) == 0)
				% No data elements => null value is plotted.
				CNVplot{chr}(pos) = 0;
			else
				% Sum of data elements is divided by the number of data elements.
				normalizationFactor = chr_refCNVdata{chr,1}(pos)/chr_refCNVdata{chr,2}(pos);
				if (normalizationFactor == 0)
					CNVplot{chr}(pos) = 0;
				else
					CNVplot{chr}(pos) = chr_CNVdata{chr,1}(pos)/chr_CNVdata{chr,2}(pos)/normalizationFactor;
				end;
			end;
		end;
		chr_max(chr) = max(CNVplot{chr});
		chr_med(chr) = median(CNVplot{chr});
		all_CNV_data = [all_CNV_data CNVplot{chr}];
	end;
end;

%max_count     = max(chr_max);
%median_count  = sum(chr_med)/length(chr_med);
min_count     = min(all_CNV_data);
max_count     = max(all_CNV_data);
mean_count    = mean(all_CNV_data);
median_count  = median(all_CNV_data);
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		CNVplot2{chr} = CNVplot{chr}/median_count;
	end;
end;

fprintf( 'CNV copy number stats:\n');
fprintf(['\tmin    = ' num2str(min_count)    '\n']);
fprintf(['\tmax    = ' num2str(max_count)    '\n']);
fprintf(['\tmean   = ' num2str(mean_count)   '\n']);
fprintf(['\tmedian = ' num2str(median_count) '\n']);

ploidy = str2num(ploidyString);
fprintf(['Ploidy string = "' ploidyString '"\n']);

% Test chr break point.
% Aneuploidy(1).chr   = 3;
% Aneuploidy(1).break = 0.5;
[chr_breaks, chrCopyNum, ploidyAdjust] = FindChrSizes_4(Aneuploidy,CNVplot2,ploidy,num_chrs,chr_in_use);

fprintf('\n');
largestChr = find(chr_width == max(chr_width));

%% Determine average copy number.
data_length = 0;
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		data_length = data_length + length(chr_SNPdata{chr,1});
	end;
end;
copy_num_data = zeros(1,data_length);
j = 0;
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		for i = 1:length(chr_SNPdata{chr,1})
			j = j+1;
			copy_num_data(j) = chr_SNPdata{chr,1}(i) + chr_SNPdata{chr,2}(i) + chr_SNPdata{chr,3}(i);
		end;
	end;
end;

min_copy_num    = min(copy_num_data);
max_copy_num    = max(copy_num_data);
mean_copy_num   = mean(copy_num_data);
median_copy_num = median(copy_num_data);

fprintf( 'SNP copy number stats:\n');
fprintf(['\tmin    = ' num2str(min_copy_num)    '\n']);
fprintf(['\tmax    = ' num2str(max_copy_num)    '\n']);
fprintf(['\tmean   = ' num2str(mean_copy_num)   '\n']);
fprintf(['\tmedian = ' num2str(median_copy_num) '\n']);

%% SNP pre-figure calculations.
data_mode = 3;
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		fprintf(['chr = ' num2str(chr) ':' '\n']);
		if (data_mode == 1)
		elseif (data_mode == 2)
		elseif (data_mode == 3)
			% Details from LOH_v3a.m :
			TOTplot{chr}                                  = chr_SNPdata{chr,1}+chr_SNPdata{chr,2}+chr_SNPdata{chr,3};  % TOT data
			TOTave{chr}                                   = sum(TOTplot{chr})/length(TOTplot{chr});
			TOTplot2{chr}                                 = TOTplot{chr}/median_copy_num;
			TOTplot2{chr}(TOTplot2{chr} > 1)              = 1;
			TOTave2{chr}                                  = sum(TOTplot2{chr})/length(TOTplot2{chr});

			HETplot{chr}                                  = chr_SNPdata{chr,1};  % HET data
			HETave{chr}                                   = sum(HETplot{chr})/length(HETplot{chr});
			HETplot2{chr}                                 = HETplot{chr}/median_copy_num;
%			HETplot2{chr}(HETplot2{chr} > 1)              = 1;
			HETave2{chr}                                  = sum(HETplot2{chr})/length(HETplot2{chr});

			oddHETplot{chr}                               = chr_SNPdata{chr,2};  % oddHET data
			oddHETave{chr}                                = sum(oddHETplot{chr})/length(oddHETplot{chr});
			oddHETplot2{chr}                              = oddHETplot{chr}/median_copy_num;
%			oddHETplot2{chr}(oddHETplot2{chr} > 1)        = 1;

			HOMplot{chr}                                  = chr_SNPdata{chr,3};  % HOM data
			HOMave{chr}                                   = sum(HOMplot{chr})/length(HOMplot{chr});  
			HOMplot2{chr}                                 = HOMplot{chr}/median_copy_num;
%			HOMplot2{chr}(HOMplot2{chr} > 1)              = 1;
		end;
	end;
end;
fprintf('\n');
largestChr = find(chr_width == max(chr_width));


%% -----------------------------------------------------------------------------------------
% Setup for linear-view figure generation. 
%-------------------------------------------------------------------------------------------
if (Linear_display == true)
	Linear_fig           = figure(2);
	Linear_genome_size   = sum(chr_size);

	Linear_Chr_max_width = 0.91;               % width for all chromosomes across figure.  1.00 - leftMargin - rightMargin - subfigure gaps.
	Linear_left_start    = 0.01;               % left margin (also right margin).
	Linear_left_chr_gap  = 0.07/(num_chrs-1);  % gaps between chr subfigures.

	Linear_height        = 0.6;
	Linear_base          = 0.1;
	Linear_TickSize      = -0.01;  %negative for outside, percentage of longest chr figure.
	Linear_maxY          = 10;
	Linear_left          = Linear_left_start;
end;

%% Initialize copy number string.
stringChrCNVs = '';


%% -----------------------------------------------------------------------------------------
% Median normalize CNV data before figure generation.
%-------------------------------------------------------------------------------------------
% Gather CGH data for LOWESS fitting.
CNVdata_all = [];
for chr = 1:num_chrs
    if (chr_in_use(chr) == 1)
        CNVdata_all = [CNVdata_all   CNVplot2{chr}];
    end;
end;
CNVdata_all(CNVdata_all == 0) = [];
medianCNV = median(CNVdata_all)
for chr = 1:num_chrs
    if (chr_in_use(chr) == 1)
        CNVplot2{chr} = CNVplot2{chr}/medianCNV;
    end;
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
		fprintf(['figposition = [' num2str(left) ' | ' num2str(bottom) ' | ' num2str(width) ' | ' num2str(height) ']\t']);
		hold on;

		%% SNP(LOH) plot section.
		c_prev = colorInit;
		c_post = colorInit;
		c_     = c_prev;

		infill = zeros(1,length(HETplot2{chr}));
		colors = [];

		% determines the color of each bin.
		for i = 1:length(TOTplot2{chr})+1;
			if (i-1 < length(TOTplot2{chr}))
				c_tot_post = TOTplot2{chr}(i)+TOTplot2{chr}(i);
				if (c_tot_post == 0)
					c_post = colorNoData;
				else
					%	c_post = 	colorHET*HETplot2{chr}(i) + ...
					%				colorHOM*HOMplot2{chr}(i) + ...
					%				colorNoData*(1-min([HETplot2{chr}(i)+HOMplot2{chr}(i) 1]));
					%	colorMix =	colorHET   *HETplot2   {chr}(i)/TOTplot2{chr}(i) + ...
					%				colorOddHET*oddHETplot2{chr}(i)/TOTplot2{chr}(i) + ...
					%				colorHOM   *HOMplot2   {chr}(i)/TOTplot2{chr}(i);
					%	colorMix =	colorHET   *   HETplot2{chr}(i)/TOTplot2{chr}(i) + ...
					%				colorOddHET*oddHETplot2{chr}(i)/TOTplot2{chr}(i) + ...
					%				colorHOM   *   HOMplot2{chr}(i)/TOTplot2{chr}(i);
					%	c_post =	colorMix   *   min(1,TOTplot2{chr}(i)) + ...
					%				colorNoData*(1-min(1,TOTplot2{chr}(i)));
					%	colorNoData*(1-min([HETplot2{chr}(i)+oddHETplot2{chr}(i)+HOMplot2{chr}(i) 1]));

					% color equations for low density RNAseq data.
					% colorMix is a weighted average of the three colors for each category.

					colorMix =	colorHET   *   HETplot2{chr}(i)/TOTplot2{chr}(i) + ...
								colorOddHET*oddHETplot2{chr}(i)/TOTplot2{chr}(i) + ...
								colorHOM   *   HOMplot2{chr}(i)/TOTplot2{chr}(i);

					%	if (SNPhistPlot == true )
					%		threshold = median_copy_num*4;

					threshold = 4;
					c_post =	colorMix   *   max(0,min(1,HETplot2{chr}(i)/threshold)) + ...
								colorNoData*(1-max(0,min(1,HETplot2{chr}(i)/threshold)));

					%	elseif (CNVhistPlot == true)
					%		threshold = 25;
					%		c_post =	colorMix   *   max(0,min(1,HETplot2{chr}(i)/threshold)) + ...
					%					colorNoData*(1-max(0,min(1,HETplot2{chr}(i)/threshold)));
					%	else
					%		c_post =	colorMix   *   max(0,min(1,TOTplot2{chr}(i))) + ...
					%					colorNoData*(1-max(0,min(1,TOTplot2{chr}(i))));
					%	end;
					% color range is [0..1]; min/max statements limits to range.
				end;
			else
				c_post = colorInit;
			end;
			colors(i,1) = c_post(1);
			colors(i,2) = c_post(2);
			colors(i,3) = c_post(3);
		end;

	    % draw colorbars.
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
			if (c_(1) < 0); c_(1) = 0; end;
			if (c_(2) < 0); c_(2) = 0; end;
			if (c_(3) < 0); c_(3) = 0; end;
			if isnan(c_(1)); c_(1) = 1; end;
			if isnan(c_(2)); c_(2) = 1; end;
			if isnan(c_(3)); c_(3) = 1; end;
			if (blendColorBars == false)
				f = fill(x_,y_,c_);
			else
				f = fill(x_,y_,c_/2+c_prev/4+c_post/4);
			end;
			c_prev = c_;
			c_     = c_post;
			set(f,'linestyle','none');
		end;

		%% CGH(CNV) plot section.
		c_ = [0 0 0];
		fprintf(['chr' num2str(chr) ':' num2str(length(CNVplot2{chr})) '\n']);
		for i = 1:length(CNVplot2{chr});
			x_ = [i i i-1 i-1];
			% ploidy_factor : is used to adjust default ploidy in figure generation.  Haploid default ploidy position is set to 1/4 of the
			% scale of 1..4 while diploid position is set to 1/2 of the scale.
			if (ploidy_default == 1)
				ploidy_factor = 4;
			elseif (ploidy_default == 2)
				ploidy_factor = 2;
			end;
			if (ploidy == 0) % no value; no scaling.
				y_ = [maxY/ploidy_factor CNVplot2{chr}(i)*maxY/2 CNVplot2{chr}(i)*maxY/2 maxY/ploidy_factor];
			else
				if (Yscale_nearest_even_ploidy == true)
					if (abs(2-ploidy) < abs(4-ploidy)) % nearer to diploid
						y_ = [maxY/ploidy_factor CNVplot2{chr}(i)*ploidy/2*maxY/2 CNVplot2{chr}(i)*ploidy/2*maxY/2 maxY/ploidy_factor];
					else
						y_ = [maxY/ploidy_factor CNVplot2{chr}(i)*ploidy/4*maxY/2 CNVplot2{chr}(i)*ploidy/4*maxY/2 maxY/ploidy_factor];
					end;
				else
					y_ = [maxY/ploidy_factor CNVplot2{chr}(i)*ploidy/2*maxY/2 CNVplot2{chr}(i)*ploidy/2*maxY/2 maxY/ploidy_factor];
				end;
			end;

			% makes a blackbar for each bin.
			f = fill(x_,y_,c_);
			set(f,'linestyle','none');
		end;
		x2 = chr_size(chr)*chr_length_scale_multiplier;
		plot([0; x2], [maxY/2; maxY/2],'color',[0 0 0]);  % 2n line.
    
		%% draw lines across plots for easier interpretation of CNV regions.
		if (ploidy_default == 2)
			line([0 x2], [maxY/4*3 maxY/4*3],'Color',[0.85 0.85 0.85]);
			line([0 x2], [maxY/4   maxY/4],  'Color',[0.85 0.85 0.85]);
		elseif (ploidy_default == 1)
			line([0 x2], [maxY/4*3 maxY/4*3],'Color',[0.85 0.85 0.85]);
			line([0 x2], [maxY/2   maxY/2],  'Color',[0.85 0.85 0.85]);
		end;
		%% end cgh plot section.
    
		%axes labels etc.
		hold off;
		xlim([0,chr_size(chr)*chr_length_scale_multiplier]);
    
		%% modify y axis limits to show annotation locations if any are provided.
		if (length(annotations) > 0)
			ylim([-1.5,maxY]);
		else
			ylim([0,maxY]);
		end;

		set(gca,'YTick',[]);
		set(gca,'TickLength',[(TickSize*chr_size(largestChr)/chr_size(chr)) 0]); %ensures same tick size on all subfigs.
		ylabel(chr_label{chr}, 'Rotation', 90, 'HorizontalAlign', 'center', 'VerticalAlign', 'bottom');
		set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
		set(gca,'XTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8','2.0','2.2','2.4','2.6','2.8','3.0','3.2'});
    
		set(gca,'YTick',[maxY/4 maxY/2 maxY/4*3 maxY]);
		if (Yscale_nearest_even_ploidy == true)
			if (abs(2-ploidy) < abs(4-ploidy)) % nearer to diploid
				set(gca,'YTickLabel',{'1','2','3','4'});
			else % nearer to tetraploid
				set(gca,'YTickLabel',{'2','4','6','8'});
			end;
		else
			set(gca,'YTickLabel',{'1','2','3','4'});
		end;

		set(gca,'FontSize',6);
		if (chr == find(chr_posY == max(chr_posY)))
			title([ projectName2_child ' CNV map'],'Interpreter','none','FontSize',12);
		end;

		hold on;
		%end axes labels etc.

		%show segmental anueploidy breakpoints.
		if (displayBREAKS == true)
			for segment = 2:length(chr_breaks{chr})-1
				bP = chr_breaks{chr}(segment)*length(CNVplot2{chr});
				c_ = [0 0 1];
				x_ = [bP bP bP-1 bP-1];
				y_ = [0 maxY maxY 0];
				f = fill(x_,y_,c_);
				set(f,'linestyle','none');
			end;
		end;
    
		%show centromere.
		x1 = cen_start(chr)*chr_length_scale_multiplier;
		x2 = cen_end(chr)*chr_length_scale_multiplier;
		leftEnd  = 0.5*(5000/bases_per_bin);
		rightEnd = chr_size(chr)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
		if (Centromere_format == 0)
			% standard chromosome cartoons in a way which will not cause segfaults when running via commandline.
			dx = cen_tel_Xindent;
			dy = cen_tel_Yindent;
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
		%end show centromere.
    
		%show annotation locations
		if (show_annotations) && (length(annotations) > 0)
			plot([leftEnd rightEnd], [-1.5 -1.5],'color',[0 0 0]);
			hold on;
			annotation_location = (annotation_start+annotation_end)./2;
			for i = 1:length(annotation_location)
				if (annotation_chr(i) == chr)
					annotationloc = annotation_location(i)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
					plot(annotationloc,-1.5,'k:o','MarkerEdgeColor',annotation_edgecolor{i},'MarkerFaceColor',annotation_fillcolor{i},'MarkerSize',annotation_size(i));
					% plot(annotationloc,-1.5,'k:o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);
				end;
			end;
			hold off;
		end;
		%end show annotation locations.

		% make CGH/SNP histograms to the right of the main chr cartoons.
		if (CNVhistPlot == true) || (SNPhistPlot == true)
			width     = 0.020;
			height    = chr_height(chr);
			bottom    = chr_posY(chr);
			chr_CNVdata;
			histAll   = [];
			histAll2  = [];
			smoothed  = [];
			smoothed2 = [];

			if (SNPhistPlot == true)
				fprintf(['\n' num2str(chr_SNPdata{chr,1}) '\n']);
				fprintf(['\n' num2str(chr_SNPdata{chr,2}) '\n']);
				fprintf(['\n' num2str(chr_SNPdata{chr,3}) '\n']);
			end;

			for segment = 1:length(chrCopyNum{chr})
				subplot('Position',[(left+chr_width(chr)+0.005)+width*(segment-1) bottom width height]);

				if (CNVhistPlot == true)
					for i = round(1+length(CNVplot2{chr})*chr_breaks{chr}(segment)):round(length(CNVplot2{chr})*chr_breaks{chr}(segment+1))
						if (ploidy == 0) % no value; no scaling.
							y_ = CNVplot2{chr}(i);
						else
							if (Yscale_nearest_even_ploidy == true)
								if (abs(2-ploidy) < abs(4-ploidy)) % nearer to diploid
									y_ = CNVplot2{chr}(i)*ploidy/2*2;
								else
									y_ = CNVplot2{chr}(i)*ploidy/4*2;
								end;
							else
								y_ = CNVplot2{chr}(i)*ploidy/2*2;
							end;
						end;
						histAll{segment}(i) = y_*ploidyAdjust;
					end;
				elseif (SNPhistPlot == true)
					fprintf(['segment[' num2str(1+length(HETplot2{chr})*chr_breaks{chr}(segment)) '..' num2str(length(HETplot2{chr})*chr_breaks{chr}(segment+1)) ']\n']);
					for i = round(1+length(HETplot2{chr})*chr_breaks{chr}(segment)):round(length(HETplot2{chr})*chr_breaks{chr}(segment+1))
						colorMix =	colorHET   *   HETplot2{chr}(i)/TOTplot2{chr}(i) + ...
									colorOddHET*oddHETplot2{chr}(i)/TOTplot2{chr}(i) + ...
									colorHOM   *   HOMplot2{chr}(i)/TOTplot2{chr}(i);
						c_post =	colorMix   *   min(1,TOTplot2{chr}(i)) + ...
									colorNoData*(1-min(1,TOTplot2{chr}(i)));

						y_ = chr_SNPdata{chr,1}(i);
						histAll{segment}(i) = y_;
					end;
				end;

				histAll_ = histAll{segment};
				histAll_(isnan(histAll_)) = [];
				fprintf(['chr = ' num2str(chr)  '; min = ' num2str(min(histAll_)) ...
				         '; max = ' num2str(max(histAll_)) ...
				         '; median = ' num2str(median(histAll_)) '\n']);
				fprintf(['\tlength = ' num2str(length(histAll{segment})) '\n']);
				fprintf(['\tsum    = ' num2str(sum(histAll{segment})) '\n']);
				fprintf([num2str(histAll{segment}) '\n']);

				if (CNVhistPlot == true)
					threshold = 6;
				elseif (SNPhistPlot == true)
					threshold = 25;
				end;

				% make a histogram of CNV or SNP data, then smooth it for display.
				histAll{segment}(histAll{segment}<=0) = [];
				histAll{segment}(histAll{segment}>threshold) = threshold;
				histAll{segment}(length(histAll{segment})+1) = 0;   % endpoints added to ensure histogram bounds.
				histAll{segment}(length(histAll{segment})+1) = threshold;
				smoothed{segment} = smooth_gaussian(hist(histAll{segment},300),5,20);
				% make a smoothed version of just the endpoints used to ensure histogram bounds.
				histAll2{segment}(1) = 0;
				histAll2{segment}(2) = threshold;
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
			% subplot to show chr copy number value.
			width  = 0.020;
			height = chr_height(chr);
			bottom = chr_posY(chr);

			if (CNVhistPlot == true) || (SNPhistPlot == true)
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
				% chr_string = num2str(chrCopyNum{chr}(1));
				for i = 2:length(chrCopyNum{chr})
					chr_string = [chr_string ',' num2str(chrCopyNum{chr}(i))];
				end;
			end;
			text(0.1,0.5, chr_string,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',12);

			stringChrCNVs = [stringChrCNVs ';' chr_string];
		end;

		%% Linear figure draw section
		if (Linear_display == true)
			figure(Linear_fig);
			Linear_width = Linear_Chr_max_width*chr_size(chr)/Linear_genome_size;
			subplot('Position',[Linear_left Linear_base Linear_width Linear_height]);
			Linear_left = Linear_left + Linear_width + Linear_left_chr_gap;
			hold on;
			title(chr_label{chr},'Interpreter','none','FontSize',10);

 			% draw colorbars.
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
				if (c_(1) < 0); c_(1) = 0; end;
				if (c_(2) < 0); c_(2) = 0; end;
				if (c_(3) < 0); c_(3) = 0; end;
				if isnan(c_(1)); c_(1) = 1; end;
				if isnan(c_(2)); c_(2) = 1; end;
				if isnan(c_(3)); c_(3) = 1; end;
				if (blendColorBars == false)
					f = fill(x_,y_,c_);
				else
					f = fill(x_,y_,c_/2+c_prev/4+c_post/4);
				end;
				c_prev = c_;
				c_     = c_post;
				set(f,'linestyle','none');
			end;

			%% cgh plot section.
			c_ = [0 0 0];
			fprintf(['chr' num2str(chr) ':' num2str(length(CNVplot2{chr})) '\n']);
			for i = 1:length(CNVplot2{chr}); 
				x_ = [i i i-1 i-1];
				% ploidy_factor : is used to adjust default ploidy in figure generation.  Haploid default ploidy position is set to 1/4 of the
				% scale of 1..4 while diploid position is set to 1/2 of the scale.
				if (ploidy_default == 1)
					ploidy_factor = 4;
				elseif (ploidy_default == 2)
					ploidy_factor = 2;
				end;
				if (ploidy == 0) % no value; no scaling.
					y_ = [maxY/ploidy_factor CNVplot2{chr}(i)*maxY/2 CNVplot2{chr}(i)*maxY/2 maxY/ploidy_factor];
				else
					if (Yscale_nearest_even_ploidy == true)
						if (abs(2-ploidy) < abs(4-ploidy)) % nearer to diploid   
							y_ = [maxY/ploidy_factor CNVplot2{chr}(i)*ploidy/2*maxY/2 CNVplot2{chr}(i)*ploidy/2*maxY/2 maxY/ploidy_factor];
						else
							y_ = [maxY/ploidy_factor CNVplot2{chr}(i)*ploidy/4*maxY/2 CNVplot2{chr}(i)*ploidy/4*maxY/2 maxY/ploidy_factor];
						end;
					else
						y_ = [maxY/ploidy_factor CNVplot2{chr}(i)*ploidy/2*maxY/2 CNVplot2{chr}(i)*ploidy/2*maxY/2 maxY/ploidy_factor];
					end; 
				end;

				% makes a blackbar for each bin.
				f = fill(x_,y_,c_);
				set(f,'linestyle','none');
			end;
			x2 = chr_size(chr)*chr_length_scale_multiplier;
			plot([0; x2], [maxY/2; maxY/2],'color',[0 0 0]);  % 2n line.

			%% draw lines across plots for easier interpretation of CNV regions.
			if (ploidy_default == 2)
				line([0 x2], [maxY/4*3 maxY/4*3],'Color',[0.85 0.85 0.85]);
				line([0 x2], [maxY/4   maxY/4],  'Color',[0.85 0.85 0.85]);
			elseif (ploidy_default == 1)
				line([0 x2], [maxY/4*3 maxY/4*3],'Color',[0.85 0.85 0.85]);      
				line([0 x2], [maxY/2   maxY/2],  'Color',[0.85 0.85 0.85]);
			end;
			%% end cgh plot section.

			%show segmental anueploidy breakpoints.
			if (displayBREAKS == true)
				for segment = 2:length(chr_breaks{chr})-1
					bP = chr_breaks{chr}(segment)*length(CNVplot2{chr});
					c_ = [0 0 1];
					x_ = [bP bP bP-1 bP-1];
					y_ = [0 maxY maxY 0];
					f = fill(x_,y_,c_);   
					set(f,'linestyle','none');
				end;
			end;

			%show centromere.
			x1 = cen_start(chr)*chr_length_scale_multiplier;
			x2 = cen_end(chr)*chr_length_scale_multiplier;
			leftEnd  = 0.5*(5000/bases_per_bin);
			rightEnd = chr_size(chr)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);

			if (Centromere_format == 0)
				% standard chromosome cartoons in a way which will not cause segfaults when running via commandline.
				dx = cen_tel_Xindent;
				dy = cen_tel_Yindent;
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
				plot([leftEnd   leftEnd   leftEnd+dx   x1-dx   x1        x2        x2+dx   rightEnd-dx   rightEnd   rightEnd   rightEnd-dx   x2+dx   x2   x1   x1-dx   leftEnd+dx  leftEnd],...
				     [dy        maxY-dy   maxY         maxY    maxY-dy   maxY-dy   maxY    maxY          maxY-dy    dy         0             0       dy   dy   0       0           dy  ],...
				      'Color',[0 0 0]);
			end;
			%end show centromere.

			%show annotation locations
			if (show_annotations) && (length(annotations) > 0)
				plot([leftEnd rightEnd], [-1.5 -1.5],'color',[0 0 0]);
				hold on;
				annotation_location = (annotation_start+annotation_end)./2;
				for i = 1:length(annotation_location)
					if (annotation_chr(i) == chr)
						annotationloc = annotation_location(i)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
						plot(annotationloc,-1.5,'k:o','MarkerEdgeColor',annotation_edgecolor{i},'MarkerFaceColor',annotation_fillcolor{i},'MarkerSize',annotation_size(i));
						% plot(annotationloc,-1.5,'k:o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);
					end;
				end;
				hold off;
			end;
			%end show annotation locations.

			%% Final formatting stuff.
			xlim([0,chr_size(chr)*chr_length_scale_multiplier]);
			% modify y axis limits to show annotation locations if any are provided.
			if (length(annotations) > 0)
				ylim([-1.5,maxY]);
			else
				ylim([0,maxY]);
			end;
			set(gca,'TickLength',[(Linear_TickSize*chr_size(1)/chr_size(chr)) 0]); %ensures same tick size on all subfigs.
			set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
			set(gca,'XTickLabel',[]);
			if (chr == 1)
				ylabel(projectName2_child, 'Rotation', 0, 'HorizontalAlign', 'right', 'VerticalAlign', 'bottom','Interpreter','none','FontSize',5);
				set(gca,'YTick',[Linear_maxY/4 Linear_maxY/2 Linear_maxY/4*3 Linear_maxY]);
				if (Yscale_nearest_even_ploidy == true)
					if (abs(2-ploidy) < abs(4-ploidy)) % nearer to diploid
						set(gca,'YTickLabel',{'1','2','3','4'});
					else % nearer to tetraploid
						set(gca,'YTickLabel',{'2','4','6','8'});
					end;
				else
					set(gca,'YTickLabel',{'1','2','3','4'});
				end;
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

if (strcmp(projectName1_parent,projectName2_child) == 1)
	saveas(fig, [figureDir projectName1_parent '.normalized_CNV-SNP-map.1.eps'], 'epsc');
	saveas(fig, [figureDir projectName1_parent '.normalized_CNV-SNP-map.1.png'], 'png');
	set(Linear_fig,'PaperPosition',[0 0 8 0.62222222]*2);
	saveas(Linear_fig, [figureDir projectName1_parent '.normalized_CNV-SNP-map.2.eps'], 'epsc');
	saveas(Linear_fig, [figureDir projectName1_parent '.normalized_CNV-SNP-map.2.png'], 'png');
else
	saveas(fig, [figureDir projectName1_parent '->' projectName2_child '.normalized_CNV-LOH-map.1.eps'], 'epsc');
	saveas(fig, [figureDir projectName1_parent '->' projectName2_child '.normalized_CNV-LOH-map.1.png'], 'png');
	set(Linear_fig,'PaperPosition',[0 0 8 0.62222222]*2);
	saveas(fig, [figureDir projectName1_parent '->' projectName2_child '.normalized_CNV-LOH-map.2.eps'], 'epsc');
	saveas(fig, [figureDir projectName1_parent '->' projectName2_child '.normalized_CNV-LOH-map.2.png'], 'png');
end;

% Output chromosome copy number estimates.
textFileName = [figureDir projectName2_child '.CNV-map.3.txt'];
fprintf(['Text output of CNVs : "' textFileName '"\n']);
textFileID = fopen(textFileName,'w');
fprintf(textFileID,stringChrCNVs);
fclose(textFileID);

end
