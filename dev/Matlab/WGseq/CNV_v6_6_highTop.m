function [] = CNV_v6_6_highTop(main_dir,user,genomeUser,project,genome,ploidyEstimateString,ploidyBaseString, ...
                               CNV_verString,rDNA_verString,displayBREAKS, referenceCHR);

%% ========================================================================
Centromere_format_default   = 0;
Yscale_nearest_even_ploidy  = true;
HistPlot                    = true;
ChrNum                      = true;
show_annotations            = true;
analyze_rDNA                = true;
Linear_display              = true;
Low_quality_ploidy_estimate = true;


%%=========================================================================
% Load FASTA file name from 'reference.txt' file for project.
%--------------------------------------------------------------------------
Reference    = [main_dir 'users/' genomeUser '/genomes/' genome '/reference.txt'];
FASTA_string = strtrim(fileread(Reference));
[FastaPath,FastaName,FastaExt] = fileparts(FASTA_string);


%%=========================================================================
% Control variables.
%--------------------------------------------------------------------------
projectDir = [main_dir 'users/' user '/projects/' project '/'];
genomeDir  = [main_dir 'users/' genomeUser '/genomes/' genome '/'];

fprintf(['\n$$ projectDir : ' projectDir '\n']);
fprintf([  '$$ genomeDir  : ' genomeDir  '\n']);
fprintf([  '$$ genome     : ' genome     '\n']);
fprintf([  '$$ project    : ' project    '\n']);
[centromeres, chr_sizes, figure_details, annotations, ploidy_default] = Load_genome_information_1(projectDir,genomeDir, genome);
Aneuploidy = [];

num_chrs  = length(chr_sizes);

for i = 1:num_chrs
	chr_size(i)  = 0;
	cen_start(i) = 0;
	cen_end(i)   = 0;
end;
for i = 1:num_chrs
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
        fprintf(['\t[' num2str(annotations(i).chr) ':' annotations(i).type ':' num2str(annotations(i).start) ':' num2str(annotations(i).end) ':' annotations(i).fillcolor ':' ...
            annotations(i).edgecolor ':' num2str(annotations(i).size) ']\n']);
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
		chr_posY  (figure_details(i).chr) = figure_details(i).posY + 0.1 + 0.025*figure_details(i).chr;
		chr_width (figure_details(i).chr) = figure_details(i).width;
		chr_height(figure_details(i).chr) = figure_details(i).height;
		chr_in_use(figure_details(i).chr) = str2num(figure_details(i).useChr);
	end;
end;


%%=========================================================================
%%= No further control variables below. ===================================
%%=========================================================================

% Sanitize user input of euploid state base for species.
ploidyBase = round(str2num(ploidyBaseString));
if (ploidyBase > 4);   ploidyBase = 4;   end;
if (ploidyBase < 1);   ploidyBase = 1;   end; 
fprintf(['\nEuploid base = "' num2str(ploidyBase) '"\n']);

% basic plot parameters not defined per genome.
TickSize         = 0;  % -0.005;  %negative for outside, percentage of longest chr figure.
bases_per_bin    = max(chr_size)/700;
maxY             = ploidyBase*2;
cen_tel_Xindent  = 5;
cen_tel_Yindent  = maxY/5;

fprintf(['\nGenerating horizontal CNV highTop figure from ''' project ''' sequence data.\n']);


%%================================================================================================
% Load corrected CGH data for display.
%-------------------------------------------------------------------------------------------------
fprintf('\nCommon_CNV data file found, loading.\n');
load([projectDir 'Common_CNV.mat']);   %% 'CNVplot2','genome_CNV'.


%% -----------------------------------------------------------------------------------------
% Calculate chromosome copy number from ploidy estimatre.
%-------------------------------------------------------------------------------------------
ploidy = str2num(ploidyEstimateString);
[chr_breaks, chrCopyNum, ploidyAdjust] = FindChrSizes_4(Aneuploidy,CNVplot2,ploidy,num_chrs,chr_in_use);
fprintf('\n');
largestChr = find(chr_width == max(chr_width));


%% -----------------------------------------------------------------------------------------
% Setup for main figure generation.
%------------------------------------------------------------------------------------------
% threshold for full color saturation in SNP/LOH figure.
% synced to bases_per_bin as below, or defaulted to 50.
full_data_threshold = floor(bases_per_bin/100);

Standard_fig = figure();
set(gcf, 'Position', [0 70 1024 600]);
largestChr = find(chr_width == max(chr_width));


%% -----------------------------------------------------------------------------------------
% Setup for linear-view figure generation.
%-------------------------------------------------------------------------------------------
if (Linear_display == true)
	Linear_fig           = figure();
	Linear_genome_size   = sum(chr_size);
	Linear_Chr_max_width = 0.91;               % width for all chromosomes across figure.  1.00 - leftMargin - rightMargin - subfigure gaps.
	Linear_left_start    = 0.02;               % left margin (also right margin).  (formerly 0.01)
	Linear_left_chr_gap  = 0.07/(num_chrs-1);  % gaps between chr subfigures.
	Linear_height        = 0.6;                % height of highPlot subfigures.
	Linear_base          = 0.1;                % base position of subfigures in figure.
	Linear_TickSize      = -0.01;              % negative for outside, percentage of longest chr figure.
	maxY                 = ploidyBase*2;       % maximum y-axis of chromosome cartoons.
	maxY_highTop         = ploidyBase*2*3;     % maximum y-axis of region above chromosome cartoons.
	Linear_left          = Linear_left_start;  % used to track left end of current chromosome.

	axisLabelPosition_horiz = -50000/bases_per_bin;
	axisLabelPosition_horiz = 0.01125;
end;

axisLabelPosition_vert = -50000/bases_per_bin;
axisLabelPosition_vert = 0.01125;

maxY_highTop           = ploidyBase*2*3;

%% Initialize copy numbers string.
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
medianCNV = median(CNVdata_all)
for chr = 1:num_chrs
    if (chr_in_use(chr) == 1)
        CNVplot2{chr} = CNVplot2{chr}/medianCNV;
    end;
end;


%% -----------------------------------------------------------------------------------------
% Make figures
%-------------------------------------------------------------------------------------------
first_chr = true;
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		%% Standard figure draw section.
		figure(Standard_fig);
		left   = chr_posX(chr);
		bottom = chr_posY(chr);
		width  = chr_width(chr);
		height = chr_height(chr);
		subplot('Position',[left bottom width height]);
		fprintf(['figposition = [' num2str(left) ' | ' num2str(bottom) ' | ' num2str(width) ' | ' num2str(height) ']\t']);
		hold on;

		%% cgh plot section.
		c_ = [0 0 0];
		fprintf(['chr' num2str(chr) ':' num2str(length(CNVplot2{chr})) '\n']);
		for i = 1:length(CNVplot2{chr});
			x_ = [i i i-1 i-1];
			if (CNVplot2{chr}(i) == 0)
				CNVhistValue = 1;
			else
				CNVhistValue = CNVplot2{chr}(i);
			end;

			% The CNV-histogram values were normalized to a median value of 1.
			% The ratio of 'ploidy' to 'ploidyBase' determines where the data is displayed relative to the
			% median line.
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

		x2 = chr_size(chr)/bases_per_bin;
		plot([0; x2], [maxY/2; maxY/2],'color',[0 0 0]);  % 2n line.

		%% draw lines across plots for easier interpretation of CNV regions.
		switch ploidyBase
		case 1
			% above chr bounds.
			for lineNum = 3:6
				line([0 x2], [maxY/2*lineNum  maxY/2*lineNum ],'Color',[0.85 0.85 0.85]);
			end;
		case 2
			% inside chr bounds.
			line([0 x2], [maxY/4*1 maxY/4*1],'Color',[0.85 0.85 0.85]);
			line([0 x2], [maxY/4*3 maxY/4*3],'Color',[0.85 0.85 0.85]);
			% above chr bounds.
			for lineNum = 5:12
				line([0 x2], [maxY/4*lineNum  maxY/4*lineNum ],'Color',[0.85 0.85 0.85]);
			end;
		case 3
			% inside chr bounds.
			line([0 x2], [maxY/6*1 maxY/6*1],'Color',[0.85 0.85 0.85]);
			line([0 x2], [maxY/6*2 maxY/6*2],'Color',[0.85 0.85 0.85]);
			line([0 x2], [maxY/6*4 maxY/6*4],'Color',[0.85 0.85 0.85]);
			line([0 x2], [maxY/6*5 maxY/6*5],'Color',[0.85 0.85 0.85]);
			% above chr bounds.
			for lineNum = 7:18
				line([0 x2], [maxY/6*lineNum  maxY/6*lineNum ],'Color',[0.85 0.85 0.85]);
			end;
		case 4
			% inside chr bounds.
			line([0 x2], [maxY/8*1 maxY/8*1],'Color',[0.85 0.85 0.85]);
			line([0 x2], [maxY/8*2 maxY/8*2],'Color',[0.85 0.85 0.85]);
			line([0 x2], [maxY/8*3 maxY/8*3],'Color',[0.85 0.85 0.85]);
			line([0 x2], [maxY/8*5 maxY/8*5],'Color',[0.85 0.85 0.85]);
			line([0 x2], [maxY/8*6 maxY/8*6],'Color',[0.85 0.85 0.85]);
			line([0 x2], [maxY/8*7 maxY/8*7],'Color',[0.85 0.85 0.85]);
			% above chr bounds.
			for lineNum = 9:24
				line([0 x2], [maxY/8*lineNum  maxY/8*lineNum ],'Color',[0.85 0.85 0.85]);
			end;
		end;
		%% end cgh plot section.

		%axes labels etc.
		hold off;
		% limit x-axis to range of chromosome.
		xlim([0,chr_size(chr)/bases_per_bin]);

		% modify y axis limits to show annotation locations if any are provided.
		if (length(annotations) > 0)
			ylim([-maxY/10*1.5,maxY_highTop]);
		else
			ylim([0,maxY_highTop]);
		end;

		set(gca,'YTick',[]);
		set(gca,'TickLength',[(TickSize*chr_size(largestChr)/chr_size(chr)) 0]); %ensures same tick size on all subfigs.

		% ylabel(chr_label{chr}, 'Rotation', 90, 'HorizontalAlign', 'center', 'VerticalAlign', 'bottom');
		text(-50000/5000/2*3, maxY/2,     chr_label{chr}, 'Rotation',90, 'HorizontalAlignment','center', 'VerticalAlign','bottom', 'Fontsize',20);

		set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
		set(gca,'XTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8','2.0','2.2','2.4','2.6','2.8','3.0','3.2'});

		% This section sets the Y-axis labelling.
		switch ploidyBase
			case 1
				set(gca,'YTick',[0 maxY*1/2 maxY*2/2 maxY*3/2 maxY*4/2 maxY*5/2 maxY*6/2]);
				set(gca,'YTickLabel',{'','','','','','',''});
				text(axisLabelPosition_vert, maxY*1/2,    '1' ,'HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*2/2,    '2' ,'HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*3/2,    '3' ,'HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*4/2,    '4' ,'HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*5/2,    '5' ,'HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*6/2,    '6' ,'HorizontalAlignment','right','Fontsize',10);
			case 2
				set(gca,'YTick',[0 maxY*1/4 maxY*2/4 maxY*3/4 maxY*4/4 maxY*5/4 maxY*6/4 maxY*7/4 maxY*8/4 ...
				                 maxY*9/4 maxY*10/4 maxY*11/4 maxY*12/4]);
				set(gca,'YTickLabel',{'','','','','','','','','','','','',''});
				text(axisLabelPosition_vert, maxY*1/4,    '1' ,'HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*2/4,    '2' ,'HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*3/4,    '3' ,'HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*4/4,    '4' ,'HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*5/4,    '5' ,'HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*6/4,    '6' ,'HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*7/4,    '7' ,'HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*8/4,    '8' ,'HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*9/4,    '9' ,'HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*10/4,   '10','HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*11/4,   '11','HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*12/4,   '12','HorizontalAlignment','right','Fontsize',10);
			case 3
				set(gca,'YTick',[0 maxY*1/6 maxY*2/6 maxY*3/6 maxY*4/6 maxY*5/6 maxY*6/6 maxY*7/6 maxY*8/6 ...
				                 maxY*9/6 maxY*10/6 maxY*11/6 maxY*12/6 maxY*13/6 maxY*14/6 maxY*15/6 maxY*16/6 ...
				                 maxY*17/6 maxY*18/6]);
				set(gca,'YTickLabel',{'','','','','','','','','','','','',','','','','','','});
				text(axisLabelPosition_vert, maxY*3/6,    '3' ,'HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*6/6,    '6' ,'HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*9/6,    '9' ,'HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*12/6,   '12','HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*15/6,   '15','HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*18/6,   '18','HorizontalAlignment','right','Fontsize',10);
			case 4
				set(gca,'YTick',[0 maxY*1/8 maxY*2/8 maxY*3/8 maxY*4/8 maxY*5/8 maxY*6/8 maxY*7/8 maxY*8/8 ...
				                 maxY*9/8 maxY*10/8 maxY*11/8 maxY*12/8 maxY*13/8 maxY*14/8 maxY*15/8 maxY*16/8 ...
				                 maxY*17/8 maxY*18/8 maxY*19/8 maxY*20/8 maxY*21/8 maxY*22/8 maxY*23/8 maxY*24/8]);
				set(gca,'YTickLabel',{'','','','','','','','','','','','','','','','','','','','','','','','',''});
				text(axisLabelPosition_vert, maxY*2/8,    '2' ,'HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*4/8,    '4' ,'HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*6/8,    '6' ,'HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*8/8,    '8' ,'HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*10/8,   '10','HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*12/8,   '12','HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*14/8,   '14','HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*16/8,   '16','HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*18/8,   '18','HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*20/8,   '20','HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*22/8,   '22','HorizontalAlignment','right','Fontsize',10);
				text(axisLabelPosition_vert, maxY*24/8,   '24','HorizontalAlignment','right','Fontsize',10);
		end;

		set(gca,'FontSize',6);
		if (chr == find(chr_posY == max(chr_posY)))
			title([ project ' CNV map'],'Interpreter','none','FontSize',12);
		end;

		hold on;
		% end axes labels etc.

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
		if (chr_size(chr) < 100000)
			Centromere_format = 1;
		else
			Centromere_format = Centromere_format_default;
		end;
		x1 = cen_start(chr)/bases_per_bin;
		x2 = cen_end(chr)/bases_per_bin;
		leftEnd  = 0.5*(5000/bases_per_bin);
		rightEnd = chr_size(chr)/bases_per_bin-0.5*(5000/bases_per_bin);
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

			% draw outlines of chromosome cartoon.   (drawn after horizontal lines to that cartoon edges are
			% not interrupted by horiz lines.

			plot([leftEnd   leftEnd   leftEnd+dx   x1-dx   x1        x2        x2+dx   rightEnd-dx   rightEnd   rightEnd   rightEnd-dx   x2+dx   x2   x1   x1-dx   leftEnd+dx   leftEnd],...
				 [dy        maxY-dy   maxY         maxY    maxY-dy   maxY-dy   maxY    maxY          maxY-dy    dy         0             0       dy   dy   0       0            dy     ],...
				'Color',[0 0 0]);
		elseif (Centromere_format == 1)
			leftEnd  = 0;
			rightEnd = chr_size(chr)/bases_per_bin;

			% Minimal outline for examining very small sequence regions, such as C.albicans MTL locus.
			plot([leftEnd   leftEnd   rightEnd   rightEnd   leftEnd], [0   maxY   maxY   0   0], 'Color',[0 0 0]);
		end;
		%end show centromere.

		%show annotation locations
		if (show_annotations) && (length(annotations) > 0)
			plot([leftEnd rightEnd], [-maxY/10*1.5 -maxY/10*1.5],'color',[0 0 0]);
			hold on;
			annotation_location = (annotation_start+annotation_end)./2;
			for i = 1:length(annotation_location)
				if (annotation_chr(i) == chr)
					annotationLoc   = annotation_location(i)/bases_per_bin-0.5*(5000/bases_per_bin);
					annotationStart = annotation_start(i)/bases_per_bin-0.5*(5000/bases_per_bin);
					annotationEnd   = annotation_end(i)/bases_per_bin-0.5*(5000/bases_per_bin);
					if (strcmp(annotation_type{i},'dot') == 1)
						plot(annotationLoc,-maxY/10*1.5,'k:o','MarkerEdgeColor',annotation_edgecolor{i}, ...
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
		%end show annotation locations.

		%% Linear figure draw section.
		if (Linear_display == true)
			figure(Linear_fig);
			Linear_width = Linear_Chr_max_width*chr_size(chr)/Linear_genome_size;
			subplot('Position',[Linear_left Linear_base Linear_width Linear_height]);
			Linear_left = Linear_left + Linear_width + Linear_left_chr_gap;
			hold on;
			title(chr_label{chr},'Interpreter','none','FontSize',20);

			%% cgh plot section.
			c_ = [0 0 0];
			fprintf(['chr' num2str(chr) ':' num2str(length(CNVplot2{chr})) '\n']);
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
			x2 = chr_size(chr)/bases_per_bin;
			plot([0; x2], [maxY/2; maxY/2],'color',[0 0 0]);  % 2n line.

			%% draw lines across plots for easier interpretation of CNV regions.
			switch ploidyBase
				case 1
					% above chr bounds.
					for lineNum = 3:6
						line([0 x2], [maxY/2*lineNum  maxY/2*lineNum ],'Color',[0.85 0.85 0.85]);
					end;
				case 2
					% inside chr bounds.
					line([0 x2], [maxY/4*1  maxY/4*1 ],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/4*3  maxY/4*3 ],'Color',[0.85 0.85 0.85]);
					% above chr bounds.
					for lineNum = 5:12
						line([0 x2], [maxY/4*lineNum  maxY/4*lineNum ],'Color',[0.85 0.85 0.85]);
					end;
				case 3
					% inside chr bounds.
					line([0 x2], [maxY/6*1 maxY/6*1],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/6*2 maxY/6*2],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/6*4 maxY/6*4],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/6*5 maxY/6*5],'Color',[0.85 0.85 0.85]);
					% above chr bounds.
					for lineNum = 7:18   
						line([0 x2], [maxY/6*lineNum  maxY/6*lineNum ],'Color',[0.85 0.85 0.85]);
					end;
				case 4
					% inside chr bounds.
					line([0 x2], [maxY/8*1 maxY/8*1],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/8*2 maxY/8*2],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/8*3 maxY/8*3],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/8*5 maxY/8*5],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/8*6 maxY/8*6],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/8*7 maxY/8*7],'Color',[0.85 0.85 0.85]);
					% above chr bounds.
					for lineNum = 9:24   
						line([0 x2], [maxY/8*lineNum  maxY/8*lineNum ],'Color',[0.85 0.85 0.85]);
					end;
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
			if (chr_size(chr) < 100000)
				Centromere_format = 1;
			else
				Centromere_format = Centromere_format_default;
			end;
			x1 = cen_start(chr)/bases_per_bin;
			x2 = cen_end(chr)/bases_per_bin;
			leftEnd  = 0.5*(5000/bases_per_bin);
			rightEnd = chr_size(chr)/bases_per_bin-0.5*(5000/bases_per_bin);
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
			elseif (Centromere_format == 1)
				leftEnd  = 0;
				rightEnd = chr_size(chr)/bases_per_bin;

				% Minimal outline for examining very small sequence regions, such as C.albicans MTL locus.
				plot([leftEnd   leftEnd   rightEnd   rightEnd   leftEnd], [0   maxY   maxY   0   0], 'Color',[0 0 0]);
			end;
			%end show centromere.

			%show annotation locations
			if (show_annotations) && (length(annotations) > 0)
				plot([leftEnd rightEnd], [-maxY/10*1.5 -maxY/10*1.5],'color',[0 0 0]);
				hold on;
				annotation_location = (annotation_start+annotation_end)./2;
				for i = 1:length(annotation_location)
					if (annotation_chr(i) == chr)
						annotationLoc   = annotation_location(i)/bases_per_bin-0.5*(5000/bases_per_bin);
						annotationStart = annotation_start(i)/bases_per_bin-0.5*(5000/bases_per_bin);
						annotationEnd   = annotation_end(i)/bases_per_bin-0.5*(5000/bases_per_bin);
						if (strcmp(annotation_type{i},'dot') == 1)
							plot(annotationLoc,-maxY/10*1.5,'k:o','MarkerEdgeColor',annotation_edgecolor{i}, ...
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
			%end show annotation locations.

			%% Final formatting stuff.
			xlim([0,chr_size(chr)/bases_per_bin]);

			% modify y axis limits to show annotation locations if any are provided.
			if (length(annotations) > 0)
				ylim([-maxY/10*1.5,maxY_highTop]);
			else
				ylim([0,maxY_highTop]);
			end;

			set(gca,'TickLength',[(Linear_TickSize*chr_size(1)/chr_size(chr)) 0]); %ensures same tick size on all subfigs.
			set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
			set(gca,'XTickLabel',[]);
			%if (first_chr)
				% This section sets the Y-axis labelling.
				switch ploidyBase
				case 1
					set(gca,'YTick',[0 maxY*1/2 maxY*2/2 maxY*3/2 maxY*4/2 maxY*5/2 maxY*6/2]);
					set(gca,'YTickLabel',{'','','','','','',''});
					text(axisLabelPosition_horiz, maxY*1/2,    '1' ,'HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*2/2,    '2' ,'HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*3/2,    '3' ,'HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*4/2,    '4' ,'HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*5/2,    '5' ,'HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*6/2,    '6' ,'HorizontalAlignment','right','Fontsize',10);
				case 2
					set(gca,'YTick',[0 maxY*1/4 maxY*2/4 maxY*3/4 maxY*4/4 maxY*5/4 maxY*6/4 maxY*7/4 maxY*8/4 ...
					                 maxY*9/4 maxY*10/4 maxY*11/4 maxY*12/4]);
					set(gca,'YTickLabel',{'','','','','','','','','','','','',''});
					text(axisLabelPosition_horiz, maxY*1/4,    '1' ,'HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*2/4,    '2' ,'HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*3/4,    '3' ,'HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*4/4,    '4' ,'HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*5/4,    '5' ,'HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*6/4,    '6' ,'HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*7/4,    '7' ,'HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*8/4,    '8' ,'HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*9/4,    '9' ,'HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*10/4,   '10','HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*11/4,   '11','HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*12/4,   '12','HorizontalAlignment','right','Fontsize',10);
				case 3
					set(gca,'YTick',[0 maxY*1/6 maxY*2/6 maxY*3/6 maxY*4/6 maxY*5/6 maxY*6/6 maxY*7/6 maxY*8/6 ...
					                 maxY*9/6 maxY*10/6 maxY*11/6 maxY*12/6 maxY*13/6 maxY*14/6 maxY*15/6 maxY*16/6 ...
					                 maxY*17/6 maxY*18/6]);
					set(gca,'YTickLabel',{'','','','','','','','','','','','',','','','','','','});
					text(axisLabelPosition_horiz, maxY*3/6,    '3' ,'HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*6/6,    '6' ,'HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*9/6,    '9' ,'HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*12/6,   '12','HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*15/6,   '15','HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*18/6,   '18','HorizontalAlignment','right','Fontsize',10);
				case 4
					set(gca,'YTick',[0 maxY*1/8 maxY*2/8 maxY*3/8 maxY*4/8 maxY*5/8 maxY*6/8 maxY*7/8 maxY*8/8 ...
					                 maxY*9/8 maxY*10/8 maxY*11/8 maxY*12/8 maxY*13/8 maxY*14/8 maxY*15/8 maxY*16/8 ...
					                 maxY*17/8 maxY*18/8 maxY*19/8 maxY*20/8 maxY*21/8 maxY*22/8 maxY*23/8 maxY*24/8]);
					set(gca,'YTickLabel',{'','','','','','','','','','','','','','','','','','','','','','','','',''});
					text(axisLabelPosition_horiz, maxY*2/8,    '2' ,'HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*4/8,    '4' ,'HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*6/8,    '6' ,'HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*8/8,    '8' ,'HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*10/8,   '10','HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*12/8,   '12','HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*14/8,   '14','HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*16/8,   '16','HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*18/8,   '18','HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*20/8,   '20','HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*22/8,   '22','HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_horiz, maxY*24/8,   '24','HorizontalAlignment','right','Fontsize',10);
				end;
			%else
			%	set(gca,'YTick',[]);
			%	set(gca,'YTickLabel',[]);
			%end;
			set(gca,'FontSize',12);
			%end final reformatting.

			% shift back to main figure generation.
			figure(Standard_fig);
			hold on;
	
			first_chr = false;
		end;
	end; 
end;
	

% Save primary genome figure.
set(Standard_fig,'PaperPosition',[0 0 8 6*3]*2);
saveas(Standard_fig, [projectDir 'fig.CNV-map.highTop.1.eps'], 'epsc');
saveas(Standard_fig, [projectDir 'fig.CNV-map.highTop.1.png'], 'png');
delete(Standard_fig);

% Save horizontal aligned genome figure.
set(Linear_fig,'PaperPosition',[0 0 8 0.62222222*3]*2);
saveas(Linear_fig,   [projectDir 'fig.CNV-map.highTop.2.eps'], 'epsc');
saveas(Linear_fig,   [projectDir 'fig.CNV-map.highTop.2.png'], 'png');
delete(Linear_fig);


end
