function [] = CNV_SNP_hapmap_v5_RADseq(main_dir,user,genomeUser,project,parent,hapmap,genome,ploidyEstimateString,ploidyBaseString, ...
                                       SNP_verString,LOH_verString,CNV_verString,displayBREAKS);
addpath('../');


fprintf('\n');
fprintf('#########################################\n');
fprintf('## CNV_SNP_hapmap_v5_RADseq_RedGreen.m ##\n');
fprintf('#########################################\n');


%% =========================================================================================
% Load workspace variables saved in "CNV_SNP_hapmap_v4.m"
%-------------------------------------------------------------------------------------------
projectDir  = [main_dir 'users/' user '/projects/' project '/'];
load([projectDir 'CNV_SNP_hapmap_v5_RADseq.workspace_variables.mat']);

if ((useHapmap) || (useParent))
	fprintf(['\n##\n## Hapmap in use, so "allelic_ratios_ddRADseq_D.m" is being processed.\n##\n']);


	%%================================================================================================
	% Load SNP/LOH data.
	%-------------------------------------------------------------------------------------------------
	LOH_file = [projectDir 'SNP_' SNP_verString '.reduced_RedGreen.mat'];
	load(LOH_file);
	% 'chr_SNPdata','new_bases_per_bin','chr_SNPdata_colorsC', 'chr_SNPdata_colorsP'.


	%% =========================================================================================
	% Setup for main figure generation.
	%-------------------------------------------------------------------------------------------
	Main_fig = figure(1);
	set(gcf, 'Position', [0 70 1024 600]);
	largestChr = find(chr_width == max(chr_width));
	largestChr = largestChr(1);


	%% =========================================================================================
	% Setup for linear-view figure generation.
	%-------------------------------------------------------------------------------------------
	if (Linear_display == true)
		Linear_fig              = figure(2);
		Linear_genome_size      = sum(chr_size);
		Linear_Chr_max_width    = 0.91;               % width for all chromosomes across figure.  1.00 - leftMargin - rightMargin - subfigure gaps.
		Linear_left_start       = 0.02;               % left margin (also right margin).
		Linear_left_chr_gap     = 0.07/(num_chrs-1);  % gaps between chr subfigures.
		Linear_height           = 0.6;
		Linear_base             = 0.1;
		Linear_TickSize         = -0.01;  %negative for outside, percentage of longest chr figure.
		maxY                    = ploidyBase*2;
		Linear_left             = Linear_left_start;
		axisLabelPosition_horiz = 0.01125;
	end;
	axisLabelPosition_vert = 0.01125;


	%% =========================================================================================
	% Make figures
	%-------------------------------------------------------------------------------------------
	first_chr = true;
	for chr = 1:num_chrs
		if (chr_in_use(chr) == 1)
			fprintf(['-------------------------------- drawing subfigures for [chr' num2str(chr) '] --------------------------------\n']);
			figure(Main_fig);
			% make standard chr cartoons.
			left   = chr_posX(chr);
			bottom = chr_posY(chr);
			width  = chr_width(chr);
			height = chr_height(chr);
			subplot('Position',[left bottom width height]);
			fprintf(['\tfigposition = [' num2str(left) ' | ' num2str(bottom) ' | ' num2str(width) ' | ' num2str(height) ']\n']);
			hold on;

			% standard : draw colorbars.
			if (useHapmap)
				for i = 1:ceil(chr_size(chr)/new_bases_per_bin)
					colorR   = chr_SNPdata_colorsC{chr,1}(i);
					colorG   = chr_SNPdata_colorsC{chr,2}(i);
					colorB   = chr_SNPdata_colorsC{chr,3}(i);
					plot([i i], [0 maxY],'Color',[colorR colorG colorB]);
				end;
			elseif (useParent)
				for i = 1:ceil(chr_size(chr)/new_bases_per_bin)
					datumY_C = chr_SNPdata{chr,2}(i)*maxY;
					datumY_P = chr_SNPdata{chr,4}(i)*maxY;
					plot([i/2 i/2], [maxY datumY_C     ],'Color',[1.0 0.0 0.0]);
					plot([i/2 i/2], [0    maxY-datumY_P],'Color',[1/3 1/3 1/3]);
				end;
			end;
			% end standard : draw color bars.


			%% standard : cgh plot section.
			c_ = [0 0 0];
			fprintf(['\nmain-plot : chr' num2str(chr) ':' num2str(length(CNVplot2{chr})) '\n']);
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
			% standard : end cgh plot section.


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
			% standard : end draw lines across plots for easier interpretation of CNV regions.


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
			set(gca,'YTickLabel',[]);
			set(gca,'TickLength',[(TickSize*chr_size(largestChr)/chr_size(chr)) 0]); %ensures same tick size on all subfigs.
			text(-50000/5000/2*3, maxY/2,     chr_label{chr}, 'Rotation',90, 'HorizontalAlignment','center', 'VerticalAlign','bottom', 'Fontsize',20);
			set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
			set(gca,'XTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8','2.0','2.2','2.4','2.6','2.8','3.0','3.2'});
			% This section sets the Y-axis labelling.
			switch ploidyBase
				case 1
					text(axisLabelPosition_vert, maxY/2,   '1','HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_vert, maxY,     '2','HorizontalAlignment','right','Fontsize',10);
				case 2
					text(axisLabelPosition_vert, maxY/4,   '1','HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_vert, maxY/2,   '2','HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_vert, maxY/4*3, '3','HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_vert, maxY,     '4','HorizontalAlignment','right','Fontsize',10);
				case 3
					text(axisLabelPosition_vert, maxY/2,   '3','HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_vert, maxY,     '6','HorizontalAlignment','right','Fontsize',10);
				case 4
					text(axisLabelPosition_vert, maxY/4,   '2','HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_vert, maxY/2,   '4','HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_vert, maxY/4*3, '6','HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_vert, maxY,     '8','HorizontalAlignment','right','Fontsize',10);
			end;
			set(gca,'FontSize',12);
			if (chr == find(chr_posY == max(chr_posY)))
				title([ project ' CNV map'],'Interpreter','none','FontSize',24);
			end;
			hold on;
			% standard : end axes labels etc.


			% standard : add red tick-marks beneath 
			if ((displayBREAKS == true) && (show_annotations == true))
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
					histogram_end                                    = 15;             % end point in copy numbers for the histogram, this should be way outside the expected range.
					histAll{segment}(histAll{segment}<=0)            = [];
					histAll{segment}(length(histAll{segment})+1)     = 0;              % endpoints added to ensure histogram bounds.
					histAll{segment}(length(histAll{segment})+1)     = histogram_end;
					histAll{segment}(histAll{segment}<0)             = [];             % crop off any copy data outside the range.
					histAll{segment}(histAll{segment}>histogram_end) = [];
					smoothed{segment}                                = smooth_gaussian(hist(histAll{segment},histogram_end*20),2,10);

					% make a smoothed version of just the endpoints used to ensure histogram bounds.
					histAll2{segment}(1)                             = 0;
					histAll2{segment}(2)                             = histogram_end;
					smoothed2{segment}                               = smooth_gaussian(hist(histAll2{segment},histogram_end*20),2,10);

					% subtract the smoothed endpoints from the histogram to remove the influence of the added endpoints.
					smoothed{segment}                                = (smoothed{segment}-smoothed2{segment});
					smoothed{segment}                                = smoothed{segment}/max(smoothed{segment});

					% draw lines to mark whole copy number changes.
					plot([0;       0      ],[0; 1],'color',[0.00 0.00 0.00]);
					hold on;
					for i = 1:15
						plot([20*i;  20*i],[0; 1],'color',[0.75 0.75 0.75]);
					end;

					% draw histogram.
					area(smoothed{segment},'FaceColor',[0 0 0]);

					% Draw red ticks between histplot segments
					if (displayBREAKS == true) && (show_annotations == true)
						if (segment > 1)
							plot([-maxY*20/10*1.5 0],[0 0],  'Color',[1 0 0],'LineWidth',2);
						end;
					end;

					% Flip subfigure around the origin.
					view(-90,90);
					set(gca,'YDir','Reverse');

					% ensure subplot axes are consistent with main chr plots.
					hold off;
					axis off;
					set(gca,'YTick',[]);    set(gca,'XTick',[]);
					ylim([0,1]);            xlim([0,maxY*20]);
					if (show_annotations == true)
						xlim([-maxY*20/10*1.5,maxY*20]);
					else
						xlim([0,maxY*20]);
					end;
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
				if (length(chrCopyNum{chr}) > 0)
					if (length(chrCopyNum{chr}) == 1)
						chr_string = num2str(chrCopyNum{chr}(1));
					else
						chr_string = num2str(chrCopyNum{chr}(1));
						for i = 2:length(chrCopyNum{chr})
							chr_string = [chr_string ',' num2str(chrCopyNum{chr}(i))];
						end;
					end;
					text(0.1,0.5, chr_string,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',24);
				end;
			end;
			% standard : end copy number at right.


			%% =========================================================================================
			% Draw angleplots to left of main chromosome cartoons.
			%-------------------------------------------------------------------------------------------
			apply_phasing = false;
			angle_plot_subfigures;


%%%%%%%%%%%%%%%%%%%%%%% END standard draw section.


			%% Linear figure draw section
			if (Linear_display == true)
				figure(Linear_fig);
				Linear_width = Linear_Chr_max_width*chr_size(chr)/Linear_genome_size;
				subplot('Position',[Linear_left Linear_base Linear_width Linear_height]);
				Linear_left = Linear_left + Linear_width + Linear_left_chr_gap;
				hold on;
				title(chr_label{chr},'Interpreter','none','FontSize',20);

				% linear : draw colorbars.
				if (useHapmap)
					for i = 1:ceil(chr_size(chr)/new_bases_per_bin)
						colorR   = chr_SNPdata_colorsC{chr,1}(i);
						colorG   = chr_SNPdata_colorsC{chr,2}(i);
						colorB   = chr_SNPdata_colorsC{chr,3}(i);
						if (colorR < 1) || (colorG < 1) || (colorB < 1)
							plot([i i], [0 maxY],'Color',[colorR colorG colorB]);
						end;
					end;
				elseif (useParent)
					for i = 1:ceil(chr_size(chr)/new_bases_per_bin)
						datumY_C = chr_SNPdata{chr,2}(i)*maxY;
						datumY_P = chr_SNPdata{chr,4}(i)*maxY;
						plot([i/2 i/2], [maxY datumY_C     ],'Color',[1.0 0.0 0.0]);
						plot([i/2 i/2], [0    maxY-datumY_P],'Color',[1/3 1/3 1/3]);
					end;
				else
					for i = 1:ceil(chr_size(chr)/new_bases_per_bin)
						datumY_C = chr_SNPdata{chr,2}(i)*maxY;
						datumY_P = chr_SNPdata{chr,4}(i)*maxY;
						plot([i/2 i/2], [maxY datumY_P     ],'Color',[1/3 1/3 1/3]);
						plot([i/2 i/2], [0    maxY-datumY_P],'Color',[1/3 1/3 1/3]);
					end;
				end;
				% end linear : draw colorbars.

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
				if (Linear_displayBREAKS == true) && (show_annotations == true)
					chr_length = ceil(chr_size(chr)/bases_per_bin);
                        	        for segment = 2:length(chr_breaks{chr})-1
                                	        bP = chr_breaks{chr}(segment)*chr_length;
                                        	plot([bP bP], [(-maxY/10*2.5) 0],  'Color',[1 0 0],'LineWidth',2);
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
				set(gca,'YTick',[]);
				set(gca,'YTickLabel',[]);
				set(gca,'TickLength',[(Linear_TickSize*chr_size(largestChr)/chr_size(chr)) 0]); %ensures same tick size on all subfigs.
				set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
				set(gca,'XTickLabel',[]);
				if (first_chr)
					% This section sets the Y-axis labelling.
					switch ploidyBase
						case 1
							text(axisLabelPosition_horiz, maxY/2,   '1','HorizontalAlignment','right','Fontsize',10);
							text(axisLabelPosition_horiz, maxY,     '2','HorizontalAlignment','right','Fontsize',10);
						case 2
							text(axisLabelPosition_horiz, maxY/4,   '1','HorizontalAlignment','right','Fontsize',10);
							text(axisLabelPosition_horiz, maxY/2,   '2','HorizontalAlignment','right','Fontsize',10);
							text(axisLabelPosition_horiz, maxY/4*3, '3','HorizontalAlignment','right','Fontsize',10);
							text(axisLabelPosition_horiz, maxY,     '4','HorizontalAlignment','right','Fontsize',10);
						case 3
							text(axisLabelPosition_horiz, maxY/2,   '3','HorizontalAlignment','right','Fontsize',10);
							text(axisLabelPosition_horiz, maxY,     '6','HorizontalAlignment','right','Fontsize',10);
						case 4
							text(axisLabelPosition_horiz, maxY/4,   '2','HorizontalAlignment','right','Fontsize',10);
							text(axisLabelPosition_horiz, maxY/2,   '4','HorizontalAlignment','right','Fontsize',10);
							text(axisLabelPosition_horiz, maxY/4*3, '6','HorizontalAlignment','right','Fontsize',10);
							text(axisLabelPosition_horiz, maxY,     '8','HorizontalAlignment','right','Fontsize',10);
					end;
				end;
				set(gca,'FontSize',12);
				%end final reformatting.
	        
				% shift back to main figure generation.
				figure(Main_fig);
				hold on;

				first_chr = false;
			end;
		end;
	end;


	%% ========================================================================
	% end stuff
	%==========================================================================
	%% Save figures.
	set(Main_fig,'PaperPosition',[0 0 8 6]*2);
	saveas(Main_fig,        [projectDir 'fig.CNV-SNP-map.RedGreen.1.eps'], 'epsc');
	saveas(Main_fig,        [projectDir 'fig.CNV-SNP-map.RedGreen.1.png'], 'png');
	delete(Main_fig);

	if (Linear_display == true)
		set(Linear_fig,'PaperPosition',[0 0 8 0.62222222]*2);
		saveas(Linear_fig, [projectDir 'fig.CNV-SNP-map.RedGreen.2.eps'], 'epsc');
		saveas(Linear_fig, [projectDir 'fig.CNV-SNP-map.RedGreen.2.png'], 'png');
		delete(Linear_fig);
	end;
elseif (useParent)
	% Dataset was compared to a parent, so don't draw a Red/Green alternate colors plot.
	fprintf(['\n##\n## Parent in use, so "CNV_SNP_hapmap_v5_RADseq_RedGreen.m" is being skipped.\n##\n']);
else
	% Dataset was not compared to a hapmap or parent, so don't draw a Red/Green alternate colors plot.
	fprintf(['\n##\n## Neither parent or hapmap in use, so "CNV_SNP_hapmap_v5_RADseq_RedGreen.m" is being skipped.\n##\n']);
end;

end
