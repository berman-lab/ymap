function [] = allelic_ratios_ddRADseq_D(main_dir,user,genomeUser,project,parent,hapmap,genome,ploidyEstimateString,ploidyBaseString,SNP_verString,LOH_verString,CNV_verString,displayBREAKS);
addpath('../');


%% =========================================================================================
% Load workspace variables saved in "allelic_ratios_ddRADseq_D.m"
%-------------------------------------------------------------------------------------------
projectDir  = [main_dir 'users/' user '/projects/' project '/'];
load([projectDir 'allelic_ratios_ddRADseq_B.workspace_variables.mat']);

if (useHapmap)
%
% Only run when compared vs. a hapmap.
%
	% Dataset was compared to a hapmap, so draw a Red/Green alternate colors plot.
	fprintf(['\n##\n## Hapmap in use, so "allelic_ratios_ddRADseq_D.m" is being processed.\n##\n']);

	%% ===============================================================================================
	% Setup for main figure generation.
	%-------------------------------------------------------------------------------------------------
	fig = figure(1);
	set(gcf, 'Position', [0 70 1024 600]);
	largestChr = find(chr_width == max(chr_width));


	%% ===============================================================================================
	% Setup for linear-view figure generation.
	%-------------------------------------------------------------------------------------------------
	if (Linear_display == true)
		Linear_fig = figure(2);
		Linear_genome_size   = sum(chr_size);
		Linear_Chr_max_width = 0.91;               % width for all chromosomes across figure.  1.00 - leftMargin - rightMargin - subfigure gaps.
		Linear_left_start    = 0.02;               % left margin (also right margin).
		Linear_left_chr_gap  = 0.07/(num_chrs-1);  % gaps between chr subfigures.
		Linear_height        = 0.6;
		Linear_base          = 0.1;
		Linear_TickSize      = -0.01;  %negative for outside, percentage of longest chr figure.
		maxY                 = 1; % ploidyBase*2;
		Linear_left          = Linear_left_start;
		axisLabelPosition_horiz = 0.01125;
	end;
	axisLabelPosition_vert = 0.01125;


	%% ===============================================================================================
	% Define colors for figure generation.
	%-------------------------------------------------------------------------------------------------
	% define colors for colorBars plot
	colorNoData     = [1.0   1.0   1.0  ]; %used when no data is available for the bin.
	colorInit       = [0.5   0.5   0.5  ]; %external; used in blending at ends of chr.
	% The 'project' is different than the 'hapmap'('parent').
	hom_color       = [1.0     0.0     0.0    ]; % homozygous, unphased.
	het_color       = [0.66667 0.66667 0.66667]; % heterozygous.
	oddHet_color    = [0.0     1.0     0.0    ]; % non-heterozygous data that isn't 100 hom.
	% haploid colors.
	color_1of1      = hom_color;
	% diploid colors.
	color_2of2      = hom_color;
	color_1of2      = het_color;
	% triploid colors.
	color_3of3      = hom_color;
	color_2of3      = oddHet_color;
	% tetraploid colors.
	color_4of4      = hom_color;
	color_3of4      = oddHet_color;
	color_2of4      = het_color;
	% pentaploid colors.
	color_5of5      = hom_color;
	color_4of5      = oddHet_color;
	color_3of5      = oddHet_color;
	% hexaploid colors.
	color_6of6      = hom_color;
	color_5of6      = oddHet_color;
	color_4of6      = oddHet_color;
	color_3of6      = het_color;
	% heptaploid colors.
	color_7of7      = hom_color;
	color_6of7      = oddHet_color;
	color_5of7      = oddHet_color;
	color_4of7      = oddHet_color;
	% octaploid colors.
	color_8of8      = hom_color;
	color_7of8      = oddHet_color;
	color_6of8      = oddHet_color;
	color_5of8      = oddHet_color;
	color_4of8      = het_color;
	% nonaploid colors.
	color_9of9      = hom_color;
	color_8of9      = oddHet_color;
	color_7of9      = oddHet_color;
	color_6of9      = oddHet_color;
	color_5of9      = oddHet_color;


	%%================================================================================================
	% Process SNP/hapmap data to determine colors for presentation.
	%-------------------------------------------------------------------------------------------------
	for chr = 1:num_chrs
		if (chr_in_use(chr) == 1)
			if (length(C_chr_count{chr}) > 1)
				%
				% Determining colors for each SNP coordinate.
				%
				for i = 1:length(C_chr_count{chr})
					pos                             = ceil(C_chr_SNP_data_positions{chr}(i)/new_bases_per_bin);
					localCopyEstimate               = round(CNVplot2{chr}(pos)*ploidy*ploidyAdjust);
					chr_SNPdata{chr,2}(pos)         = C_chr_SNP_data_ratios{ chr}(i);
					baseCall                        = C_chr_baseCall{        chr}{i};
					homologA                        = C_chr_SNP_homologA{    chr}{i};
					homologB                        = C_chr_SNP_homologB{    chr}{i};
					flipper                         = C_chr_SNP_flipHomologs{chr}(i);
					if (flipper == 1)                          % Variable 'flipper' value of '1' means the homologs are phased wrong.
						temp                    = homologA;
						homologA                = homologB;
						homologB                = temp;
					elseif (flipper == 10)                     % Variable 'flipper' value of '10' indicates no phasing information is available in the hapmap.
						baseCall                = 'Z';     % Variable 'baseCall' value of 'Z' will prevent either hapmap allele from matching and so unphased ratio colors will be used in the following section.
					end;
					allelicFraction                 = C_chr_SNP_data_ratios{chr}(i);
					if (localCopyEstimate <= 0);                colorList = colorNoData;
					elseif (localCopyEstimate == 1);            colorList = color_1of1;
					elseif (localCopyEstimate == 2)
						if (allelicFraction > 3/4);         colorList = color_2of2;
						else;                               colorList = color_1of2;
						end;
					elseif (localCopyEstimate == 3)
						if (allelicFraction > 5/6);         colorList = color_3of3;
						else;                               colorList = color_2of3;
						end;
					elseif (localCopyEstimate == 4)
						if (allelicFraction > 7/8);         colorList = color_4of4;
						elseif (allelicFraction > 5/8);     colorList = color_3of4;
						else;                               colorList = color_2of4;
						end;
					elseif (localCopyEstimate == 5)
						if (allelicFraction > 9/10);        colorList = color_5of5;
						elseif (allelicFraction > 7/10);    colorList = color_4of5;
						else;                               colorList = color_3of5;
						end;
					elseif (localCopyEstimate == 6)
						if (allelicFraction > 11/12);       colorList = color_6of6;
						elseif (allelicFraction > 9/12);    colorList = color_5of6;
						elseif (allelicFraction > 7/12);    colorList = color_4of6;
						else;                               colorList = color_3of6;
						end;
					elseif (localCopyEstimate == 7)
						if (allelicFraction > 13/14);       colorList = color_7of7;
						elseif (allelicFraction > 11/14);   colorList = color_6of7;
						elseif (allelicFraction > 9/14);    colorList = color_5of7;
						else;                               colorList = color_4of7;
						end;
					elseif (localCopyEstimate == 8)
						if (allelicFraction > 15/16);       colorList = color_8of8;
						elseif (allelicFraction > 13/16);   colorList = color_7of8;
						elseif (allelicFraction > 11/16);   colorList = color_6of8;
						elseif (allelicFraction > 9/16);    colorList = color_5of8;
						else;                               colorList = color_4of8;
						end;
					elseif (localCopyEstimate >= 9)
						if (allelicFraction > 17/18);       colorList = color_9of9;
						elseif (allelicFraction > 15/18);   colorList = color_8of9;
						elseif (allelicFraction > 13/18);   colorList = color_7of9;
						elseif (allelicFraction > 11/18);   colorList = color_6of9;
						else;                               colorList = color_5of9;
						end;
					end;

					chr_SNPdata_colorsC{chr,1}(pos) = chr_SNPdata_colorsC{chr,1}(pos) + colorList(1);
					chr_SNPdata_colorsC{chr,2}(pos) = chr_SNPdata_colorsC{chr,2}(pos) + colorList(2);
					chr_SNPdata_colorsC{chr,3}(pos) = chr_SNPdata_colorsC{chr,3}(pos) + colorList(3);
					chr_SNPdata_countC{ chr  }(pos) = chr_SNPdata_countC{ chr  }(pos) + 1;
				end;

				%
				% Average color per bin.
				%
				for pos = 1:length(chr_SNPdata_countC{chr})
					if (chr_SNPdata_countC{chr}(pos) > 0)
						chr_SNPdata_colorsC{chr,1}(pos) = chr_SNPdata_colorsC{chr,1}(pos)/chr_SNPdata_countC{chr}(pos);
						chr_SNPdata_colorsC{chr,2}(pos) = chr_SNPdata_colorsC{chr,2}(pos)/chr_SNPdata_countC{chr}(pos);
						chr_SNPdata_colorsC{chr,3}(pos) = chr_SNPdata_colorsC{chr,3}(pos)/chr_SNPdata_countC{chr}(pos);
					else
						chr_SNPdata_colorsC{chr,1}(pos) = 1.0;
						chr_SNPdata_colorsC{chr,2}(pos) = 1.0;
						chr_SNPdata_colorsC{chr,3}(pos) = 1.0;
					end;
				end;
			end;
		end;
	end;


	save([projectDir 'SNP_' SNP_verString '.reduced_RedGreen.mat'],'chr_SNPdata','new_bases_per_bin','chr_SNPdata_colorsC','chr_SNPdata_colorsP');

	%%================================================================================================
	% Make figures
	%-------------------------------------------------------------------------------------------------
	first_chr = true;
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

			% standard : This section sets the Y-axis labelling.
			axisLabelPosition = -50000/bases_per_bin;
			set(gca,'FontSize',12);
			if (chr == find(chr_posY == max(chr_posY)))
				title([ project ' allelic fraction map'],'Interpreter','none','FontSize',24);
			end;
			% standard : end axes labels etc.

			% standard : draw colorbars.
			for i = 1:ceil(chr_size(chr)/new_bases_per_bin)
				colorR   = chr_SNPdata_colorsC{chr,1}(i);
				colorG   = chr_SNPdata_colorsC{chr,2}(i);
				colorB   = chr_SNPdata_colorsC{chr,3}(i);
				if (colorR < 1) || (colorG < 1) || (colorB < 1)
					plot([i i], [0 maxY],'Color',[colorR colorG colorB]);
				end;
			end;
			% standard : end draw colorbars

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

				% linear : draw colorbars
				for i = 1:ceil(chr_size(chr)/new_bases_per_bin)
					colorR   = chr_SNPdata_colorsC{chr,1}(i);
					colorG   = chr_SNPdata_colorsC{chr,2}(i);
					colorB   = chr_SNPdata_colorsC{chr,3}(i);
					if (colorR < 1) || (colorG < 1) || (colorB < 1)
						plot([i i], [0 maxY],'Color',[colorR colorG colorB]);
					end;
				end;
				% linear : end draw colorbars

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
	saveas(fig,        [projectDir 'fig.allelic_ratio-map.RedGreen.c1.eps'], 'epsc');
	saveas(fig,        [projectDir 'fig.allelic_ratio-map.RedGreen.c1.png'], 'png');
	delete(fig);

	set(Linear_fig,'PaperPosition',[0 0 8 0.62222222]*2);
	saveas(Linear_fig, [projectDir 'fig.allelic_ratio-map.RedGreen.c2.eps'], 'epsc');
	saveas(Linear_fig, [projectDir 'fig.allelic_ratio-map.RedGreen.c2.png'], 'png');
	delete(Linear_fig);
elseif (useParent)
%
% Only run when compared vs. a parent.
%
	% Dataset was compared to a parent, so don't draw a Red/Green alternate colors plot.
	fprintf(['\n##\n## Parent in use, so "allelic_ratios_ddRADseq_D.m" is being skipped.\n##\n']);
else
%
% Only run when compared vs. itself.
%
	% Dataset was not compared to a hapmap or parent, so don't draw a Red/Green alternate colors plot.
	fprintf(['\n##\n## Neither parent or hapmap in use, so "allelic_ratios_ddRADseq_D.m" is being skipped.\n##\n']);
end;

end
