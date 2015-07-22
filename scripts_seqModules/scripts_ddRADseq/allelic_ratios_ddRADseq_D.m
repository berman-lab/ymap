function [] = allelic_ratios_ddRADseq_D(main_dir,user,genomeUser,project,parent,hapmap,genome,ploidyEstimateString,ploidyBaseString,SNP_verString,LOH_verString,CNV_verString,displayBREAKS);
addpath('../');


%% =========================================================================================
% Load workspace variables saved in "allelic_ratios_ddRADseq_D.m"
%-------------------------------------------------------------------------------------------
projectDir  = [main_dir 'users/' user '/projects/' project '/'];
load([projectDir 'allelic_ratios_ddRADseq_B.workspace_variables.mat']);

if ((useHapmap) || (useParent))
	fprintf(['\n##\n## Hapmap in use, so "allelic_ratios_ddRADseq_D.m" is being processed.\n##\n']);

	%% =========================================================================================
	% Define colors for figure generation.
	%-------------------------------------------------------------------------------------------
	fprintf('\t|\tDefine colors used in figure generation.\n');
	phased_and_unphased_color_definitions;

	%% =========================================================================================
	% Calculate allelic fraction cutoffs for each chromosome and chromosome segment.
	%-------------------------------------------------------------------------------------------
	calculate_allelic_ratio_cutoffs;


	%%================================================================================================
	% Process SNP/hapmap data to determine colors for presentation.
	%-------------------------------------------------------------------------------------------------
	%%%% chr_SNPdata{chr,1}{pos} = phased SNP ratio data.
	%%%% chr_SNPdata{chr,2}{pos} = unphased SNP ratio data.
	%%%% chr_SNPdata{chr,3}{pos} = phased SNP position data.
	%%%% chr_SNPdata{chr,4}{pos} = unphased SNP position data.
	%%%% chr_SNPdata{chr,5}{pos} = flipper value for phased SNP.
	%%%% chr_SNPdata{chr,6}{pos} = flipper value for unphased SNP.
	for chr = 1:num_chrs
		if (chr_in_use(chr) == 1)
			if (length(C_chr_count{chr}) > 1)
				%
				% Determining colors for each SNP coordinate.
				%
				for SNP = 1:length(C_chr_count{chr})
					coordinate                      = C_chr_SNP_data_positions{chr}(SNP);
					pos                             = ceil(coordinate/new_bases_per_bin);
					localCopyEstimate               = round(CNVplot2{chr}(pos)*ploidy*ploidyAdjust);
					baseCall                        = C_chr_baseCall{        chr}{SNP};
					homologA                        = C_chr_SNP_homologA{    chr}{SNP};
					homologB                        = C_chr_SNP_homologB{    chr}{SNP};
					flipper                         = C_chr_SNP_flipHomologs{chr}(SNP);
					allelic_ratio                   = C_chr_SNP_data_ratios{ chr}(SNP);
					% Allelic ratio here is the ratio of the majority read call to all reads.
					% The consequence of this is that it will always be on the range [0.5 .. 1.0].
					if (flipper == 10)                         % Variable 'flipper' value of '10' indicates no phasing information is available in the hapmap.
						baseCall                = 'Z';     % Variable 'baseCall' value of 'Z' will prevent either hapmap allele from matching and so unphased ratio colors will be used in the following section.
						chr_SNPdata{chr,2}{pos} = [chr_SNPdata{chr,2}{pos} allelic_ratio 1-allelic_ratio];
						chr_SNPdata{chr,4}{pos} = [chr_SNPdata{chr,4}{pos} coordinate    coordinate     ];
						chr_SNPdata{chr,6}{pos} = [chr_SNPdata{chr,6}{pos} flipper       flipper        ];
					elseif (flipper == 1)
						temp                    = homologA;
						homologA                = homologB;
						homologB                = temp;
						if (baseCall == homologA)
							allelic_ratio = 1-allelic_ratio;
						end;
						chr_SNPdata{chr,1}{pos} = [chr_SNPdata{chr,1}{pos} allelic_ratio];
						chr_SNPdata{chr,3}{pos} = [chr_SNPdata{chr,3}{pos} coordinate   ];
						chr_SNPdata{chr,5}{pos} = [chr_SNPdata{chr,5}{pos} flipper      ];
					else % (flipper == 0)
						if (baseCall == homologA)
							allelic_ratio = 1-allelic_ratio;
						end;
						chr_SNPdata{chr,1}{pos} = [chr_SNPdata{chr,1}{pos} allelic_ratio];
						chr_SNPdata{chr,3}{pos} = [chr_SNPdata{chr,3}{pos} coordinate   ];
						chr_SNPdata{chr,5}{pos} = [chr_SNPdata{chr,5}{pos} flipper      ];
					end;

					% identify the segment containing the SNP.
					segmentID                       = 0;
					for segment = 1:(length(chrCopyNum{chr}))
						segment_start           = chr_breaks{chr}(segment  )*chr_size(chr);
						segment_end             = chr_breaks{chr}(segment+1)*chr_size(chr);
						if (coordinate > segment_start) && (coordinate <= segment_end)
							segmentID       = segment;
						end;
					end;

					% Load cutoffs between Gaussian fits performed earlier.
					segment_copyNum                 = round(chrCopyNum{              chr}(segmentID));
					actual_cutoffs                  = chrSegment_actual_cutoffs{     chr}{segmentID};
					mostLikelyGaussians             = chrSegment_mostLikelyGaussians{chr}{segmentID};
					SNPratio_int                    = (allelic_ratio)*199+1;

					% Identify the allelic ratio region containing the SNP.
					cutoffs                         = [1 actual_cutoffs 200];
					ratioRegionID                   = 0;
					for GaussianRegionID = 1:length(mostLikelyGaussians)
						cutoff_start            = cutoffs(GaussianRegionID  );
						cutoff_end              = cutoffs(GaussianRegionID+1);
						if (GaussianRegionID == 1)
							if (SNPratio_int >= cutoff_start) && (SNPratio_int <= cutoff_end)
								ratioRegionID   = mostLikelyGaussians(GaussianRegionID);
							end;
						else
							if (SNPratio_int > cutoff_start) && (SNPratio_int <= cutoff_end)
								ratioRegionID   = mostLikelyGaussians(GaussianRegionID);
							end;
						end;
					end;

					allelicFraction                = C_chr_SNP_data_ratios{chr}(SNP);
					if (segment_copyNum <= 0);              colorList = colorNoData;
					elseif (segment_copyNum == 1)
															colorList = alternate_color_1of1;
					elseif (segment_copyNum == 2)
						if (ratioRegionID == 3);            colorList = alternate_color_2of2;
						elseif (ratioRegionID == 2);        colorList = alternate_color_1of2;
						else                                colorList = alternate_color_2of2;
						end;
					elseif (segment_copyNum == 3)
						if (ratioRegionID == 4);            colorList = alternate_color_3of3;
						elseif (ratioRegionID == 3);        colorList = alternate_color_2of3;
						elseif (ratioRegionID == 2);        colorList = alternate_color_2of3;
						else                                colorList = alternate_color_3of3;
						end;
					elseif (segment_copyNum == 4)
						if (ratioRegionID == 5);            colorList = alternate_color_4of4;
						elseif (ratioRegionID == 4);        colorList = alternate_color_3of4;
						elseif (ratioRegionID == 3);        colorList = alternate_color_2of4;
						elseif (ratioRegionID == 2);        colorList = alternate_color_3of4;
						else                                colorList = alternate_color_4of4;
						end;
					elseif (segment_copyNum == 5)
						if (ratioRegionID == 6);            colorList = alternate_color_5of5;
						elseif (ratioRegionID == 5);        colorList = alternate_color_4of5;
						elseif (ratioRegionID == 4);        colorList = alternate_color_3of5;
						elseif (ratioRegionID == 3);        colorList = alternate_color_3of5;
						elseif (ratioRegionID == 2);        colorList = alternate_color_4of5;
						else                                colorList = alternate_color_5of5;
						end;
					elseif (segment_copyNum == 6)
						if (ratioRegionID == 7);            colorList = alternate_color_6of6;
						elseif (ratioRegionID == 6);        colorList = alternate_color_5of6;
						elseif (ratioRegionID == 5);        colorList = alternate_color_4of6;
						elseif (ratioRegionID == 4);        colorList = alternate_color_3of6;
						elseif (ratioRegionID == 3);        colorList = alternate_color_4of6;
						elseif (ratioRegionID == 2);        colorList = alternate_color_5of6;
						else                                colorList = alternate_color_6of6;
						end;
					elseif (segment_copyNum == 7)
						if (ratioRegionID == 8);            colorList = alternate_color_7of7;
						elseif (ratioRegionID == 7);        colorList = alternate_color_6of7;
						elseif (ratioRegionID == 6);        colorList = alternate_color_5of7;
						elseif (ratioRegionID == 5);        colorList = alternate_color_4of7;
						elseif (ratioRegionID == 3);        colorList = alternate_color_4of7;
						elseif (ratioRegionID == 3);        colorList = alternate_color_5of7;
						elseif (ratioRegionID == 2);        colorList = alternate_color_6of7;
						else                                colorList = alternate_color_7of7;
						end;
					elseif (segment_copyNum == 8)
						if (ratioRegionID == 9);            colorList = alternate_color_8of8;
						elseif (ratioRegionID == 8);        colorList = alternate_color_7of8;
						elseif (ratioRegionID == 7);        colorList = alternate_color_6of8;
						elseif (ratioRegionID == 6);        colorList = alternate_color_5of8;
						elseif (ratioRegionID == 5);        colorList = alternate_color_4of8;
						elseif (ratioRegionID == 4);        colorList = alternate_color_5of8;
						elseif (ratioRegionID == 3);        colorList = alternate_color_6of8;
						elseif (ratioRegionID == 2);        colorList = alternate_color_7of8;
						else                                colorList = alternate_color_8of8;
						end;
					elseif (segment_copyNum >= 9)
						if (ratioRegionID == 10);           colorList = alternate_color_9of9;
						elseif (ratioRegionID == 9);        colorList = alternate_color_8of9;
						elseif (ratioRegionID == 8);        colorList = alternate_color_7of9;
						elseif (ratioRegionID == 7);        colorList = alternate_color_6of9;
						elseif (ratioRegionID == 6);        colorList = alternate_color_5of9;
						elseif (ratioRegionID == 5);        colorList = alternate_color_5of9;
						elseif (ratioRegionID == 4);        colorList = alternate_color_6of9;
						elseif (ratioRegionID == 3);        colorList = alternate_color_7of9;
						elseif (ratioRegionID == 2);        colorList = alternate_color_8of9;
						else                                colorList = alternate_color_9of9;
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


	%% ===============================================================================================
	% Setup for main figure generation.
	%-------------------------------------------------------------------------------------------------
	Main_fig = figure();
	set(gcf, 'Position', [0 70 1024 600]);
	largestChr = find(chr_width == max(chr_width));
	largestChr = largestChr(1);


	%% ===============================================================================================
	% Setup for linear-view figure generation.
	%-------------------------------------------------------------------------------------------------
	if (Linear_display == true)
		Linear_fig = figure();
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


	%%================================================================================================
	% Make figures
	%-------------------------------------------------------------------------------------------------
	first_chr = true;
	for chr = 1:num_chrs
		if (chr_in_use(chr) == 1)
			figure(Main_fig);

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
				figure(Main_fig);
				first_chr = false;
			end;
		end;
	end;

	%% Save figures.
	set(Main_fig,'PaperPosition',[0 0 8 6]*2);
	saveas(Main_fig,   [projectDir 'fig.allelic_ratio-map.RedGreen.c1.eps'], 'epsc');
	saveas(Main_fig,   [projectDir 'fig.allelic_ratio-map.RedGreen.c1.png'], 'png');
	delete(Main_fig);

	set(Linear_fig,'PaperPosition',[0 0 8 0.62222222]*2);
	saveas(Linear_fig, [projectDir 'fig.allelic_ratio-map.RedGreen.c2.eps'], 'epsc');
	saveas(Linear_fig, [projectDir 'fig.allelic_ratio-map.RedGreen.c2.png'], 'png');
	delete(Linear_fig);
elseif (useParent)
	% Dataset was compared to a parent, so don't draw a Red/Green alternate colors plot.
	fprintf(['\n##\n## Parent in use, so "allelic_ratios_ddRADseq_D.m" is being skipped.\n##\n']);
else
	% Dataset was not compared to a hapmap or parent, so don't draw a Red/Green alternate colors plot.
	fprintf(['\n##\n## Neither parent or hapmap in use, so "allelic_ratios_ddRADseq_D.m" is being skipped.\n##\n']);
end;

end
