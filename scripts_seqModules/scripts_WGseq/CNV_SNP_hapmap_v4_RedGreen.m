function [] = CNV_SNP_hapmap_v4_RedGreen(main_dir,user,genomeUser,project,hapmap,genome,ploidyEstimateString,ploidyBaseString, ...
                                         SNP_verString,LOH_verString,CNV_verString,displayBREAKS);
addpath('../');

%% =========================================================================================
% Load workspace variables saved in "CNV_SNP_hapmap_v4.m"
%-------------------------------------------------------------------------------------------
projectDir  = [main_dir 'users/' user '/projects/' project '/'];
load([projectDir 'CNV_SNP_hapmap_v4.workspace_variables.mat']);

if (strcmp(project,hapmap) == 1)
	% Either a hapmap was in use or the project and parent are the same, so the
	% figure drawn by "CNV_SNP_hapmap_v4.m" is all that needs to be done.
	fprintf(['\n##\n## CNV_SNP_hapmap_v4_RedGreen.m is being skipped...\n']);
	fprintf(['##\tbecause the dataset is not being compared to another dataset.\n']);
else
	fprintf(['\n##\n## CNV_SNP_hapmap_v4_RedGreen.m is being processed.\n##\n']);

	%% =========================================================================================
	% Setup for figure generation.
	%-------------------------------------------------------------------------------------------
	fig = figure(1);
	set(gcf, 'Position', [0 70 1024 600]);


	%% =========================================================================================
	% Setup for linear-view figure generation.
	%-------------------------------------------------------------------------------------------
	if (Linear_display == true)
		Linear_fig           = figure(2);
		Linear_genome_size   = sum(chr_size);
		Linear_Chr_max_width = 0.91;               % width for all chromosomes across figure.  1.00 - leftMargin - rightMargin - subfigure gaps.
		Linear_left_start    = 0.02;               % left margin (also right margin).  (formerly 0.01)
		Linear_left_chr_gap  = 0.07/(num_chrs-1);  % gaps between chr subfigures.
		Linear_height        = 0.6;
		Linear_base          = 0.1;
		Linear_TickSize      = -0.01;  %negative for outside, percentage of longest chr figure.
		maxY                 = ploidyBase*2;
		Linear_left          = Linear_left_start;
		axisLabelPosition_horiz = -50000/bases_per_bin;
		axisLabelPosition_horiz = 0.01125;
	end;
	axisLabelPosition_vert = -50000/5000/2;
	axisLabelPosition_vert = 0.01125;


	%% =========================================================================================
	% Define colors for figure generation.
	%-------------------------------------------------------------------------------------------
	% define colors for colorBars plot
	colorNoData     = [1.0   1.0   1.0  ]; %used when no data is available for the bin.
	colorInit       = [0.5   0.5   0.5  ]; %external; used in blending at ends of chr.
	% The 'project' is different than the 'hapmap'('parent').
	hom_color       = [1.0     0.0     0.0    ]; % homozygous, unphased.
	het_color       = [0.66667 0.66667 0.66667]; % heterozygous.
	oddhet_color    = [0.0     1.0     0.0    ]; % non-heterozygous data that isn't 100 hom.
	% haploid colors.
	color_1of1      = hom_color;
	% diploid colors.
	color_2of2      = hom_color;
	color_1of2      = het_color;
	% triploid colors.
	color_3of3      = hom_color;
	color_2of3      = oddhet_color;
	% tetraploid colors.
	color_4of4      = hom_color;
	color_3of4      = oddhet_color;
	color_2of4      = het_color;
	% pentaploid colors.
	color_5of5      = hom_color;
	color_4of5      = oddhet_color;
	color_3of5      = oddhet_color;
	% hexaploid colors.
	color_6of6      = hom_color;
	color_5of6      = oddhet_color;
	color_4of6      = oddhet_color;
	color_3of6      = het_color;
	% heptaploid colors.
	color_7of7      = hom_color;
	color_6of7      = oddhet_color;
	color_5of7      = oddhet_color;
	color_4of7      = oddhet_color;
	% octaploid colors.
	color_8of8      = hom_color;
	color_7of8      = oddhet_color;
	color_6of8      = oddhet_color;
	color_5of8      = oddhet_color;
	color_4of8      = het_color;
	% nonaploid colors.
	color_9of9      = hom_color;
	color_8of9      = oddhet_color;
	color_7of9      = oddhet_color;
	color_6of9      = oddhet_color;
	color_5of9      = oddhet_color;


	%% =========================================================================================
	% Make figures
	%-------------------------------------------------------------------------------------------
	first_chr = true;
	for chr = 1:num_chrs
		if (chr_in_use(chr) == 1)
			figure(fig);
			% make standard chr cartoons.
			left          = chr_posX(chr);
			bottom        = chr_posY(chr);
			width         = chr_width(chr);
			height        = chr_height(chr);
			subPlotHandle = subplot('Position',[left bottom width height]);
			fprintf(['\tfigposition = [' num2str(left) ' | ' num2str(bottom) ' | ' num2str(width) ' | ' num2str(height) ']\n']);
			hold on;

			c_prev = colorInit;
			c_post = colorInit;
			c_     = c_prev;
			infill = zeros(1,length(unphased_plot2{chr}));
			colors = [];

			% standard : determines the color of each bin.
			for chr_bin = 1:length(SNPs_to_fullData_ratio{chr})+1;
				if (chr_bin-1 < length(SNPs_to_fullData_ratio{chr}))
					c_tot_post = SNPs_to_fullData_ratio{chr}(chr_bin)+SNPs_to_fullData_ratio{chr}(chr_bin);
					if (c_tot_post == 0)
						c_post = colorNoData;
						fprintf('.');
						if (mod(chr_bin,100) == 0);   fprintf('\n');   end;
					else
						% Define colorMix using localized copy number estimate to define SNP cutoff thresholds,
						%     then the ratio of SNP data in each SNP ratio bin.
						% For testing, consider all loci haploid, so only two ratio bins.
						ratioData_all       = [chr_SNPdata{chr,1}{chr_bin} chr_SNPdata{chr,2}{chr_bin}];
						coordinate_all      = [chr_SNPdata{chr,3}{chr_bin} chr_SNPdata{chr,4}{chr_bin}];
						colors_all          = cell(1,length(coordinate_all));

						% Determine localized copy number estimate, per bin.
						localCopyEstimate   = round(CNVplot2{chr}(chr_bin)*ploidy*ploidyAdjust);

						fprintf(num2str(localCopyEstimate));
						if (mod(chr_bin,100) == 0);   fprintf('\n');   end;

						if (localCopyEstimate <= 0)
							% Copy number indicates a deletion, no SNPs should exist.
							% Spurious SNP data will be drawn in white, the color used to represent a lack of SNP data.
							if (length(ratioData_all) > 0)
								for i = 1:length(ratioData_all)
									colors_all{i} = colorNoData;
								end;
							end;
						elseif (localCopyEstimate == 1)
							binCounts_all = 0;
							if (length(ratioData_all) > 0)
								for i = 1:length(ratioData_all)
									binCounts_all = binCounts_all+1;
									colors_all{i} = color_1of1;
								end;
							end;
							colorMix = color_1of1 * binCounts_all / sum(binCounts_all);
							% 1of1 contains [aa, bb].
						elseif (localCopyEstimate == 2)
							binCounts_all = zeros(1,2);
							if (length(ratioData_all) > 0)
								for i = 1:length(ratioData_all)
									if (ratioData_all(i) < 1/4);
										binCounts_all(1) = binCounts_all(1)+1;
										colors_all{i}    = color_2of2;
									elseif (ratioData_all(i) > 3/4);
										binCounts_all(1) = binCounts_all(1)+1;
										colors_all{i}    = color_2of2;
									else
										binCounts_all(2) = binCounts_all(2)+1;
										colors_all{i}    = color_1of2;
									end;
								end;
							end;
							colorMix = color_2of2 * binCounts_all(1)/ sum(binCounts_all) + ...
							           color_1of2 * binCounts_all(2)/ sum(binCounts_all);
							% 2of2 contains [aa, bb].
							% 1of2 contains [ab].
						elseif (localCopyEstimate == 3)
							binCounts_all = zeros(1,2);
							if (length(ratioData_all) > 0)
								for i = 1:length(ratioData_all)
									if (ratioData_all(i) < 1/6);
										binCounts_all(1) = binCounts_all(1)+1;
										colors_all{i}    = color_3of3;
									elseif (ratioData_all(i) > 5/6);
										binCounts_all(1) = binCounts_all(1)+1;
										colors_all{i}    = color_3of3;
									else
										binCounts_all(2) = binCounts_all(2)+1;
										colors_all{i}    = color_2of3;
									end;
								end;
							end;
							colorMix = color_3of3 * binCounts_all(1)/ sum(binCounts_all) + ...
							           color_2of3 * binCounts_all(2)/ sum(binCounts_all);
							% 3of3 contains [aaa, bbb].
							% 2of3 contains [abb, abb]
						elseif (localCopyEstimate == 4)
							binCounts_all = zeros(1,3);
							if (length(ratioData_all) > 0)
								for i = 1:length(ratioData_all)
									if (ratioData_all(i) < 1/8);
										binCounts_all(1) = binCounts_all(1)+1;
										colors_all{i}    = color_4of4;
									elseif (ratioData_all(i) < 3/8);
										binCounts_all(2) = binCounts_all(2)+1;
										colors_all{i}    = color_3of4;
									elseif (ratioData_all(i) > 7/8);
										binCounts_all(1) = binCounts_all(1)+1;
										colors_all{i}    = color_4of4;
									elseif (ratioData_all(i) > 5/8);
										binCounts_all(2) = binCounts_all(2)+1;
										colors_all{i}    = color_3of4;
									else
										binCounts_all(3) = binCounts_all(3)+1;
										colors_all{i}    = color_2of4;
									end;
								end;
							end;
							colorMix = color_4of4 * binCounts_all(1)/ sum(binCounts_all) + ...
							           color_3of4 * binCounts_all(2)/ sum(binCounts_all) + ...
							           color_2of4 * binCounts_all(3)/ sum(binCounts_all);
							% 4of4 contains [aaaa, bbbb].
							% 3of4 contains [aaab, abbb].
							% 2of4 contains [aabb].
						elseif (localCopyEstimate == 5)
							binCounts_all = zeros(1,3);
							if (length(ratioData_all) > 0)
								for i = 1:length(ratioData_all)
									if (ratioData_all(i) < 1/10);
										binCounts_all(1) = binCounts_all(1)+1;
										colors_all{i}    = color_5of5;
									elseif (ratioData_all(i) < 3/10);
										binCounts_all(2) = binCounts_all(2)+1;
										colors_all{i}    = color_4of5;
									elseif (ratioData_all(i) > 9/10);
										binCounts_all(1) = binCounts_all(1)+1;
										colors_all{i}    = color_5of5;
									elseif (ratioData_all(i) > 7/10);
										binCounts_all(2) = binCounts_all(2)+1;
										colors_all{i}    = color_4of5;
									else
										binCounts_all(3) = binCounts_all(3)+1;
										colors_all{i}    = color_3of5;
									end;
								end;
							end;
							colorMix = color_5of5 * binCounts_all(1)/ sum(binCounts_all) + ...
							           color_4of5 * binCounts_all(2)/ sum(binCounts_all) + ...
							           color_3of5 * binCounts_all(3)/ sum(binCounts_all);
							% 5of5 contains [aaaaa, bbbbb].
							% 4of5 contains [aaaab, abbbb].
							% 3of5 contains [aaabb, aabbb].
						elseif (localCopyEstimate == 6)
							binCounts_all = zeros(1,4);
							if (length(ratioData_all) > 0)
								for i = 1:length(ratioData_all)
									if (ratioData_all(i) < 1/12);
										binCounts_all(1) = binCounts_all(1)+1;
										colors_all{i}    = color_6of6;
									elseif (ratioData_all(i) < 3/12);
										binCounts_all(2) = binCounts_all(2)+1;
										colors_all{i}    = color_5of6;
									elseif (ratioData_all(i) < 5/12);
										binCounts_all(3) = binCounts_all(3)+1;
										colors_all{i}    = color_4of6;
									elseif (ratioData_all(i) > 11/12);
										binCounts_all(1) = binCounts_all(1)+1;
										colors_all{i}    = color_6of6;
									elseif (ratioData_all(i) > 9/12);
										binCounts_all(2) = binCounts_all(2)+1;
										colors_all{i}    = color_5of6;
									elseif (ratioData_all(i) > 7/12);
										binCounts_all(3) = binCounts_all(3)+1;
										colors_all{i}    = color_4of6;
									else
										binCounts_all(4) = binCounts_all(4)+1;
										colors_all{i}    = color_3of6;
									end;
								end;
							end;
							colorMix = color_6of6 * binCounts_all(1)/ sum(binCounts_all) + ...
							           color_5of6 * binCounts_all(2)/ sum(binCounts_all) + ...
							           color_4of6 * binCounts_all(3)/ sum(binCounts_all) + ...
							           color_3of6 * binCounts_all(4)/ sum(binCounts_all);
							% 6of6 contains [aaaaaa, bbbbbb].
							% 5of6 contains [aaaaab, abbbbb].
							% 4of6 contains [aaaabb, aabbbb].
							% 3of6 contains [aaabbb].
						elseif (localCopyEstimate == 7)
							binCounts_all = zeros(1,4);
							if (length(ratioData_all) > 0)
								for i = 1:length(ratioData_all)
									if (ratioData_all(i) < 1/14);
										binCounts_all(1) = binCounts_all(1)+1;
										colors_all{i}    = color_7of7;
									elseif (ratioData_all(i) < 3/14);
										binCounts_all(2) = binCounts_all(2)+1;
										colors_all{i}    = color_6of7;
									elseif (ratioData_all(i) < 5/14);
										binCounts_all(3) = binCounts_all(3)+1;
										colors_all{i}    = color_5of7;
									elseif (ratioData_all(i) > 13/14);
										binCounts_all(1) = binCounts_all(1)+1;
										colors_all{i}    = color_7of7;
									elseif (ratioData_all(i) > 11/14);
										binCounts_all(2) = binCounts_all(2)+1;
										colors_all{i}    = color_6of7;
									elseif (ratioData_all(i) > 9/14);
										binCounts_all(3) = binCounts_all(3)+1;
										colors_all{i}    = color_5of7;
									else
										binCounts_all(4) = binCounts_all(4)+1;
										colors_all{i}    = color_4of7;
									end;
								end;
							end;
							colorMix = color_7of7 * binCounts_all(1)/ sum(binCounts_all) + ...
							           color_6of7 * binCounts_all(2)/ sum(binCounts_all) + ...
							           color_5of7 * binCounts_all(3)/ sum(binCounts_all) + ...
							           color_4of7 * binCounts_all(4)/ sum(binCounts_all);
							% 7of7 contains [aaaaaaa, bbbbbbb].
							% 6of7 contains [aaaaaab, abbbbbb].
							% 5of7 contains [aaaaabb, aabbbbb].
							% 4of7 contains [aaaabbb, aaabbbb].
						elseif (localCopyEstimate == 8)
							binCounts_all = zeros(1,5);
							if (length(ratioData_all) > 0)
								for i = 1:length(ratioData_all)
									if (ratioData_all(i) < 1/16);
										binCounts_all(1) = binCounts_all(1)+1;
										colors_all{i}    = color_8of8;
									elseif (ratioData_all(i) < 3/16);
										binCounts_all(2) = binCounts_all(2)+1;
										colors_all{i}    = color_7of8;
									elseif (ratioData_all(i) < 5/16);
										binCounts_all(3) = binCounts_all(3)+1;
										colors_all{i}    = color_6of8;
									elseif (ratioData_all(i) < 6/16);
										binCounts_all(4) = binCounts_all(4)+1;
										colors_all{i}    = color_5of8;
									elseif (ratioData_all(i) > 15/16);
										binCounts_all(1) = binCounts_all(1)+1;
										colors_all{i}    = color_8of8;
									elseif (ratioData_all(i) > 13/16);
										binCounts_all(2) = binCounts_all(2)+1;
										colors_all{i}    = color_7of8;
									elseif (ratioData_all(i) > 11/16);
										binCounts_all(3) = binCounts_all(3)+1;
										colors_all{i}    = color_6of8;
									elseif (ratioData_all(i) > 9/16);
										binCounts_all(4) = binCounts_all(4)+1;
										colors_all{i}    = color_5of8;
									else
										binCounts_all(5) = binCounts_all(5)+1;
										colors_all{i}    = color_4of8;
									end;
								end;
							end;
							colorMix = color_8of8 * binCounts_all(1)/ sum(binCounts_all) + ...
							           color_7of8 * binCounts_all(2)/ sum(binCounts_all) + ...
							           color_6of8 * binCounts_all(3)/ sum(binCounts_all) + ...
							           color_5of8 * binCounts_all(4)/ sum(binCounts_all) + ...
							           color_4of8 * binCounts_all(5)/ sum(binCounts_all);
							% 8of8 contains [aaaaaaaa, bbbbbbbb].
							% 7of8 contains [aaaaaaab, abbbbbbb].
							% 6of8 contains [aaaaaabb, aabbbbbb].
							% 5of8 contains [aaaaabbb, aaabbbbb].
							% 4of8 contains [aaaabbbb].
						elseif (localCopyEstimate >= 9)
							binCounts_all = zeros(1,5);
							if (length(ratioData_all) > 0)
								for i = 1:length(ratioData_all)
									if (ratioData_all(i) < 1/18);
										binCounts_all(1)  = binCounts_all(1) +1;
										colors_all{i}     = color_9of9;
									elseif (ratioData_all(i) < 3/18);
										binCounts_all(2)  = binCounts_all(2) +1;
										colors_all{i}     = color_8of9;
									elseif (ratioData_all(i) < 5/18);
										binCounts_all(3)  = binCounts_all(3) +1;
										colors_all{i}     = color_7of9;
									elseif (ratioData_all(i) < 7/18);
										binCounts_all(4)  = binCounts_all(4) +1;
										colors_all{i}     = color_6of9;
									elseif (ratioData_all(i) > 17/18);
										binCounts_all(1)  = binCounts_all(1) +1;
										colors_all{i}     = color_9of9;
									elseif (ratioData_all(i) > 15/18);
										binCounts_all(2)  = binCounts_all(2) +1;
										colors_all{i}     = color_8of9;
									elseif (ratioData_all(i) > 13/18);
										binCounts_all(3)  = binCounts_all(3) +1;
										colors_all{i}     = color_7of9;
									elseif (ratioData_all(i) > 11/18);
										binCounts_all(4)  = binCounts_all(4) +1;
										colors_all{i}     = color_6of9;
									else
										binCounts_all(5) = binCounts_all(5)+1;
										colors_all{i}     = color_5of9;
									end;
								end;
							end;
							colorMix = color_9of9 * binCounts_all(1)/ sum(binCounts_all) + ...
							           color_8of9 * binCounts_all(2)/ sum(binCounts_all) + ...
							           color_7of9 * binCounts_all(3)/ sum(binCounts_all) + ...
							           color_6of9 * binCounts_all(4)/ sum(binCounts_all) + ...
							           color_5of9 * binCounts_all(5)/ sum(binCounts_all);
							% 9of9 contains [aaaaaaaaa, bbbbbbbbb].
							% 8of9 contains [aaaaaaaab, abbbbbbbb].
							% 7of9 contains [aaaaaaabb, aabbbbbbb].
							% 6of9 contains [aaaaaabbb, aaabbbbbb].
							% 5of9 contains [aaaaabbbb, aaaabbbbb].
						end;
						c_post =   colorMix   *   min(1,SNPs_to_fullData_ratio{chr}(chr_bin)) + colorNoData*(1-min(1,SNPs_to_fullData_ratio{chr}(chr_bin)));
					end;
				else
					c_post = colorInit;
				end;
				colors(chr_bin,1) = c_post(1);
				colors(chr_bin,2) = c_post(2);
				colors(chr_bin,3) = c_post(3);
			end;
		
			% standard : draw colorbars.
			for chr_bin = 1:length(phased_plot2{chr})+1;
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
			fprintf(['\nmain-plot : chr' num2str(chr) ':' num2str(length(CNVplot2{chr})) '\n']);
			fprintf(['ploidy     = ' num2str(ploidy)     '\n']);
			fprintf(['ploidyBase = ' num2str(ploidyBase) '\n']);
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
					endY = min(maxY,CNVhistValue*ploidy*ploidyAdjust);
				else
					endY = min(maxY,CNVhistValue*ploidy);
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
					line([0 x2], [maxY/4*1   maxY/4*1  ],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/4*3   maxY/4*3  ],'Color',[0.85 0.85 0.85]);
				case 3
					line([0 x2], [maxY/6*1   maxY/6*1  ],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/6*2   maxY/6*2  ],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/6*4   maxY/6*4  ],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/6*5   maxY/6*5  ],'Color',[0.85 0.85 0.85]);
				case 4
					line([0 x2], [maxY/8*1   maxY/8*1  ],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/8*2   maxY/8*2  ],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/8*3   maxY/8*3  ],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/8*5   maxY/8*5  ],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/8*6   maxY/8*6  ],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/8*7   maxY/8*7  ],'Color',[0.85 0.85 0.85]);
				case 5
					line([0 x2], [maxY/10*2  maxY/10*2 ],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/10*4  maxY/10*4 ],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/10*6  maxY/10*6 ],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/10*8  maxY/10*8 ],'Color',[0.85 0.85 0.85]);
				case 6
					line([0 x2], [maxY/12*2  maxY/12*2 ],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/12*4  maxY/12*4 ],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/12*8  maxY/12*8 ],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/12*10 maxY/12*10],'Color',[0.85 0.85 0.85]);
				case 7
					line([0 x2], [maxY/14*2  maxY/14*2 ],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/14*4  maxY/14*4 ],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/14*6  maxY/14*6 ],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/14*8  maxY/14*8 ],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/14*10 maxY/14*10],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/14*12 maxY/14*12],'Color',[0.85 0.85 0.85]);
				case 8
					line([0 x2], [maxY/16*2  maxY/16*2 ],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/16*4  maxY/16*4 ],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/16*6  maxY/16*6 ],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/16*10 maxY/16*10],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/16*12 maxY/16*12],'Color',[0.85 0.85 0.85]);
					line([0 x2], [maxY/16*14 maxY/16*14],'Color',[0.85 0.85 0.85]);
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
			set(gca,'TickLength',[(TickSize*chr_size(largestChr)/chr_size(chr)) 0]); %ensures same tick size on all subfigs.
			set(gca,'YTick',[]);
			set(gca,'YTickLabel',[]);
			set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
			set(gca,'XTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8','2.0','2.2','2.4','2.6','2.8','3.0','3.2'});
			text(-50000/5000/2*3, maxY/2,     chr_label{chr}, 'Rotation',90, 'HorizontalAlignment','center', 'VerticalAlign','bottom', 'Fontsize',20);

			switch ploidyBase
				case 1
					text(axisLabelPosition_vert, maxY/2,     '1','HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_vert, maxY,       '2','HorizontalAlignment','right','Fontsize',10);
				case 2
					text(axisLabelPosition_vert, maxY/4,     '1','HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_vert, maxY/4*2,   '2','HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_vert, maxY/4*3,   '3','HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_vert, maxY,       '4','HorizontalAlignment','right','Fontsize',10);
				case 3
					text(axisLabelPosition_vert, maxY/6*3,   '3','HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_vert, maxY,       '6','HorizontalAlignment','right','Fontsize',10);
				case 4
					text(axisLabelPosition_vert, maxY/8*2,   '2','HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_vert, maxY/8*4,   '4','HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_vert, maxY/8*6,   '6','HorizontalAlignment','right','Fontsize',10);
					text(axisLabelPosition_vert, maxY,       '8','HorizontalAlignment','right','Fontsize',10);
			end;
			set(gca,'FontSize',12);
			if (chr == find(chr_posY == max(chr_posY)))
				title([ project ' vs. (hapmap)' hapmap ' SNP/LOH map'],'Interpreter','none','FontSize',24);
			end;
			hold on;

			% standard : end axes labels etc.
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
					set(gca,'YTick',[]);
					set(gca,'XTick',[]);
					ylim([0,1]);
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
			% standard : end of chr copy number at right of the main chr cartons.

			% standard : places allelic fraction histogram to the left of the main chr cartoons.
			if (AnglePlot == true)
				width      = 0.075;
				height     = chr_height(chr);
				bottom     = chr_posY(chr);
				chr_length = chr_size(chr);
				for segment = 1:length(chrCopyNum{chr})
					fprintf(['^^^     segment#            = ' num2str(segment) ':' num2str(length(chrCopyNum{chr})) '\n']);
					if (segment == 1) % generate sublot for each segment.
						subplot('Position',[0.03 bottom width (height/length(chrCopyNum{chr}))]);
					else
						subplot('Position',[0.03 (bottom+height/length(chrCopyNum{chr})*(segment-1)) width (height/length(chrCopyNum{chr}))]);
					end;
					peaks               = chrSegment_peaks{chr,segment};
					mostLikelyGaussians = chrSegment_mostLikelyGaussians{chr,segment};
					actual_cutoffs      = chrSegment_actual_cutoffs{chr,segment};
					smoothed            = chrSegment_smoothed{chr,segment};

					hold on;
					segment_copyNum           = round(chrCopyNum{chr}(segment));  % copy number estimate of this segment.
					segment_chrBreaks         = chr_breaks{chr}(segment);         % break points of this segment.
					segment_smoothedHistogram = smoothed;                         % whole chromosome allelic ratio histogram smoothed.

					fprintf(['^^^     copyNum             = ' num2str(segment_copyNum)     '\n']);
					fprintf(['^^^     peaks               = ' num2str(peaks)               '\n']);
					fprintf(['^^^     mostLikelyGaussians = ' num2str(mostLikelyGaussians) '\n']);
					fprintf(['^^^     actual_cutoffs      = ' num2str(actual_cutoffs)      '\n']);

					copynum = round(chrCopyNum{chr}(segment));
					region_ = 0;
					for region = mostLikelyGaussians
						region_ = region_+1;

						% Define color of the histogram region.
						if (FillColors == true)			
							fprintf(['region_ #                = ' num2str(region_) '\n']);
							if (show_uncalibrated == true)
								color = het_color;
							else
								fprintf(['    copyNum              = ' num2str(copynum) '\n']);
								    if (copynum == 0) %deletion or error
								elseif (copynum == 1) %monosomy
									color = color_1of1;
									if (segment == 1)
										set(gca,'XTick',[0 200]);
										set(gca,'XTickLabel',{'a','b'});
									end;
								elseif (copynum == 2) %disomy
									if (region == 1);     color = color_2of2;
									elseif (region == 2); color = color_1of2;
									else                  color = color_2of2;
									end;
									if (segment == 1)
										set(gca,'XTick',0:100:200);
										set(gca,'XTickLabel',{'a','ab','b'});
									end;
								elseif (copynum == 3) %trisomy
									if (region == 1);     color = color_3of3;
									elseif (region == 2); color = color_2of3;
									elseif (region == 3); color = color_2of3;
									else                  color = color_3of3;
									end;
									if (segment == 1)
										set(gca,'XTick',[0 66.667 133.333 200]);
										set(gca,'XTickLabel',{'a','aab','abb','b'});
									end;
								elseif (copynum == 4) %tetrasomy
									if (region == 1);     color = color_4of4;
									elseif (region == 2); color = color_3of4;
									elseif (region == 3); color = color_2of4;
									elseif (region == 4); color = color_3of4;
									else                  color = color_4of4;
									end;
									if (segment == 1)
										set(gca,'XTick',0:50:200);
										set(gca,'XTickLabel',{'a', '3:1', '2:2', '1:3' 'b'});
									end;
								elseif (copynum == 5) %pentasomy
									if (region == 1);     color = color_5of5;
									elseif (region == 2); color = color_4of5;
									elseif (region == 3); color = color_3of5;
									elseif (region == 4); color = color_3of5;
									elseif (region == 5); color = color_4of5;
									else                  color = color_5of5;
									end;
									if (segment == 1)
										set(gca,'XTick',0:40:200);
										set(gca,'XTickLabel',{'a', '4:!', '3:2', '2:3', '1:4' 'b'});
									end;
								elseif (copynum == 6) %hexasomy
									if (region == 1);     color = color_6of6;
									elseif (region == 2); color = color_5of6;
									elseif (region == 3); color = color_4of6;
									elseif (region == 4); color = color_3of6;
									elseif (region == 5); color = color_4of6;
									elseif (region == 6); color = color_5of6;
									else                  color = color_6of6;
									end;
									if (segment == 1)
										set(gca,'XTick',0:33.333:200);
										set(gca,'XTickLabel',{'a', '5:1', '4:2', '3:3', '2:4', '1:5' 'b'});
									end;
								elseif (copynum == 7) %heptasomy
									if (region == 1);     color = color_7of7;
									elseif (region == 2); color = color_6of7;
									elseif (region == 3); color = color_5of7;
									elseif (region == 4); color = color_4of7;
									elseif (region == 5); color = color_4of7;
									elseif (region == 6); color = color_5of7;
									elseif (region == 7); color = color_6of7;
									else                  color = color_7of7;
									end;
									if (segment == 1)
										set(gca,'XTick',0:28.571:200);
										set(gca,'XTickLabel',{'a', '', '5:2', '', '', '2:5', '' 'b'});
									end;
								else % if (copynum == 8) %octasomy
									if (region == 1);     color = color_8of8;
									elseif (region == 2); color = color_7of8;
									elseif (region == 3); color = color_6of8;
									elseif (region == 4); color = color_5of8;
									elseif (region == 5); color = color_4of8;
									elseif (region == 6); color = color_5of8;
									elseif (region == 7); color = color_6of8;
									elseif (region == 8); color = color_7of8;
									else                  color = color_8of8;
									end;
									if (segment == 1)
										set(gca,'XTick',0:22.222:200);
										set(gca,'XTickLabel',{'a', '', '6:2', '', '4:4', '', '2:6', '' 'b'});
									end;
								end;
							end;
						else
							color = colorAB;
						end;

						fprintf(['    mostLikelyGaussian   = ' num2str(region) '\n']);
						if (length(mostLikelyGaussians) <= 1)
							% draw entire smoothed histogram.
							area(1:200,smoothed(1:200),'FaceColor',color,'EdgeColor',color);
						else
							% draw segment of smoothed histogram corresponding to region.
							if (region_ == 1) % first region in list.
								coord1 = round(200-actual_cutoffs(region_))+1;
								area(coord1:200, smoothed(coord1:200), 'FaceColor',color,'EdgeColor',color);
								fprintf(['    angleplotCoordinates = 200:' num2str(coord1) '\n']);
							elseif (region_ == length(mostLikelyGaussians)) % last region in list.
								coord2 = round(200-actual_cutoffs(region_-1))+1;
								area(1:coord2, smoothed(1:coord2), 'FaceColor',color,'EdgeColor',color);
								fprintf(['    angleplotCoordinate  = ' num2str(coord2) ':1\n']);
							else
								coord3 = round(200-actual_cutoffs(region_  ))+1;
								coord4 = round(200-actual_cutoffs(region_-1))+1;
								area(coord3:coord4, smoothed(coord3:coord4), 'FaceColor',color,'EdgeColor',color);
								fprintf(['    angleplotCoordinates = ' num2str(coord4) ':' num2str(coord3) '\n']);
							end;
						end;
						fprintf(['    color = ' num2str(color) '   (colorA = [' num2str(colorA) ']; colorB = [' num2str(colorB) '])\n']);
					end;

					colorPeak   = [0.5 0.5 0.5]; % color of lines drawn at peak locations.
					colorCutoff = [1.0 0.0 0.0]; % color of lines drawn at cutoffs between Gaussian fits. 
					for peak = 1:length(peaks)
						plot([200-peaks(peak); 200-peaks(peak)],[0; 1],'color',colorPeak);
					end;
					for cutoff = 1:length(actual_cutoffs)
						plot([200-actual_cutoffs(cutoff); 200-actual_cutoffs(cutoff)],[0; 1],'color',colorCutoff);
					end;
					set(gca,'FontSize',10);
					hold off;
					set(gca,'YTick',[]);
					if (segment ~= 1)
						set(gca,'XTick',[]);
					end;
					xlim([0,200]);
					ylim([0,1]);
				end;
			end;
			% standard : end of allelic fraction histogram at the left end of main chr cartoons.

			%% END of standard figure draw section.

			%% Linear figure draw section
			if (Linear_display == true)
				figure(Linear_fig);
				Linear_width = Linear_Chr_max_width*chr_size(chr)/Linear_genome_size;
				subplot('Position',[Linear_left Linear_base Linear_width Linear_height]);
				Linear_left = Linear_left + Linear_width + Linear_left_chr_gap;
				hold on;
				title(chr_label{chr},'Interpreter','none','FontSize',20);

				% linear : draw colorbars.
				for chr_bin = 1:length(unphased_plot2{chr})+1;
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
				% linear : end draw colorbars.

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
						endY = min(maxY,CNVhistValue*ploidy*ploidyAdjust);
					else
						endY = min(maxY,CNVhistValue*ploidy);
					end;
					y_ = [startY endY endY startY];
					% makes a blackbar for each bin.
					f = fill(x_,y_,c_);
					set(f,'linestyle','none');
				end;
				% linear : end CGH plot section.

				% linear : draw lines across plots for easier interpretation of CNV regions.
				x2 = chr_size(chr)/bases_per_bin;
				plot([0; x2], [maxY/2; maxY/2],'color',[0 0 0]);  % 2n line.
				switch ploidyBase
					case 1
					case 2
						line([0 x2], [maxY/4*1   maxY/4*1  ],'Color',[0.85 0.85 0.85]);
						line([0 x2], [maxY/4*3   maxY/4*3  ],'Color',[0.85 0.85 0.85]);
					case 3
						line([0 x2], [maxY/6*1   maxY/6*1  ],'Color',[0.85 0.85 0.85]);
						line([0 x2], [maxY/6*2   maxY/6*2  ],'Color',[0.85 0.85 0.85]);
						line([0 x2], [maxY/6*4   maxY/6*4  ],'Color',[0.85 0.85 0.85]);
						line([0 x2], [maxY/6*5   maxY/6*5  ],'Color',[0.85 0.85 0.85]);
					case 4
						line([0 x2], [maxY/8*1   maxY/8*1  ],'Color',[0.85 0.85 0.85]);
						line([0 x2], [maxY/8*2   maxY/8*2  ],'Color',[0.85 0.85 0.85]);
						line([0 x2], [maxY/8*3   maxY/8*3  ],'Color',[0.85 0.85 0.85]);
						line([0 x2], [maxY/8*5   maxY/8*5  ],'Color',[0.85 0.85 0.85]);
						line([0 x2], [maxY/8*6   maxY/8*6  ],'Color',[0.85 0.85 0.85]);
						line([0 x2], [maxY/8*7   maxY/8*7  ],'Color',[0.85 0.85 0.85]);
					case 5
						line([0 x2], [maxY/10*2  maxY/10*2 ],'Color',[0.85 0.85 0.85]);
						line([0 x2], [maxY/10*4  maxY/10*4 ],'Color',[0.85 0.85 0.85]);
						line([0 x2], [maxY/10*6  maxY/10*6 ],'Color',[0.85 0.85 0.85]);
						line([0 x2], [maxY/10*8  maxY/10*8 ],'Color',[0.85 0.85 0.85]);
					case 6
						line([0 x2], [maxY/12*2  maxY/12*2 ],'Color',[0.85 0.85 0.85]);
						line([0 x2], [maxY/12*4  maxY/12*4 ],'Color',[0.85 0.85 0.85]);
						line([0 x2], [maxY/12*8  maxY/12*8 ],'Color',[0.85 0.85 0.85]);
						line([0 x2], [maxY/12*10 maxY/12*10],'Color',[0.85 0.85 0.85]);
					case 7
						line([0 x2], [maxY/14*2  maxY/14*2 ],'Color',[0.85 0.85 0.85]);
						line([0 x2], [maxY/14*4  maxY/14*4 ],'Color',[0.85 0.85 0.85]);
						line([0 x2], [maxY/14*6  maxY/14*6 ],'Color',[0.85 0.85 0.85]);
						line([0 x2], [maxY/14*8  maxY/14*8 ],'Color',[0.85 0.85 0.85]);
						line([0 x2], [maxY/14*10 maxY/14*10],'Color',[0.85 0.85 0.85]);
						line([0 x2], [maxY/14*12 maxY/14*12],'Color',[0.85 0.85 0.85]);
					case 8
						line([0 x2], [maxY/16*2  maxY/16*2 ],'Color',[0.85 0.85 0.85]);
						line([0 x2], [maxY/16*4  maxY/16*4 ],'Color',[0.85 0.85 0.85]);
						line([0 x2], [maxY/16*6  maxY/16*6 ],'Color',[0.85 0.85 0.85]);
						line([0 x2], [maxY/16*10 maxY/16*10],'Color',[0.85 0.85 0.85]);
						line([0 x2], [maxY/16*12 maxY/16*12],'Color',[0.85 0.85 0.85]);
						line([0 x2], [maxY/16*14 maxY/16*14],'Color',[0.85 0.85 0.85]);
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
				% linear : end segmental aneuploidy breakpoint section.

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
				set(gca,'TickLength',[(Linear_TickSize*chr_size(largestChr)/chr_size(chr)) 0]); %ensures same tick size on all subfigs.
				set(gca,'YTick',[]);
				set(gca,'YTickLabel',[]);
				set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
				set(gca,'XTickLabel',[]);
				if (first_chr)
					% This section sets the Y-axis labelling.
					switch ploidyBase
						case 1
							text(axisLabelPosition_horiz, maxY/2,     '1','HorizontalAlignment','right','Fontsize',10);
							text(axisLabelPosition_horiz, maxY,       '2','HorizontalAlignment','right','Fontsize',10);
						case 2
							text(axisLabelPosition_horiz, maxY/4,     '1','HorizontalAlignment','right','Fontsize',10);
							text(axisLabelPosition_horiz, maxY/4*2,   '2','HorizontalAlignment','right','Fontsize',10);
							text(axisLabelPosition_horiz, maxY/4*3,   '3','HorizontalAlignment','right','Fontsize',10);
							text(axisLabelPosition_horiz, maxY,       '4','HorizontalAlignment','right','Fontsize',10);
						case 3
							text(axisLabelPosition_horiz, maxY/6*3,   '3','HorizontalAlignment','right','Fontsize',10);
							text(axisLabelPosition_horiz, maxY,       '6','HorizontalAlignment','right','Fontsize',10);
						case 4
							text(axisLabelPosition_horiz, maxY/8*2,   '2','HorizontalAlignment','right','Fontsize',10);
							text(axisLabelPosition_horiz, maxY/8*4,   '4','HorizontalAlignment','right','Fontsize',10);
							text(axisLabelPosition_horiz, maxY/8*6,   '6','HorizontalAlignment','right','Fontsize',10);
							text(axisLabelPosition_horiz, maxY,       '8','HorizontalAlignment','right','Fontsize',10);
					end;
				end;
				set(gca,'FontSize',12);
				% linear : end final reformatting.
	        
				% shift back to main figure generation.
				figure(fig);
				hold on;

				first_chr = false;
			end;
		end;
	end;

	%% ========================================================================
	% end stuff
	%==========================================================================

	fprintf('\n###\n### Saving main figure.\n###\n');
	set(   fig,        'PaperPosition',[0 0 8 6]*2);
	saveas(fig,        [projectDir 'fig.CNV-SNP-map.RedGreen.1.eps'], 'epsc');
	saveas(fig,        [projectDir 'fig.CNV-SNP-map.RedGreen.1.png'], 'png' );
	delete(fig);

	fprintf('\n###\n### Saving linear figure.\n###\n');
	set(   Linear_fig, 'PaperPosition',[0 0 8 0.62222222]*2);
	saveas(Linear_fig, [projectDir 'fig.CNV-SNP-map.RedGreen.2.eps'], 'epsc');
	saveas(Linear_fig, [projectDir 'fig.CNV-SNP-map.RedGreen.2.png'], 'png' );
	delete(Linear_fig);
end;

end
