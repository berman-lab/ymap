function [] = CNV_SNP_hapmap_v4_RedGreen(main_dir,user,genomeUser,project,hapmap,genome,ploidyEstimateString,ploidyBaseString,SNP_verString,LOH_verString,CNV_verString,displayBREAKS);
addpath('../');
workingDir      = [main_dir 'users/' user '/projects/' project '/'];


fprintf('\n');
fprintf('##################################\n');
fprintf('## CNV_SNP_hapmap_v4_RedGreen.m ##\n');
fprintf('##################################\n');


%% =========================================================================================
% Load workspace variables saved in "CNV_SNP_hapmap_v4.m"
%-------------------------------------------------------------------------------------------
projectDir  = [main_dir 'users/' user '/projects/' project '/'];
load([projectDir 'CNV_SNP_hapmap_v4.workspace_variables.mat']);


if (~useHapmap) && (~useParent)
	fprintf(['\n##\n## CNV_SNP_hapmap_v4_RedGreen.m is being skipped...\n']);
	fprintf(['##\tbecause the dataset is not being compared to another dataset.\n']);
else
	fprintf(['\n##\n## CNV_SNP_hapmap_v4_RedGreen.m is being processed.\n##\n']);


	%% =========================================================================================
	% Define alternate color scheme for figure generation.
	%-------------------------------------------------------------------------------------------
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
	% Process SNP/hapmap data to determine colors to be presented for each SNP locus.
	%-------------------------------------------------------------------------------------------------

	%% =========================================================================================
	% Calculate allelic fraction cutoffs for each segment and populate data structure containing
	% SNP phasing information.
	%       chr_SNPdata{chr,1}{chr_bin} = phased SNP ratio data.
	%       chr_SNPdata{chr,2}{chr_bin} = unphased SNP ratio data.
	%       chr_SNPdata{chr,3}{chr_bin} = phased SNP position data.
	%       chr_SNPdata{chr,4}{chr_bin} = unphased SNP position data.
	%       chr_SNPdata{chr,5}{chr_bin} = phased SNP allele strings.   (baseCall:alleleA/alleleB)
	%       chr_SNPdata{chr,6}{chr_bin} = unphased SNP allele strings.
	%-------------------------------------------------------------------------------------------
	% Prepare data for "calculate_allelic_ratio_cutoffs.m".
	temp_holding = chr_SNPdata;
	calculate_allelic_ratio_cutoffs;
	chr_SNPdata = temp_holding;

	%% =========================================================================================
	% Define new colors for SNPs, using Gaussian fitting crossover points as ratio cutoffs.
	%-------------------------------------------------------------------------------------------
	for chr = 1:num_chrs
		if (chr_in_use(chr) == 1)
			for chr_bin = 1:ceil(chr_size(chr)/bases_per_bin);
				%
				% Determining colors for each SNP coordinate from calculated cutoffs.
				%
				localCopyEstimate                       = round(CNVplot2{chr}(chr_bin)*ploidy*ploidyAdjust);
				allelic_ratios                          = [chr_SNPdata{chr,1}{chr_bin} chr_SNPdata{chr,2}{chr_bin}];
				coordinates                             = [chr_SNPdata{chr,3}{chr_bin} chr_SNPdata{chr,4}{chr_bin}];
				allele_strings                          = [chr_SNPdata{chr,5}{chr_bin} chr_SNPdata{chr,6}{chr_bin}];

				if (length(allelic_ratios) > 0)
					for SNP = 1:length(allelic_ratios)
						% Load phased SNP data from earlier defined structure.
						allelic_ratio                   = allelic_ratios(SNP);
						coordinate                      = coordinates(SNP);
						if (length(allelic_ratios) > 1)
							allele_string                   = allele_strings{SNP};
						else
							allele_string                   = allele_strings;
						end;
						baseCall                        = allele_string(1);
						homologA                        = allele_string(3);
						homologB                        = allele_string(5);

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

						% Calculate allelic ratio on range of [1..200].
						SNPratio_int                    = (allelic_ratio)*199+1;

						% Identify the allelic ratio region containing the SNP.
						cutoffs                         = [1 actual_cutoffs 200];
						ratioRegionID                   = 0;
						for GaussianRegionID = 1:length(mostLikelyGaussians)
							cutoff_start            = cutoffs(GaussianRegionID  );
							cutoff_end              = cutoffs(GaussianRegionID+1);
							if (GaussianRegionID == 1)
								if (SNPratio_int >= cutoff_start) && (SNPratio_int <= cutoff_end)
									ratioRegionID = mostLikelyGaussians(GaussianRegionID);
								end;
							else
								if (SNPratio_int > cutoff_start) && (SNPratio_int <= cutoff_end)
									ratioRegionID = mostLikelyGaussians(GaussianRegionID);
								end;
							end;
						end;

						if (segment_copyNum <= 0);                          colorList = colorNoData;
						elseif (segment_copyNum == 1)
							% allelic fraction cutoffs: [0.50000] => [A B]
							colorList = unphased_color_1of1;
						elseif (segment_copyNum == 2)
							%   allelic fraction cutoffs: [0.25000 0.75000] => [AA AB BB]
							if (ratioRegionID == 3);            colorList = unphased_color_2of2;
							elseif (ratioRegionID == 2);        colorList = unphased_color_1of2;
							else                                colorList = unphased_color_2of2;
							end;
						elseif (segment_copyNum == 3)
							% allelic fraction cutoffs: [0.16667 0.50000 0.83333] => [AAA AAB ABB BBB]
							if (ratioRegionID == 4);            colorList = unphased_color_3of3;
							elseif (ratioRegionID == 3);        colorList = unphased_color_2of3;
							elseif (ratioRegionID == 2);        colorList = unphased_color_2of3;
							else                                colorList = unphased_color_3of3;
							end;
						elseif (segment_copyNum == 4)
							% allelic fraction cutoffs: [0.12500 0.37500 0.62500 0.87500] => [AAAA AAAB AABB ABBB BBBB]
							if (ratioRegionID == 5);            colorList = unphased_color_4of4;
							elseif (ratioRegionID == 4);        colorList = unphased_color_3of4;
							elseif (ratioRegionID == 3);        colorList = unphased_color_2of4;
							elseif (ratioRegionID == 2);        colorList = unphased_color_3of4;
							else                                colorList = unphased_color_4of4;
							end;
						elseif (segment_copyNum == 5)
							% allelic fraction cutoffs: [0.10000 0.30000 0.50000 0.70000 0.90000] => [AAAAA AAAAB AAABB AABBB ABBBB BBBBB]
							if (ratioRegionID == 6);            colorList = unphased_color_5of5;
							elseif (ratioRegionID == 5);        colorList = unphased_color_4of5;
							elseif (ratioRegionID == 4);        colorList = unphased_color_3of5;
							elseif (ratioRegionID == 3);        colorList = unphased_color_3of5;
							elseif (ratioRegionID == 2);        colorList = unphased_color_4of5;
							else                                colorList = unphased_color_5of5;
							end;
						elseif (segment_copyNum == 6)
							% allelic fraction cutoffs: [0.08333 0.25000 0.41667 0.58333 0.75000 0.91667] => [AAAAAA AAAAAB AAAABB AAABBB AABBBB ABBBBB BBBBBB]
							if (ratioRegionID == 7);            colorList = unphased_color_6of6;
							elseif (ratioRegionID == 6);        colorList = unphased_color_5of6;
							elseif (ratioRegionID == 5);        colorList = unphased_color_4of6;
							elseif (ratioRegionID == 4);        colorList = unphased_color_3of6;
							elseif (ratioRegionID == 3);        colorList = unphased_color_4of6;
							elseif (ratioRegionID == 2);        colorList = unphased_color_5of6;
							else                                colorList = unphased_color_6of6;
							end;
						elseif (segment_copyNum == 7)
							% allelic fraction cutoffs: [0.07143 0.21429 0.35714 0.50000 0.64286 0.78571 0.92857] => [AAAAAAA AAAAAAB AAAAABB AAAABBB AAABBBB AABBBBB ABBBBBB BBBBBBB]
							if (ratioRegionID == 8);            colorList = unphased_color_7of7;
							elseif (ratioRegionID == 7);        colorList = unphased_color_6of7;
							elseif (ratioRegionID == 6);        colorList = unphased_color_5of7;
							elseif (ratioRegionID == 5);        colorList = unphased_color_4of7;
							elseif (ratioRegionID == 3);        colorList = unphased_color_4of7;
							elseif (ratioRegionID == 3);        colorList = unphased_color_5of7;
							elseif (ratioRegionID == 2);        colorList = unphased_color_6of7;
							else                                colorList = unphased_color_7of7;
							end;
						elseif (segment_copyNum == 8)
							% allelic fraction cutoffs: [0.06250 0.18750 0.31250 0.43750 0.56250 0.68750 0.81250 0.93750] => [AAAAAAAA AAAAAAAB AAAAAABB AAAAABBB AAAABBBB AAABBBBB AABBBBBB ABBBBBBB BBBBBBBB]
							if (ratioRegionID == 9);            colorList = unphased_color_8of8;
							elseif (ratioRegionID == 8);        colorList = unphased_color_7of8;
							elseif (ratioRegionID == 7);        colorList = unphased_color_6of8;
							elseif (ratioRegionID == 6);        colorList = unphased_color_5of8;
							elseif (ratioRegionID == 5);        colorList = unphased_color_4of8;
							elseif (ratioRegionID == 4);        colorList = unphased_color_5of8;
							elseif (ratioRegionID == 3);        colorList = unphased_color_6of8;
							elseif (ratioRegionID == 2);        colorList = unphased_color_7of8;
							else                                colorList = unphased_color_8of8;
							end;
						elseif (segment_copyNum >= 9)
							% allelic fraction cutoffs: [0.05556 0.16667 0.27778 0.38889 0.50000 0.61111 0.72222 0.83333 0.94444] => [AAAAAAAAA AAAAAAAAB AAAAAAABB AAAAAABBB AAAAABBBB AAAABBBBB AAABBBBBB AABBB$
							%                                                                                                         ABBBBBBBB BBBBBBBBB]
							if (ratioRegionID == 10);           colorList = unphased_color_9of9;
							elseif (ratioRegionID == 9);        colorList = unphased_color_8of9;
							elseif (ratioRegionID == 8);        colorList = unphased_color_7of9;
							elseif (ratioRegionID == 7);        colorList = unphased_color_6of9;
							elseif (ratioRegionID == 6);        colorList = unphased_color_5of9;
							elseif (ratioRegionID == 5);        colorList = unphased_color_5of9;
							elseif (ratioRegionID == 4);        colorList = unphased_color_6of9;
							elseif (ratioRegionID == 3);        colorList = unphased_color_7of9;
							elseif (ratioRegionID == 2);        colorList = unphased_color_8of9;
							else                                colorList = unphased_color_9of9;
							end;
						end;
						chr_SNPdata_colorsC{chr,1}(chr_bin) = chr_SNPdata_colorsC{chr,1}(chr_bin) + colorList(1);
						chr_SNPdata_colorsC{chr,2}(chr_bin) = chr_SNPdata_colorsC{chr,2}(chr_bin) + colorList(2);
						chr_SNPdata_colorsC{chr,3}(chr_bin) = chr_SNPdata_colorsC{chr,3}(chr_bin) + colorList(3);
						chr_SNPdata_countC{ chr  }(chr_bin) = chr_SNPdata_countC{ chr  }(chr_bin) + 1;

						% Troubleshooting output.
						% fprintf(['chr = ' num2str(chr) '; seg = ' num2str(segment) '; bin = ' num2str(chr_bin) '; ratioRegionID = ' num2str(ratioRegionID) '\n']);
					end;
				end;

				%
				% Average colors of SNPs found in bin.
				%
				if (length(allelic_ratios) > 0)
					if (chr_SNPdata_countC{chr}(chr_bin) > 0)
						chr_SNPdata_colorsC{chr,1}(chr_bin) = chr_SNPdata_colorsC{chr,1}(chr_bin)/chr_SNPdata_countC{chr}(chr_bin);
						chr_SNPdata_colorsC{chr,2}(chr_bin) = chr_SNPdata_colorsC{chr,2}(chr_bin)/chr_SNPdata_countC{chr}(chr_bin);
						chr_SNPdata_colorsC{chr,3}(chr_bin) = chr_SNPdata_colorsC{chr,3}(chr_bin)/chr_SNPdata_countC{chr}(chr_bin);
					else
						chr_SNPdata_colorsC{chr,1}(chr_bin) = 1.0;
						chr_SNPdata_colorsC{chr,2}(chr_bin) = 1.0;
						chr_SNPdata_colorsC{chr,3}(chr_bin) = 1.0;
					end;
				else
					chr_SNPdata_colorsC{chr,1}(chr_bin) = 1.0;
					chr_SNPdata_colorsC{chr,2}(chr_bin) = 1.0;
					chr_SNPdata_colorsC{chr,3}(chr_bin) = 1.0;
				end;
			end;
		end;
	end;


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
				for chr_bin = 1:ceil(chr_size(chr)/bases_per_bin)
					% Regenerate chr plot data if the save file does not exist.

					% the number of heterozygous data points in this bin.
					SNPs_count{chr}(chr_bin)                                     = length(chr_SNPdata{chr,1}{chr_bin}) + length(chr_SNPdata{chr,2}{chr_bin});

					% divide by the threshold for full color saturation in SNP/LOH figure.
					SNPs_to_fullData_ratio{chr}(chr_bin)                         = SNPs_count{chr}(chr_bin)/full_data_threshold;

					% any bins with more data than the threshold for full color saturation are limited to full saturation.
					SNPs_to_fullData_ratio{chr}(SNPs_to_fullData_ratio{chr} > 1) = 1;

					phased_plot{chr}(chr_bin)                    = length(chr_SNPdata{chr,1}{chr_bin});             % phased data.
					phased_plot2{chr}(chr_bin)                   = phased_plot{chr}(chr_bin)/full_data_threshold;   %
					phased_plot2{chr}(phased_plot2{chr} > 1)     = 1;                                               %

					unphased_plot{chr}(chr_bin)                  = length(chr_SNPdata{chr,2}{chr_bin});             % unphased data.
					unphased_plot2{chr}(chr_bin)                 = unphased_plot{chr}(chr_bin)/full_data_threshold; %
					unphased_plot2{chr}(unphased_plot2{chr} > 1) = 1;                                               %
				end;
			end;
		end;
	end;
	fprintf('\n');
	largestChr = find(chr_width == max(chr_width));
	largestChr = largestChr(1);


	%% =========================================================================================
	% Setup for main figure generation.
	%-------------------------------------------------------------------------------------------
	fig = figure(1);
	set(gcf, 'Position', [0 70 1024 600]);


	%% =========================================================================================
	% Setup for linear-view figure generation.
	%-------------------------------------------------------------------------------------------
	if (Linear_display == true)
		Linear_fig              = figure(2);
		Linear_genome_size      = sum(chr_size);
		Linear_Chr_max_width    = 0.91;               % width for all chromosomes across figure.  1.00 - leftMargin - rightMargin - subfigure gaps.
		Linear_left_start       = 0.02;               % left margin (also right margin).  (formerly 0.01)
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


			%% standard : determine color of each bin.
			for chr_bin = 1:ceil(chr_size(chr)/bases_per_bin)
				c_tot_post = SNPs_to_fullData_ratio{chr}(chr_bin)+SNPs_to_fullData_ratio{chr}(chr_bin);
				if (c_tot_post == 0)
					c_post = colorNoData;
					fprintf('.');
					if (mod(chr_bin,100) == 0);   fprintf('\n');   end;
				else
					% Average of SNP position colors defined earlier.
					colorMix = [chr_SNPdata_colorsC{chr,1}(chr_bin) chr_SNPdata_colorsC{chr,2}(chr_bin) chr_SNPdata_colorsC{chr,3}(chr_bin)];

					% Determine color to draw bin, accounting for limited data and data saturation.
					c_post =   colorMix   *   min(1,SNPs_to_fullData_ratio{chr}(chr_bin)) + ...
					           colorNoData*(1-min(1,SNPs_to_fullData_ratio{chr}(chr_bin)));
				end;
				colors(chr_bin,1) = c_post(1);
				colors(chr_bin,2) = c_post(2);
				colors(chr_bin,3) = c_post(3);
			end;
			% standard : end determine color of each bin.


			%% standard : draw colorbars.
			for chr_bin = 1:ceil(chr_size(chr)/bases_per_bin)
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
			% standard : end draw colorbars.

			%% standard : cgh plot section.
			c_ = [0 0 0];
			fprintf(['\nmain-plot : chr' num2str(chr) ':' num2str(length(CNVplot2{chr})) '\n']);
			fprintf(['ploidy     = ' num2str(ploidy)     '\n']);
			fprintf(['ploidyBase = ' num2str(ploidyBase) '\n']);
			for chr_bin = 1:ceil(chr_size(chr)/bases_per_bin)
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

			%% standard : draw lines across plots for easier interpretation of CNV regions.
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

			%% standard : axes labels etc.
			hold off;
			xlim([0,chr_size(chr)/bases_per_bin]);
    
			%% standard : modify y axis limits to show annotation locations if any are provided.
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

			%% standard : end axes labels etc.
			if (displayBREAKS == true) && (show_annotations == true)
				chr_length = ceil(chr_size(chr)/bases_per_bin);
				for segment = 2:length(chr_breaks{chr})-1
					bP = chr_breaks{chr}(segment)*chr_length;
					plot([bP bP], [(-maxY/10*2.5) 0],  'Color',[1 0 0],'LineWidth',2);
				end;
			end;   

			%% standard : show centromere outlines and horizontal marks.
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
    
			%% standard : show annotation locations
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

			%% standard : make CGH histograms to the right of the main chr cartoons.
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

			%% standard : places chr copy number to the right of the main chr cartoons.
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


			%% =========================================================================================
			% Draw angleplots to left of main chromosome cartoons.
			%-------------------------------------------------------------------------------------------
			apply_phasing = false;
			angle_plot_subfigures;


%%%%%%%%%%%%%%%%%%%%% END of standard figure draw section.


			%% Linear figure draw section
			if (Linear_display == true)
				figure(Linear_fig);
				Linear_width = Linear_Chr_max_width*chr_size(chr)/Linear_genome_size;
				subplot('Position',[Linear_left Linear_base Linear_width Linear_height]);
				Linear_left = Linear_left + Linear_width + Linear_left_chr_gap;
				hold on;
				title(chr_label{chr},'Interpreter','none','FontSize',20);

				%% linear : draw colorbars.
				for chr_bin = 1:ceil(chr_size(chr)/bases_per_bin)
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

				%% linear : cgh plot section.
				c_ = [0 0 0];
				fprintf(['linear-plot : chr' num2str(chr) ':' num2str(length(CNVplot2{chr})) '\n']);
				for chr_bin = 1:ceil(chr_size(chr)/bases_per_bin)
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

				%% linear : draw lines across plots for easier interpretation of CNV regions.
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

				%% linear : show segmental anueploidy breakpoints.
				if (Linear_displayBREAKS == true) && (show_annotations == true)
					chr_length = ceil(chr_size(chr)/bases_per_bin);
					for segment = 2:length(chr_breaks{chr})-1
						bP = chr_breaks{chr}(segment)*chr_length;
						plot([bP bP], [(-maxY/10*2.5) 0],  'Color',[1 0 0],'LineWidth',2);
					end;
				end;
				% linear : end segmental aneuploidy breakpoint section.

				%% linear : show centromere.
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

				%% linear : show annotation locations
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
