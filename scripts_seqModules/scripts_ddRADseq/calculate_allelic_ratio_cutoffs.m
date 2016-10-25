%% =========================================================================================
% Gather allelic ratio and coordinate data for SNPs into data structure.
%...........................................................................................
% if (useHapmap)
%       chr_SNPdata{chr,1}{chr_bin} = phased SNP ratio data.
%       chr_SNPdata{chr,2}{chr_bin} = unphased SNP ratio data.
%       chr_SNPdata{chr,3}{chr_bin} = phased SNP position data.
%       chr_SNPdata{chr,4}{chr_bin} = unphased SNP position data.
%       chr_SNPdata{chr,5}{chr_bin} = flipper value for phased SNP.
%       chr_SNPdata{chr,6}{chr_bin} = flipper value for unphased SNP.
% elseif (useParent)
%       chr_SNPdata{chr,1}{chr_bin} = parent SNP ratio data.
%       chr_SNPdata{chr,2}{chr_bin} = child SNP ratio data.
%       chr_SNPdata{chr,3}{chr_bin} = parent SNP position data.
%       chr_SNPdata{chr,4}{chr_bin} = child SNP position data.
% else
%       chr_SNPdata{chr,1}{chr_bin} = SNP ratio data.
%       chr_SNPdata{chr,3}{chr_bin} = SNP position data.
% end;
%-------------------------------------------------------------------------------------------
if (useHapmap)
	%       chr_SNPdata{chr,1}{chr_bin} = phased SNP ratio data.
	%       chr_SNPdata{chr,2}{chr_bin} = unphased SNP ratio data.
	%       chr_SNPdata{chr,3}{chr_bin} = phased SNP position data.
	%       chr_SNPdata{chr,4}{chr_bin} = unphased SNP position data.
	%       chr_SNPdata{chr,5}{chr_bin} = flipper value for phased SNP.
	%       chr_SNPdata{chr,6}{chr_bin} = flipper value for unphased SNP.
	for chr = 1:num_chrs
		if (chr_in_use(chr) == 1)
			if (length(C_chr_SNP_data_positions{chr}) > 1)
				%
				% Building SNP data sctructure 'chr_SNPdata'.
				%

				for SNP = 1:length(C_chr_SNP_data_positions{chr})  % 'length(C_chr_SNP_data_positions)' is the number of SNPs per chromosome.
					coordinate                      = C_chr_SNP_data_positions{chr}(SNP);
					pos                             = ceil(coordinate/new_bases_per_bin);
					localCopyEstimate               = round(CNVplot2{chr}(pos)*ploidy*ploidyAdjust);
					baseCall                        = C_chr_baseCall{          chr}{SNP};
					homologA                        = C_chr_SNP_homologA{      chr}{SNP};
					homologB                        = C_chr_SNP_homologB{      chr}{SNP};
					flipper                         = C_chr_SNP_flipHomologs{  chr}(SNP);
					allelic_ratio                   = C_chr_SNP_data_ratios{   chr}(SNP);
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
						chr_SNPdata{chr,1}{pos} = [chr_SNPdata{chr,1}{pos} allelic_ratio allelic_ratio];
						chr_SNPdata{chr,3}{pos} = [chr_SNPdata{chr,3}{pos} coordinate    coordinate   ];
						chr_SNPdata{chr,5}{pos} = [chr_SNPdata{chr,5}{pos} flipper       flipper      ];
					else % (flipper == 0)
						if (baseCall == homologA)
							allelic_ratio = 1-allelic_ratio;
						end;
						chr_SNPdata{chr,1}{pos} = [chr_SNPdata{chr,1}{pos} allelic_ratio allelic_ratio];
						chr_SNPdata{chr,3}{pos} = [chr_SNPdata{chr,3}{pos} coordinate    coordinate   ];
						chr_SNPdata{chr,5}{pos} = [chr_SNPdata{chr,5}{pos} flipper       flipper      ];
					end;
				end;
			end;
		end;
	end;
elseif (useParent)
	%       chr_SNPdata{chr,1}{chr_bin} = child SNP ratio data.
	%       chr_SNPdata{chr,2}{chr_bin} = parent SNP ratio data.
	%       chr_SNPdata{chr,3}{chr_bin} = child SNP position data.
	%       chr_SNPdata{chr,4}{chr_bin} = parent SNP position data.
	%
	% child data:  'C_chr_SNP_data_positions','C_chr_SNP_data_ratios','C_chr_count'
        % parent data: 'P_chr_SNP_data_positions','P_chr_SNP_data_ratios','P_chr_count'
	%
	for chr = 1:num_chrs
		if (chr_in_use(chr) == 1)
			if (length(C_chr_SNP_data_positions{chr}) > 1)
				%
				% Building SNP data sctructure 'chr_SNPdata' entries for child data.
				%

				for SNP = 1:length(C_chr_SNP_data_positions{chr})
					coordinate              = C_chr_SNP_data_positions{chr}(SNP);
					pos                     = ceil(coordinate/new_bases_per_bin);
					allelic_ratio           = C_chr_SNP_data_ratios{   chr}(SNP);
					allelic_ratio           = max(allelic_ratio, 1-allelic_ratio);

					% Allelic ratio here is the ratio of the majority read call to all reads.
					% The consequence of this is that it will always be on the range [0.5 .. 1.0].
					chr_SNPdata{chr,1}{pos} = [chr_SNPdata{chr,1}{pos} allelic_ratio];
					chr_SNPdata{chr,3}{pos} = [chr_SNPdata{chr,3}{pos} coordinate   ];
				end;
			end;
			if (length(P_chr_SNP_data_positions{chr}) > 1)
				%
				% Building SNP data sctructure 'chr_SNPdata' entries for parent data.
				%

				for SNP = 1:length(P_chr_SNP_data_positions{chr})
					coordinate              = P_chr_SNP_data_positions{chr}(SNP);
					pos                     = ceil(coordinate/new_bases_per_bin);
					allelic_ratio           = P_chr_SNP_data_ratios{   chr}(SNP);
					allelic_ratio           = max(allelic_ratio, 1-allelic_ratio);

					% Allelic ratio here is the ratio of the majority read call to all reads.
					% The consequence of this is that it will always be on the range [0.5 .. 1.0].
					chr_SNPdata{chr,2}{pos} = [chr_SNPdata{chr,2}{pos} allelic_ratio];
					chr_SNPdata{chr,4}{pos} = [chr_SNPdata{chr,4}{pos} coordinate   ];
				end;
			end;
		end;
	end;
else    % useHapmap and useParent == false.
	% child == parent, so only one set of data needs to be organized.
	%       chr_SNPdata{chr,1}{chr_bin} = SNP ratio data.
	%       chr_SNPdata{chr,3}{chr_bin} = SNP position data.
	for chr = 1:num_chrs
		if (chr_in_use(chr) == 1)
			if (length(C_chr_SNP_data_positions{chr}) > 1)
				%
				% Building SNP data sctructure 'chr_SNPdata'.
				%

				for SNP = 1:length(C_chr_SNP_data_positions{chr})
					coordinate              = C_chr_SNP_data_positions{chr}(SNP);
					pos                     = ceil(coordinate/new_bases_per_bin);
					allelic_ratio           = C_chr_SNP_data_ratios{   chr}(SNP);
					allelic_ratio           = max(allelic_ratio, 1-allelic_ratio);

					% Allelic ratio here is the ratio of the majority read call to all reads.
					% The consequence of this is that it will always be on the range [0.5 .. 1.0].
					chr_SNPdata{chr,1}{pos} = [chr_SNPdata{chr,1}{pos} allelic_ratio];
					chr_SNPdata{chr,3}{pos} = [chr_SNPdata{chr,3}{pos} coordinate   ];
				end;
			end;
		end;
	end;
end;


%% =========================================================================================
% Calculate allelic fraction cutoffs.
%-------------------------------------------------------------------------------------------
% Initialize 
for chr = num_chrs
	if (chr_in_use(chr) == 1)
		for segment = 1:length(chrCopyNum{chr})
			chrSegment_peaks{              chr}{segment} = [];
			chrSegment_mostLikelyGaussians{chr}{segment} = [];
			chrSegment_actual_cutoffs{     chr}{segment} = [];
			chrSegment_smoothed{           chr}{segment} = [];
		end;
	end;
end;
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		chr_length = chr_size(chr);
		for segment = 1:length(chrCopyNum{chr})
			histAll_a = [];
			histAll_b = [];
			histAll2  = [];
			% Look through all SNP data in every chr_bin, to determine if any are within the segment boundries.
			% Speed up by only checking possible chr_bins has not been implmented.

			fprintf( '^^^\n');
			fprintf(['^^^ chrID         = ' num2str(chr)                                   '\n']);
			fprintf(['^^^ segmentID     = ' num2str(segment)                               '\n']);
			fprintf(['^^^ segment start = ' num2str(chr_breaks{chr}(segment  )*chr_length) '\n']);
			fprintf(['^^^ segment end   = ' num2str(chr_breaks{chr}(segment+1)*chr_length) '\n']);

			%% Construct and smooth a histogram of alleleic fraction data in the segment of interest.
			% phased data is stored into arrays 'histAll_a' and 'histAll_b', since proper phasing is known.
			% unphased data is stored inverted into the second array, since proper phasing is not known.
			for chr_bin = 1:length(CNVplot2{chr})
				%   1 : phased SNP ratio data.
				%   2 : unphased SNP ratio data.
				%   3 : phased SNP position data.
				%   4 : unphased SNP position data.
				ratioData_phased        = chr_SNPdata{chr,1}{chr_bin};
				ratioData_unphased      = chr_SNPdata{chr,2}{chr_bin};
				coordinateData_phased   = chr_SNPdata{chr,3}{chr_bin};
				coordinateData_unphased = chr_SNPdata{chr,4}{chr_bin};
				if (useHapmap)
					if (length(ratioData_phased) > 0)
						for SNP_in_bin = 1:length(ratioData_phased)
							if ((coordinateData_phased(SNP_in_bin) > chr_breaks{chr}(segment)*chr_length) && (coordinateData_phased(SNP_in_bin) <= chr_breaks{chr}(segment+1)*chr_length))
								% Ratio data is phased, so it is added twice in its proper orientation (to match density of unphased data below).
								allelic_ratio = ratioData_phased(SNP_in_bin);
								histAll_a = [histAll_a allelic_ratio  ];
								histAll_b = [histAll_b allelic_ratio  ];
							end;
						end;
					end;
				end;
				if (length(ratioData_unphased) > 0)
					for SNP_in_bin = 1:length(ratioData_unphased)
						if ((coordinateData_unphased(SNP_in_bin) > chr_breaks{chr}(segment)*chr_length) && (coordinateData_unphased(SNP_in_bin) <= chr_breaks{chr}(segment+1)*chr_length))
							% Ratio data is unphased, so it is added evenly in both orientations.
							allelic_ratio = ratioData_unphased(SNP_in_bin);
							histAll_a = [histAll_a allelic_ratio  ];
							histAll_b = [histAll_b 1-allelic_ratio];
						end;
					end;
				end;
			end;

			% make a histogram of SNP allelic fractions in segment, then smooth for display.
			histAll                    = [histAll_a histAll_b];
			histAll(histAll == -1)     = [];

			% Invert histogram values;
			histAll                    = 1-histAll;

			% add bounds to the histogram values.
			histAll                    = [histAll 0 1];

			% generate the histogram.
			data_hist                  = hist(histAll,200);
			endPoints_hist             = hist([0 1],200);

			% remove the endpoints.
			data_hist                  = data_hist-endPoints_hist;

			% log-scale the histogram to minimize difference between hom & het peak heights.
			% average this with the raw histogram so the large peaks still appear visibily larger than the small peaks.
			
			%data_hist                  = (data_hist + log(data_hist+1))/2;
			data_hist                  = log(data_hist+1);
			

			% smooth the histogram.
			smoothed                   = smooth_gaussian(data_hist,10,30);

			% flip the smoothed histogram left-right to make display consistent with values.
			smoothed                   = fliplr(smoothed);

			% scale the smoothed histogram to a max of 1.
			if (max(smoothed) > 0)
				smoothed           = smoothed/max(smoothed);
			end;

			%% Calculate Gaussian fitting details for segment.
			segment_copyNum            = round(chrCopyNum{chr}(segment));  % copy number estimate of this segment.
			segment_chrBreaks          = chr_breaks{chr}(segment);         % break points of this segment.
			segment_smoothedHistogram  = smoothed;                         % whole chromosome allelic ratio histogram smoothed.

			% Define cutoffs between Gaussian fits.
			saveName = ['allelic_ratios.chr_' num2str(chr) '.seg_' num2str(segment)];
			[peaks,actual_cutoffs,mostLikelyGaussians] = FindGaussianCutoffs_3(workingDir,saveName, chr,segment, segment_copyNum,segment_smoothedHistogram, false);

			fprintf(['^^^ copyNum             = ' num2str(segment_copyNum    ) '\n']);
			fprintf(['^^^ peaks               = ' num2str(peaks              ) '\n']);
			fprintf(['^^^ mostLikelyGaussians = ' num2str(mostLikelyGaussians) '\n']);
			fprintf(['^^^ actual_cutoffs      = ' num2str(actual_cutoffs     ) '\n']);

			chrSegment_peaks{              chr}{segment} = peaks;
			chrSegment_mostLikelyGaussians{chr}{segment} = mostLikelyGaussians;
			chrSegment_actual_cutoffs{     chr}{segment} = actual_cutoffs;
			chrSegment_smoothed{           chr}{segment} = smoothed;
		end;
	end;
end;
