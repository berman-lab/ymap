function [chr_breaks, chrCopyNum] = FindChrSizes( segmental_aneuploidy,CGH_probeset_length,probeset2,chr_size,flow_ploidy )
% FindChrSizes determines chromosome sizes from
%    initial ploidy estimate and CGH data.
maxY = 10;
chr_breaks = [];
chrCopyNum = [];
chrCopyNum1 = [];
chrCopyNum2 = [];
chrCopyNum_vector = [];
% precalculation of chromosome copy numbers.
for usedChr = [8 1:7]
	% determine where the endpoints of ploidy segments are.
	chr_breaks{usedChr}(1) = 0.0;
	break_count = 1;
	if (length(segmental_aneuploidy) > 0)
		for i = 1:length(segmental_aneuploidy)
			if (segmental_aneuploidy(i).chr == usedChr)
				break_count = break_count+1;
				chr_broken = true;
				chr_breaks{usedChr}(break_count) = segmental_aneuploidy(i).break;
			end;
		end;
	end;
	chr_breaks{usedChr}(length(chr_breaks{usedChr})+1) = 1;
	chr_breaks{usedChr}                                = sort(chr_breaks{usedChr});

	for segment = 1:length(chr_breaks{usedChr})-1
		smoothed = [];
		segment_CGHdata = [];
		% find set of CGH data for this segment of this chromosome.
		for i = 1:CGH_probeset_length
			% val = ploidy estimate adjusted copy number for each CGH probe.
			val = probeset2(i).probe_Ratio;

			% Collects CGH data for segment.
			if (length(val) > 0) && (probeset2(i).probe_chromosome == usedChr) ...
			&& (probeset2(i).probe_location >= 1+chr_size(usedChr)*chr_breaks{usedChr}(segment)) ...
			&& (probeset2(i).probe_location < chr_size(usedChr)*chr_breaks{usedChr}(segment+1))
				% collected CGH data for segment.
				segment_CGHdata = [segment_CGHdata val];
			end;
		end;
		if (flow_ploidy == 0)
			segment_CGHdata = segment_CGHdata*2;
		else
			segment_CGHdata = segment_CGHdata*flow_ploidy;
		end;

		% make smoothed histogram of CGH data for this segment.
		segment_CGHdata(segment_CGHdata==0)        = [];
		segment_CGHdata(length(segment_CGHdata)+1) = 0;   % endpoints added to ensure histogram bounds.
		segment_CGHdata(length(segment_CGHdata)+1) = maxY;
		segment_CGHdata(segment_CGHdata<0)         = [];
		segment_CGHdata(segment_CGHdata>maxY)      = [];
		histogram_width                            = 200;
		smoothed                                   = smooth_gaussian(hist(segment_CGHdata,histogram_width),5,20);

		% make smoothed histogram of endpoint data.
		end_points                                 = [0, maxY];
		end_points_smoothed                        = smooth_gaussian(hist(end_points,histogram_width),5,20);

		% remove smoothed endpoint data from smoothed histogram.
		smoothed2                                  = smoothed - end_points_smoothed;
		smoothed                                   = smoothed2;

		% find initial estimage of peak location from smoothed segment CGH data.
		peakLocation                               = find(smoothed==max(smoothed));

		% fit a single Gaussian to CNV data, to estimate copy number for each segment.
		show_fitting                               = 0;
% Troubleshooting outputs.
%		fprintf(['smoothed_data_length = [' num2str(length(smoothed)) ']\n']);
%		fprintf(['smoothed_data        = [' num2str(smoothed) ']\n']);
		[CGHsegment_height, CGHsegment_location, CGHsegment_width] = fit_Gaussian_model2(smoothed, peakLocation, 'cubic',show_fitting);

		% calculate copy number from Gaussian location.
		chrCopyNum{usedChr}(segment)               = round(CGHsegment_location/(histogram_width/maxY)*10)/10;
		chrCopyNum_vector                          = [chrCopyNum_vector chrCopyNum{usedChr}(segment)];
	end;
end;
% Adjustment of ploidy estimate.
% assumes most common copy number to really be a whole number. (2.1 -> 2, etc.)
common_copyNum = mode(chrCopyNum_vector);
for chr = 1:length(chrCopyNum)
	for segment = 1:length(chrCopyNum{chr})
		% rounds to 1 decimal place.
		chrCopyNum2{chr}(segment) = round(chrCopyNum{chr}(segment)/common_copyNum*round(common_copyNum)*10)/10;
		% rounds to quarters.
		%chrCopyNum2{chr}(segment) = floor(chrCopyNum2{chr}(segment)) + round((chrCopyNum2{chr}(segment)-floor(chrCopyNum2{chr}(segment)))/0.25)*0.25;
	end;
end;
chrCopyNum1 = chrCopyNum;
chrCopyNum  = chrCopyNum2;

end

