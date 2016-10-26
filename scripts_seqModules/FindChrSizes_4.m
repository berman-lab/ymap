function [chr_breaks, chrCopyNum, ploidyAdjust] = FindChrSizes_4(Aneuploidy,CNVplot,Ploidy,num_chrs,chr_in_use)

% FindChrSizes determines chromosome sizes from
%    initial ploidy estimate and CGH data.
maxY              = 10;
chr_breaks        = [];
chrCopyNum        = [];
chrCopyNum1       = [];
chrCopyNum2       = [];
chrCopyNum_vector = [];
% precalculation of chromosome copy numbers.
for usedChr = 1:num_chrs
	if (chr_in_use(usedChr) == 1)
		% determine where the endpoints of ploidy segments are.
		chr_breaks{usedChr}(1) = 0.0;
		break_count = 1;

		if (length(Aneuploidy) > 0)
			for i = 1:length(Aneuploidy)
				%if (Aneuploidy(i).dataset == dataset) && (Aneuploidy(i).chr == usedChr)
				if (Aneuploidy(i).chr == usedChr)
					break_count = break_count+1;
					chr_broken = true;
					chr_breaks{usedChr}(break_count) = Aneuploidy(i).break;
				end;
			end;
		end;
		chr_breaks{usedChr}(length(chr_breaks{usedChr})+1) = 1;
    
		fprintf(['chr' num2str(usedChr) ' : ' num2str(length(chr_breaks{usedChr})) '\n']);
		for segment = 1:length(chr_breaks{usedChr})-1
			smoothed = [];
			smoothed2 = [];
			segment_CGHdata = [];
			segment_CGHdata2  = [];
			% find set of CGH data for this segment of this chromosome.
			for i = 1:length(CNVplot{usedChr})
				% val = ploidy estimate adjusted copy number for each CGH probe.
				if (i <= length(CNVplot{usedChr})*chr_breaks{usedChr}(segment+1)) && ...
				   (i >= length(CNVplot{usedChr})*chr_breaks{usedChr}(segment))
					val = CNVplot{usedChr}(i);
					segment_CGHdata = [segment_CGHdata val];
				end;
			end;
			if (Ploidy == 0)
				segment_CGHdata = segment_CGHdata*2;
			else
				segment_CGHdata = segment_CGHdata*Ploidy;
			end;
			% make smoothed histogram of CGH data for this segment.
			segment_CGHdata(segment_CGHdata==0) = [];
			segment_CGHdata(length(segment_CGHdata)+1) = 0;   % endpoints added to ensure histogram bounds.
			segment_CGHdata(length(segment_CGHdata)+1) = maxY;
			% clearing 
			segment_CGHdata(segment_CGHdata<0) = [];
			segment_CGHdata(segment_CGHdata>maxY) = [];
			histogram_width = 200;
			smoothed        = smooth_gaussian(hist(segment_CGHdata,histogram_width),5,20);
			
			% make a smoothed version of just the endpoints used to ensure histogram bounds.
			segment_CGHdata2(1) = 0;
			segment_CGHdata2(2) = maxY;
			smoothed2 = smooth_gaussian(hist(segment_CGHdata2,histogram_width),5,20);

			% subtract the smoothed endpoints from the histogram to remove the influence of the added endpoints.
			smoothed = smoothed - smoothed2;
			smoothed = smoothed/max(smoothed);
            
			% find initial estimage of peak location from smoothed segment CGH data.
			peakLocation = find(smoothed==max(smoothed));
			% fit Gaussian to segment CGH data.
			show_fitting = 0;
			if (peakLocation > 1)
				[CGHsegment_height, CGHsegment_location, CGHsegment_width] = fit_Gaussian_model2(smoothed, peakLocation, 'cubic',show_fitting,20);
				% calculate copy number from Gaussian location.
				chrCopyNum{usedChr}(segment) = round(CGHsegment_location/(histogram_width/maxY)*10)/10;
			else
				chrCopyNum{usedChr}(segment) = 0;
			end;
			chrCopyNum_vector = [chrCopyNum_vector chrCopyNum{usedChr}(segment)];
		end;
	end;
end;

% Adjustment of ploidy estimate.
% assumes most common copy number to really be a whole number. (2.1 -> 2, etc.)
% zero data indicate erroneous copy number estimates and are first excluded.
chrCopyNum_vector(chrCopyNum_vector == 0) = [];
common_copyNum = mode(chrCopyNum_vector);
for chr = 1:length(chrCopyNum)
	for segment = 1:length(chrCopyNum{chr})
		% avoid avoid dividing by Nan if common_copyNum is NaN (since the whole copy vector can be empty)
		if (~isnan(common_copyNum))
			% rounds to 1 decimal place.
			chrCopyNum2{chr}(segment) = round(chrCopyNum{chr}(segment)/common_copyNum*round(common_copyNum)*10)/10;
		else
			chrCopyNum2{chr}(segment) = 0;
		end;
	end;
end;
chrCopyNum1  = chrCopyNum;
chrCopyNum   = chrCopyNum2;
% set ploidy adjust according to copy num, avoid dividing by Nan if common_copyNum is NaN (since the whole copy vector can be empty)
if (~isnan(common_copyNum))
    ploidyAdjust = round(common_copyNum)/common_copyNum;
else
    ploidyAdjust = Ploidy;
end;

end

