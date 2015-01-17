function [realHomozygous_peak, disomy_fit, skew_factor] = FindRealHomozygousPeaks_2(chrCopyNum,SNP_probeset_length,probeset1,chr_breaks,chr_size,show_unnassigned,DataTypeToUse,show_fitting, workingDir)

show_fitting = true;

% FindRealHomozygousPeakLocation determines where homozygous peaks
%    are in the data for an experiment.
no_usable_ploidies = 1;
histAll = [];
for usedChr = [8 1:7]
    % Checks if any whole chromosomes are found, as only these are simple
    % to analyze for homozygous peaks.
    for segment = 1:length(chrCopyNum{usedChr})
        if (chrCopyNum{usedChr}(segment) == 1) || (chrCopyNum{usedChr}(segment) == 2) || (chrCopyNum{usedChr}(segment) == 3) || (chrCopyNum{usedChr}(segment) == 4)
            no_usable_ploidies = 0;
            for i = 1:2:SNP_probeset_length
                if (probeset1(i).probe_location <= chr_breaks{usedChr}(segment+1)*chr_size(usedChr)) && (probeset1(i).probe_chromosome == usedChr)
                    if (length(probeset1(i).probe_Ratio) > 0) && (length(probeset1(i+1).probe_Ratio) > 0)
                        % Calculate value of SNP probe pair.
                        [UsedData] = calculateValue(probeset1,i,DataTypeToUse);
                        
                        if (isfield(probeset1(1),'probe_polarity') == 1)
                            if (probeset1(i).probe_polarity == 0)
                                if (show_unnassigned == true)
                                    histAll(i) = UsedData;
                                else
                                    histAll(i) = 0;
                                end;
                            elseif (probeset1(i).probe_polarity == 4)
                                % null action for when (probe_polarity == 4)
                                % due to probe design error; probes are identical.
                            else
                                histAll(i) = UsedData;
                            end;
                        else
                            histAll(i) = UsedData;
                        end;
                    end;
                end;
            end;
        end;
    end;
end;
realHomozygous_peak = 0;
disomy_fit(1) = 10;    disomy_fit(2) = 0;      disomy_fit(3) = 10;
disomy_fit(4) = 10;    disomy_fit(5) = 100;    disomy_fit(6) = 10;
disomy_fit(7) = 10;    disomy_fit(8) = 200;    disomy_fit(9) = 10;
if (no_usable_ploidies == 0) % disomics/tetrasomics were found, so we can determine homozygous peak location.
    % make a histogram of SNP allele ratio angles, then smooth it for display.
    histAll(histAll==0) = [];
    histAll(length(histAll)+1) = 0;
    histAll(length(histAll)+1) = 1;
    % Calculation of homozygous peak location.
    % Gather SNP data from this data set.
    histAll(histAll==0) = [];
    histAll(histAll>1) = [];
    histAll(length(histAll)+1) = 0;
    histAll(length(histAll)+1) = 1;
    % make a histogram of SNP allele ratio angles in theoretical
    % homozygous regions (outer quartiles only).
    hist_a = hist(histAll,200);
    hist_b = fliplr(hist(histAll,200));
    hist_ab = hist_a+hist_b;
    smoothed_1a = smooth_gaussian(hist_ab,5,20);

    [initialHomozygous_peak, disomy_fit, skew_factor] = FindRealHomozygousPeaks_initial_2(chrCopyNum,SNP_probeset_length,probeset1,chr_breaks,chr_size,show_unnassigned,DataTypeToUse,show_fitting, workingDir);

    % fit gaussians to homozygous peaks, to fine-tune their location.
    disomy_peak_estimates = [initialHomozygous_peak*200 100 200-initialHomozygous_peak*200];
    [p1_a,p1_b,p1_c, p2_a,p2_b,p2_c, p3_a,p3_b,p3_c, skew_factor] = ...
        fit_Gaussian_model_disomy_2(smoothed_1a,disomy_peak_estimates,15,0,skew_factor,'cubic',show_fitting);
    realHomozygous_peak = p1_b/200;
    disomy_fit(1) = p1_a;    disomy_fit(2) = p1_b;    disomy_fit(3) = p1_c;
    disomy_fit(4) = p2_a;    disomy_fit(5) = p2_b;    disomy_fit(6) = p2_c;
    disomy_fit(7) = p3_a;    disomy_fit(8) = p3_b;    disomy_fit(9) = p3_c;
else
	skew_factor = 75;
end;

end
