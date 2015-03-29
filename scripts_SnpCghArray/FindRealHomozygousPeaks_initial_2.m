function [realHomozygous_peak, disomy_fit, skew_factor] = ...
    FindRealHomozygousPeaks_initial_2(chrCopyNum,SNP_probeset_length,probeset1,chr_breaks,chr_size,show_unnassigned,DataTypeToUse,show_fitting, workingDir)
% FindRealHomozygousPeakLocation determines where homozygous peaks
%    are in the data for an experiment.
histAll = [];
no_usable_ploidies = 1;
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
disomy_fit(1) = 10;    disomy_fit(2) = 1;      disomy_fit(3) = 10;
disomy_fit(4) = 10;    disomy_fit(5) = 100.5;  disomy_fit(6) = 10;
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
	hist_a(101:200) = [];
	hist_b = fliplr(hist(histAll,200));
	hist_b(101:200) = [];
	hist_ab = hist_a+hist_b;
%	hist_ab(51:100) = [];

	smoothed_1a = smooth_gaussian(hist_ab,5,20);
	% find homozygous peak.
	histMax = find(smoothed_1a == max(smoothed_1a));
	% determines final homozygous peak location.
	initialHomozygous_peak = histMax/200;
    
	% generate a smoothed dataset of all outer quartile data, for use in gaussian fitting.
	smoother = [smoothed_1a fliplr(smoothed_1a)];
%%%	smoother(51:150) = 0;

%% Fit central Gaussian peak.
	initial_peak_position = 100;
	[fitting] = fit_Gaussian_model_heterozygousPeak_initial(smoother,initial_peak_position,'cubic',show_fitting,workingDir);

%% Delete fitted central Gaussians from smoothed data.
	smoother2 = smoother-fitting;
	smoother2(51:150) = 0;   % remove central peak residual after fitting removal.

% Extend end values of data profile before smoothing, to minimize edge effects.
	smoother2_int      = [ones(1,50)*smoother2(1) smoother2 ones(1,50)*smoother2(end)];
	smoother3          = smooth_gaussian(smoother2_int,2,20);
	smoother3(1:50)    = [];
	smoother3(201:250) = [];

%% Find homozygous peaks in cleaned up data.
	% get peaks in data profile.
	[peaks,locations] = findpeaks(smoother3(1:100));
test_peaks     = peaks
test_locations = locations
	% Find most extreme (left/right) peak.
	initialHomozygous_peak = locations(1)/200;
	% fit gaussians to homozygous peaks, to fine-tune their locations.
	disomy_peak_estimates = [initialHomozygous_peak*200 100.5 201-initialHomozygous_peak*200];
	skew_factor_initial = 75;

	[p1_a,p1_b,p1_c, p2_a,p2_b,p2_c, p3_a,p3_b,p3_c, skew_factor] = ...
	   fit_Gaussian_model_disomy_initial_2(smoother,disomy_peak_estimates,15,skew_factor_initial,'cubic',show_fitting, workingDir);
	realHomozygous_peak = p1_b/200;
	disomy_fit(1) = p1_a;    disomy_fit(2) = p1_b;    disomy_fit(3) = p1_c;
	disomy_fit(4) = p2_a;    disomy_fit(5) = p2_b;    disomy_fit(6) = p2_c;
	disomy_fit(7) = p3_a;    disomy_fit(8) = p3_b;    disomy_fit(9) = p3_c;

	%----------------------------------------------------------------------
	% Generate single Gaussian fitting figure.
	fig = figure(1234);
	plot(smoother2,'-', 'color',[0.2  0.8  0.2 ],'lineWidth',1);
	hold on;
		plot(smoother3,'-', 'color',[0.8  0.2  0.2 ],'lineWidth',1);
		plot([locations(1) locations(1)],[min(smoother2) max(smoother2)],'-','color',[0.0 0.0 0.0],'lineWidth',1);
		plot([201-locations(1) 201-locations(1)],[min(smoother2) max(smoother2)],'-','color',[0.0 0.0 0.0],'lineWidth',1);
	hold off;
	title('initial minus central peak');
	% saveas(fig, [workingDir 'initGaussianFit_minusFit.eps'], 'epsc');
	  saveas(fig, [workingDir 'initGaussianFit_minusFit.png'], 'png');
	delete(fig);
	%----------------------------------------------------------------------

end;

%%    smoother = [smoothed_1a 1:100 fliplr(smoothed_1a)];
%%    smoother(51:150) = 0;
%    % fit gaussians to homozygous peaks, to fine-tune their location.
%    disomy_peak_estimates = [initialHomozygous_peak*200 100 200-initialHomozygous_peak*200];
%    skew_factor_initial = 75;
%    [p1_a,p1_b,p1_c, p2_a,p2_b,p2_c, p3_a,p3_b,p3_c, skew_factor] = ...
%        fit_Gaussian_model_disomy_initial_2(smoother,disomy_peak_estimates,15,skew_factor_initial,'cubic',show_fitting, workingDir);
%    realHomozygous_peak = p1_b/200;
%    disomy_fit(1) = p1_a;    disomy_fit(2) = p1_b;    disomy_fit(3) = p1_c;
%    disomy_fit(4) = p2_a;    disomy_fit(5) = p2_b;    disomy_fit(6) = p2_c;
%    disomy_fit(7) = p3_a;    disomy_fit(8) = p3_b;    disomy_fit(9) = p3_c;

end
