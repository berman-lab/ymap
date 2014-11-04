function [raw,smoothed,x_peak,actual_cutoffs,mostLikelyGaussians,chrCopyNum] = FindGaussianCutoffs_2(probeset1,chrCopyNum,chr_breaks,chr_size,chromosome,segment, ...
    monosomy_peak,disomy_peak,trisomy_peak,tetrasomy_peak,pentasomy_peak,hexasomy_peak,skew_factor,name,file_dir,MakeFigure,show_fitting,DataTypeToUse)

% Overlap of 1st SD before assume identity.
OverlapLimit = 0.5;

% initialize vectors for scattergram analysis
data1_0 = [];   data2_0 = [];   data1_1 = [];
data2_1 = [];   data1_2 = [];   data2_2 = [];

% Gathers SNP probe data for this chromosome segment.
for i = 1:2:length(probeset1)
    if (probeset1(i).probe_location > chr_breaks{chromosome}(segment)*chr_size(chromosome)) && (probeset1(i).probe_location <= chr_breaks{chromosome}(segment+1)*chr_size(chromosome))
        % determines if probe pair is useful; both probes have data.
        if (length(probeset1(i).probe_Ratio) > 0) && (length(probeset1(i+1).probe_Ratio) > 0)
            % swaps probes if polarity switch indicates to do so.
            if (probeset1(i).probe_chromosome == chromosome)
                val1 = probeset1(i).probe_Ratio;
                val2 = probeset1(i+1).probe_Ratio;
                if (probeset1(i).probe_polarity == 1)   % [no-swap]
                    data1_1(length(data1_1)+1) = val1;
                    data2_1(length(data2_1)+1) = val2;
                elseif (probeset1(i).probe_polarity == 2)   % [swap]
                    data1_2(length(data1_2)+1) = val2;
                    data2_2(length(data2_2)+1) = val1;
                elseif (probeset1(i).probe_polarity == 0) %unassigned probe pairs.   [no-swap]
                    data1_0(length(data1_0)+1) = val1;
                    data2_0(length(data2_0)+1) = val2;
                else
                    % null action for when (probe_polarity == 4)
                    % due to probe design error; probes are identical.
                end;
            end;
        end;
    end;
end;

% Calculate allelic_fractions/angles for this chromosome segment.
histAll_1a = [];
if (DataTypeToUse == 1)   % use AllelicFractions.
    for i = 1:length(data1_0)
        histAll_1a = [histAll_1a data2_0(i)/(data2_0(i)+data1_0(i))];
    end;
    for i = 1:length(data1_1)
        histAll_1a = [histAll_1a data2_1(i)/(data2_1(i)+data1_1(i))];
    end;
    for i = 1:length(data1_2)
        histAll_1a = [histAll_1a data2_2(i)/(data2_2(i)+data1_2(i))];
    end;
else % (DataTypeToUse == 2)   % use Angles.
    for i = 1:length(data1_0)
        histAll_1a = [histAll_1a atan2(data2_0(i),data1_0(i))/(pi/2)];
    end;
    for i = 1:length(data1_1)
        histAll_1a = [histAll_1a atan2(data2_1(i),data1_1(i))/(pi/2)];
    end;
    for i = 1:length(data1_2)
        histAll_1a = [histAll_1a atan2(data2_2(i),data1_2(i))/(pi/2)];
    end;
end;

% make a histogram of SNP alleleic fractions, then smooth it for display.
histAll_1a(histAll_1a==0) = [];
histAll_1a(histAll_1a>1)  = [];
histAll_1a                = [histAll_1a 0 1];

% smoothed version and raw version for later use.
smoothed1 = smooth_gaussian(hist(histAll_1a,200),3,20);
smoothed2 = smooth_gaussian(hist([0 1],200),3,20);
smoothed = smoothed1-smoothed2;
raw      = hist(histAll_1a,200);

%% Calculation of Gaussians against per chromosome data.
% Fits Gaussians to real data per chromomsome, then
% determines equal probability cutoffs between them.
sigma = 5;

%% FindGaussianCutoffs Finds cutoffs as intersections of Gaussians, fit to
%    the data at each peak location.
ErrorType      = 'linear';
if (MakeFigure == true)
    fig = figure(1);
    hold on;
    [~,~,~,~,~,~,~,~,~,~,~,~,~,colorPeak,colorCutoff] = DefineColors();
end;

if (chrCopyNum{chromosome}(segment) <= 1)
    fraction = 0;
    peaks = [];
    peaks(1) = monosomy_peak(1);
    peaks(2) = monosomy_peak(2);
    G = [];
    [G{1}.a,G{1}.b,G{1}.c, G{2}.a,G{2}.b,G{2}.c, skew_factor] = ...
        fit_Gaussian_model_monosomy_2(smoothed,peaks*200,sigma,fraction,skew_factor,ErrorType,show_fitting);
    G{1}.d = skew_factor;
    G{2}.d = skew_factor;
    [list] = FindHighestGaussian_2(G);

    actual_cutoffs = [];
    mostLikelyGaussians = [];
    for i = 1:199
        if (list(i) ~= list(i+1))   % we've found a boundary.
            actual_cutoffs = [actual_cutoffs FindGaussianCrossover_2(G{list(i)},G{list(i+1)},i)];
            mostLikelyGaussians = [mostLikelyGaussians list(i)];
        end;
    end;
    mostLikelyGaussians = [mostLikelyGaussians list(200)];

    if (MakeFigure == true)
        figure(fig);
        hold on;
        time1_1 =   1:floor(G{1}.b);
        time1_2 = ceil(G{1}.b):200;
        if (time1_1(end) == time1_2(1));    time1_1(end) = [];  end;
        time2_1 =   1:floor(G{2}.b);
        time2_2 = ceil(G{2}.b):200;
        if (time2_1(end) == time2_2(1));    time2_2(1) = [];    end;
        fit_curve_1 = [G{1}.a*exp(-0.5*((time1_1-G{1}.b)./G{1}.c).^2) ...
                       G{1}.a*exp(-0.5*((time1_2-G{1}.b)./G{1}.c/(skew_factor/G{1}.b)).^2)];
        fit_curve_2 = [G{2}.a*exp(-0.5*((time2_1-G{2}.b)./G{2}.c/(skew_factor/(200-G{2}.b))).^2) ...
                       G{2}.a*exp(-0.5*((time2_2-G{2}.b)./G{2}.c).^2)];

        area(           max([1 round(G{1}.b-G{1}.c)]):min([200 round(G{1}.b+G{1}.c*(skew_factor/G{1}.b))]), ...
            fit_curve_1(max([1 round(G{1}.b-G{1}.c)]):min([200 round(G{1}.b+G{1}.c*(skew_factor/G{1}.b))])),'FaceColor',[0.5000 1.0 0.5000]);
        area(           max([1 round(G{2}.b-G{2}.c*(skew_factor/(200-G{2}.b)))]):min([200 round(G{2}.b+G{2}.c)]), ...
            fit_curve_2(max([1 round(G{2}.b-G{2}.c*(skew_factor/(200-G{2}.b)))]):min([200 round(G{2}.b+G{2}.c)])),'FaceColor',[0.5000 1.0 0.5000]);
        plot(fit_curve_1 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
        plot(fit_curve_2 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
        fit_curve_tot = fit_curve_1+fit_curve_2;
        plot(fit_curve_tot ,'color',[0.0 1.0 1.0],'linestyle','-','linewidth',2);
        hold off;
    end;

	x_peak = [];
    x_peak(1) = G{1}.b;
    x_peak(2) = G{2}.b;

    if (MakeFigure == true)
        hold on;
        for i = 1:2
            plot([x_peak(i) x_peak(i)],[0 max(raw)],'color',colorPeak,'linestyle','-','linewidth',2);
        end;
        for i = 1:length(actual_cutoffs)
            plot([actual_cutoffs(i) actual_cutoffs(i)],[0 max(raw)],'color',colorCutoff,'linestyle','-','linewidth',2);
        end;
        hold off;
    end;
elseif (chrCopyNum{chromosome}(segment) <= 2)
    fraction = 2-chrCopyNum{chromosome}(segment);
    peaks = [];
    peaks(1) = disomy_peak(1);
    peaks(2) = disomy_peak(2);
    peaks(3) = disomy_peak(3);
    G = [];
    [G{1}.a,G{1}.b,G{1}.c, G{2}.a,G{2}.b,G{2}.c, G{3}.a,G{3}.b,G{3}.c, skew_factor] = ...
        fit_Gaussian_model_disomy_2(smoothed,peaks*200,sigma,fraction,skew_factor,ErrorType,show_fitting);
    G{1}.d = skew_factor;
    G{2}.d = skew_factor;
    G{3}.d = skew_factor;
    [list] = FindHighestGaussian_2(G);

	actual_cutoffs = [];
    mostLikelyGaussians = [];
    for i = 1:199
        if (list(i) ~= list(i+1))   % we've found a boundary.
            actual_cutoffs = [actual_cutoffs FindGaussianCrossover_2(G{list(i)},G{list(i+1)},i)];
            mostLikelyGaussians = [mostLikelyGaussians list(i)];
        end;
    end;
    mostLikelyGaussians = [mostLikelyGaussians list(200)];

    if (MakeFigure == true)
        figure(fig);
        hold on;
        time1_1 = 1:floor(G{1}.b);
        time1_2 = ceil(G{1}.b):200;
        if (time1_1(end) == time1_2(1));    time1_1(end) = [];  end;
        time2   = 1:200;
        time3_1 = 1:floor(G{3}.b);
        time3_2 = ceil(G{3}.b):200;
        if (time3_1(end) == time3_2(1));    time3_2(1) = [];    end;

        fit_curve_1 = [G{1}.a*exp(-0.5*((time1_1-G{1}.b)./G{1}.c).^2) ...
                       G{1}.a*exp(-0.5*((time1_2-G{1}.b)./G{1}.c/(skew_factor/G{1}.b)).^2)];
        fit_curve_2 =  G{2}.a*exp(-0.5*((time2-G{2}.b)./G{2}.c).^2);
        fit_curve_3 = [G{3}.a*exp(-0.5*((time3_1-G{3}.b)./G{3}.c/(skew_factor/(200-G{3}.b))).^2) ...
                       G{3}.a*exp(-0.5*((time3_2-G{3}.b)./G{3}.c).^2)];

		area(           max([1 round(G{1}.b-G{1}.c)]):min([200 round(G{1}.b+G{1}.c*(skew_factor/G{1}.b))]), ...
            fit_curve_1(max([1 round(G{1}.b-G{1}.c)]):min([200 round(G{1}.b+G{1}.c*(skew_factor/G{1}.b))])),'FaceColor',[0.5000 1.0 0.5000]);
        area(           max([1 round(G{2}.b-G{2}.c)]):min([200 round(G{2}.b+G{2}.c)]), ...
            fit_curve_2(max([1 round(G{2}.b-G{2}.c)]):min([200 round(G{2}.b+G{2}.c)])),'FaceColor',[0.5000 1.0 0.5000]);
        area(           max([1 round(G{3}.b-G{3}.c*(skew_factor/(200-G{3}.b)))]):min([200 round(G{3}.b+G{3}.c)]), ...
            fit_curve_3(max([1 round(G{3}.b-G{3}.c*(skew_factor/(200-G{3}.b)))]):min([200 round(G{3}.b+G{3}.c)])),'FaceColor',[0.5000 1.0 0.5000]);
        plot(fit_curve_1 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
        plot(fit_curve_2 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
        plot(fit_curve_3 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
        fit_curve_tot = fit_curve_1+fit_curve_2+fit_curve_3;
        plot(fit_curve_tot ,'color',[0.0 1.0 1.0],'linestyle','-','linewidth',2);
        hold off;
    end;

    x_peak = [];
    x_peak(1) = G{1}.b;
    x_peak(2) = G{2}.b;
    x_peak(3) = G{3}.b;

    if (MakeFigure == true)
        hold on;
        for i = 1:3
            plot([x_peak(i) x_peak(i)],[0 max(raw)],'color',colorPeak,'linestyle','-','linewidth',2);
        end;
        for i = 1:length(actual_cutoffs)
            plot([actual_cutoffs(i) actual_cutoffs(i)],[0 max(raw)],'color',colorCutoff,'linestyle','-','linewidth',2);
        end;
        hold off;
    end;
elseif (chrCopyNum{chromosome}(segment) <= 3)
	fraction = 3-chrCopyNum{chromosome}(segment);
    peaks = [];
    peaks(1) = trisomy_peak(1);
    peaks(2) = trisomy_peak(2)*(1-fraction) + disomy_peak(2)*fraction;
    peaks(3) = trisomy_peak(3)*(1-fraction) + disomy_peak(2)*fraction;
    peaks(4) = trisomy_peak(4);
    G = [];
    [G{1}.a,G{1}.b,G{1}.c, G{2}.a,G{2}.b,G{2}.c, G{3}.a,G{3}.b,G{3}.c, G{4}.a,G{4}.b,G{4}.c, skew_factor] = ...
        fit_Gaussian_model_trisomy_2(smoothed,peaks*200,sigma,fraction,skew_factor,ErrorType,show_fitting);
    G{1}.d = skew_factor;
    G{2}.d = skew_factor;
    G{3}.d = skew_factor;
    G{4}.d = skew_factor;
    [list] = FindHighestGaussian_2(G);

    actual_cutoffs = [];
    mostLikelyGaussians = [];
    for i = 1:199
        if (list(i) ~= list(i+1))   % we've found a boundary.
            actual_cutoffs = [actual_cutoffs FindGaussianCrossover_2(G{list(i)},G{list(i+1)},i)];
            mostLikelyGaussians = [mostLikelyGaussians list(i)];
        end;
    end;
    mostLikelyGaussians = [mostLikelyGaussians list(200)];

    % If central two gaussians overlap by their 1st SD, then assume one
    % central peak and recalculate.   The difference from 2n is not
    % significant.
    if (G{2}.b+G{2}.c*OverlapLimit > G{3}.b-G{3}.c*OverlapLimit)   % recalculate as if disomy.
        chrCopyNum{chromosome}(segment) = 2;
        fraction = 0;
        peaks = [];
        peaks(1) = disomy_peak(1);
		peaks(2) = disomy_peak(2);
        peaks(3) = disomy_peak(3);
        G = [];
        [G{1}.a,G{1}.b,G{1}.c, G{2}.a,G{2}.b,G{2}.c, G{3}.a,G{3}.b,G{3}.c, skew_factor] = ...
            fit_Gaussian_model_disomy_2(smoothed,peaks*200,sigma,fraction,skew_factor,ErrorType,show_fitting);
        G{1}.d = skew_factor;
        G{2}.d = skew_factor;
        G{3}.d = skew_factor;
        [list] = FindHighestGaussian_2(G);

        actual_cutoffs = [];
        mostLikelyGaussians = [];
        for i = 1:199
            if (list(i) ~= list(i+1))   % we've found a boundary.
                actual_cutoffs = [actual_cutoffs FindGaussianCrossover_2(G{list(i)},G{list(i+1)},i)];
                mostLikelyGaussians = [mostLikelyGaussians list(i)];
            end;
        end;
        mostLikelyGaussians = [mostLikelyGaussians list(200)];

        if (MakeFigure == true)
            figure(fig);
            hold on;
            time1_1 = 1:floor(G{1}.b);
            time1_2 = ceil(G{1}.b):200;
            if (time1_1(end) == time1_2(1));    time1_1(end) = [];  end;
			time2   = 1:200;
            time3_1 = 1:floor(G{3}.b);
            time3_2 = ceil(G{3}.b):200;
            if (time3_1(end) == time3_2(1));    time3_2(1) = [];    end;

            fit_curve_1 = [G{1}.a*exp(-0.5*((time1_1-G{1}.b)./G{1}.c).^2) ...
                           G{1}.a*exp(-0.5*((time1_2-G{1}.b)./G{1}.c/(skew_factor/G{1}.b)).^2)];
            fit_curve_2 =  G{2}.a*exp(-0.5*((time2-G{2}.b)./G{2}.c).^2);
            fit_curve_3 = [G{3}.a*exp(-0.5*((time3_1-G{3}.b)./G{3}.c/(skew_factor/(200-G{3}.b))).^2) ...
                           G{3}.a*exp(-0.5*((time3_2-G{3}.b)./G{3}.c).^2)];

            area(           max([1 round(G{1}.b-G{1}.c)]):min([200 round(G{1}.b+G{1}.c*(skew_factor/G{1}.b))]), ...
                fit_curve_1(max([1 round(G{1}.b-G{1}.c)]):min([200 round(G{1}.b+G{1}.c*(skew_factor/G{1}.b))])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{2}.b-G{2}.c)]):min([200 round(G{2}.b+G{2}.c)]), ...
                fit_curve_2(max([1 round(G{2}.b-G{2}.c)]):min([200 round(G{2}.b+G{2}.c)])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{3}.b-G{3}.c*(skew_factor/(200-G{3}.b)))]):min([200 round(G{3}.b+G{3}.c)]), ...
                fit_curve_3(max([1 round(G{3}.b-G{3}.c*(skew_factor/(200-G{3}.b)))]):min([200 round(G{3}.b+G{3}.c)])),'FaceColor',[0.5000 1.0 0.5000]);
            plot(fit_curve_1 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_2 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_3 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            fit_curve_tot = fit_curve_1+fit_curve_2+fit_curve_3;
            plot(fit_curve_tot ,'color',[0.0 1.0 1.0],'linestyle','-','linewidth',2);
            hold off;
        end;

        x_peak = [];
        x_peak(1) = G{1}.b;
        x_peak(2) = G{2}.b;
        x_peak(3) = G{3}.b;

		if (MakeFigure == true)
            hold on;
            for i = 1:3
                plot([x_peak(i) x_peak(i)],[0 max(raw)],'color',colorPeak,'linestyle','-','linewidth',2);
            end;
            for i = 1:length(actual_cutoffs)
                plot([actual_cutoffs(i) actual_cutoffs(i)],[0 max(raw)],'color',colorCutoff,'linestyle','-','linewidth',2);
            end;
            hold off;
        end;
    else   % procede as if approximates trisomy.
        if (MakeFigure == true)
            figure(fig);
            hold on;
            time1_1 = 1:floor(G{1}.b);
            time1_2 = ceil(G{1}.b):200;
            if (time1_1(end) == time1_2(1));    time1_1(end) = [];  end;
            time2_1 = 1:floor(G{2}.b);
            time2_2 = ceil(G{2}.b):200;
            if (time2_1(end) == time2_2(1));    time2_1(end) = [];  end;
            time3_1 = 1:floor(G{3}.b);
            time3_2 = ceil(G{3}.b):200;
            if (time3_1(end) == time3_2(1));    time3_2(1) = [];    end;
            time4_1 = 1:floor(G{4}.b);
            time4_2 = ceil(G{4}.b):200;
            if (time4_1(end) == time4_2(1));    time4_2(1) = [];    end;

			fit_curve_1 = [G{1}.a*exp(-0.5*((time1_1-G{1}.b)./G{1}.c).^2) ...
                           G{1}.a*exp(-0.5*((time1_2-G{1}.b)./G{1}.c/(skew_factor/G{1}.b)).^2)];
            fit_curve_2 = [G{2}.a*exp(-0.5*((time2_1-G{2}.b)./G{2}.c).^2) ...
                           G{2}.a*exp(-0.5*((time2_2-G{2}.b)./G{2}.c/(skew_factor/G{2}.b)).^2)];
            fit_curve_3 = [G{3}.a*exp(-0.5*((time3_1-G{3}.b)./G{3}.c/(skew_factor/(200-G{3}.b))).^2) ...
                           G{3}.a*exp(-0.5*((time3_2-G{3}.b)./G{3}.c).^2)];
            fit_curve_4 = [G{4}.a*exp(-0.5*((time4_1-G{4}.b)./G{4}.c/(skew_factor/(200-G{4}.b))).^2) ...
                           G{4}.a*exp(-0.5*((time4_2-G{4}.b)./G{4}.c).^2)];

            area(           max([1 round(G{1}.b-G{1}.c)]):min([200 round(G{1}.b+G{1}.c*(skew_factor/G{1}.b))]), ...
                fit_curve_1(max([1 round(G{1}.b-G{1}.c)]):min([200 round(G{1}.b+G{1}.c*(skew_factor/G{1}.b))])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{2}.b-G{2}.c)]):min([200 round(G{2}.b+G{2}.c*(skew_factor/G{2}.b))]), ...
                fit_curve_2(max([1 round(G{2}.b-G{2}.c)]):min([200 round(G{2}.b+G{2}.c*(skew_factor/G{2}.b))])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{3}.b-G{3}.c*(skew_factor/(200-G{3}.b)))]):min([200 round(G{3}.b+G{3}.c)]), ...
                fit_curve_3(max([1 round(G{3}.b-G{3}.c*(skew_factor/(200-G{3}.b)))]):min([200 round(G{3}.b+G{3}.c)])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{4}.b-G{4}.c*(skew_factor/(200-G{4}.b)))]):min([200 round(G{4}.b+G{4}.c)]), ...
                fit_curve_4(max([1 round(G{4}.b-G{4}.c*(skew_factor/(200-G{4}.b)))]):min([200 round(G{4}.b+G{4}.c)])),'FaceColor',[0.5000 1.0 0.5000]);
            plot(fit_curve_1 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_2 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_3 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_4 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            fit_curve_tot = fit_curve_1+fit_curve_2+fit_curve_3+fit_curve_4;
            plot(fit_curve_tot ,'color',[0.0 1.0 1.0],'linestyle','-','linewidth',2);
            hold off;
        end;

        x_peak = [];
        x_peak(1) = G{1}.b;
        x_peak(2) = G{2}.b;
		x_peak(3) = G{3}.b;
        x_peak(4) = G{4}.b;

        if (MakeFigure == true)
            hold on;
            for i = 1:4
                plot([x_peak(i) x_peak(i)],[0 max(raw)],'color',colorPeak,'linestyle','-','linewidth',2);
            end;
            for i = 1:length(actual_cutoffs)
                plot([actual_cutoffs(i) actual_cutoffs(i)],[0 max(raw)],'color',colorCutoff,'linestyle','-','linewidth',2);
            end;
            hold off;
        end;
    end;
elseif (chrCopyNum{chromosome}(segment) <= 4)
    fraction = 4-chrCopyNum{chromosome}(segment);
    peaks = [];
    peaks(1) = tetrasomy_peak(1);
    peaks(2) = tetrasomy_peak(2)*(1-fraction) + trisomy_peak(2)*fraction;
    peaks(3) = tetrasomy_peak(3);
    peaks(4) = tetrasomy_peak(4)*(1-fraction) + trisomy_peak(3)*fraction;
    peaks(5) = tetrasomy_peak(5);
    G = [];
    [G{1}.a,G{1}.b,G{1}.c, G{2}.a,G{2}.b,G{2}.c, G{3}.a,G{3}.b,G{3}.c, G{4}.a,G{4}.b,G{4}.c, G{5}.a,G{5}.b,G{5}.c, skew_factor] = ...
        fit_Gaussian_model_tetrasomy_2(smoothed,peaks*200,sigma,fraction,skew_factor,ErrorType,show_fitting);
    G{1}.d = skew_factor;
    G{2}.d = skew_factor;
    G{3}.d = skew_factor;
	G{4}.d = skew_factor;
    G{5}.d = skew_factor;
    [list] = FindHighestGaussian_2(G);

    actual_cutoffs = [];
    mostLikelyGaussians = [];
    for i = 1:199
        if (list(i) ~= list(i+1))   % we've found a boundary.
            actual_cutoffs = [actual_cutoffs FindGaussianCrossover_2(G{list(i)},G{list(i+1)},i)];
            mostLikelyGaussians = [mostLikelyGaussians list(i)];
        end;
    end;
    mostLikelyGaussians = [mostLikelyGaussians list(200)];

    % If central two gaussians overlap by their 1st SD, then assume one
    % central peak and recalculate.   The difference from 2n is not
    % significant.
    if (G{3}.b+G{3}.c*OverlapLimit > G{4}.b-G{4}.c*OverlapLimit)   % recalculate as if trisomy.
        chrCopyNum{chromosome}(segment) = 3;
        peaks = [];
        peaks(1) = trisomy_peak(1);
        peaks(2) = trisomy_peak(2);
        peaks(3) = trisomy_peak(3);
        peaks(4) = trisomy_peak(4);
        G = [];
		[G{1}.a,G{1}.b,G{1}.c, G{2}.a,G{2}.b,G{2}.c, G{3}.a,G{3}.b,G{3}.c, G{4}.a,G{4}.b,G{4}.c, skew_factor] = ...
            fit_Gaussian_model_trisomy_2(smoothed,peaks*200,sigma,fraction,skew_factor,ErrorType,show_fitting);
        G{1}.d = skew_factor;
        G{2}.d = skew_factor;
        G{3}.d = skew_factor;
        G{4}.d = skew_factor;
        [list] = FindHighestGaussian_2(G);

        actual_cutoffs = [];
        mostLikelyGaussians = [];
        for i = 1:199
            if (list(i) ~= list(i+1))   % we've found a boundary.
                actual_cutoffs = [actual_cutoffs FindGaussianCrossover_2(G{list(i)},G{list(i+1)},i)];
                mostLikelyGaussians = [mostLikelyGaussians list(i)];
            end;
        end;
        mostLikelyGaussians = [mostLikelyGaussians list(200)];

        if (MakeFigure == true)
            figure(fig);
            hold on;
            time1_1 = 1:floor(G{1}.b);
            time1_2 = ceil(G{1}.b):200;
            if (time1_1(end) == time1_2(1));    time1_1(end) = [];  end;
            time2_1 = 1:floor(G{2}.b);
            time2_2 = ceil(G{2}.b):200;
            if (time2_1(end) == time2_2(1));    time2_1(end) = [];  end;
            time3_1 = 1:floor(G{3}.b);
            time3_2 = ceil(G{3}.b):200;
            if (time3_1(end) == time3_2(1));    time3_2(1) = [];    end;
			time4_1 = 1:floor(G{4}.b);
            time4_2 = ceil(G{4}.b):200;
            if (time4_1(end) == time4_2(1));    time4_2(1) = [];    end;

            fit_curve_1 = [G{1}.a*exp(-0.5*((time1_1-G{1}.b)./G{1}.c).^2) ...
                           G{1}.a*exp(-0.5*((time1_2-G{1}.b)./G{1}.c/(skew_factor/G{1}.b)).^2)];
            fit_curve_2 = [G{2}.a*exp(-0.5*((time2_1-G{2}.b)./G{2}.c).^2) ...
                           G{2}.a*exp(-0.5*((time2_2-G{2}.b)./G{2}.c/(skew_factor/G{2}.b)).^2)];
            fit_curve_3 = [G{3}.a*exp(-0.5*((time3_1-G{3}.b)./G{3}.c/(skew_factor/(200-G{3}.b))).^2) ...
                           G{3}.a*exp(-0.5*((time3_2-G{3}.b)./G{3}.c).^2)];
            fit_curve_4 = [G{4}.a*exp(-0.5*((time4_1-G{4}.b)./G{4}.c/(skew_factor/(200-G{4}.b))).^2) ...
                           G{4}.a*exp(-0.5*((time4_2-G{4}.b)./G{4}.c).^2)];

            area(           max([1 round(G{1}.b-G{1}.c)]):min([200 round(G{1}.b+G{1}.c*(skew_factor/G{1}.b))]), ...
                fit_curve_1(max([1 round(G{1}.b-G{1}.c)]):min([200 round(G{1}.b+G{1}.c*(skew_factor/G{1}.b))])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{2}.b-G{2}.c)]):min([200 round(G{2}.b+G{2}.c*(skew_factor/G{2}.b))]), ...
                fit_curve_2(max([1 round(G{2}.b-G{2}.c)]):min([200 round(G{2}.b+G{2}.c*(skew_factor/G{2}.b))])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{3}.b-G{3}.c*(skew_factor/(200-G{3}.b)))]):min([200 round(G{3}.b+G{3}.c)]), ...
                fit_curve_3(max([1 round(G{3}.b-G{3}.c*(skew_factor/(200-G{3}.b)))]):min([200 round(G{3}.b+G{3}.c)])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{4}.b-G{4}.c*(skew_factor/(200-G{4}.b)))]):min([200 round(G{4}.b+G{4}.c)]), ...
                fit_curve_4(max([1 round(G{4}.b-G{4}.c*(skew_factor/(200-G{4}.b)))]):min([200 round(G{4}.b+G{4}.c)])),'FaceColor',[0.5000 1.0 0.5000]);
            plot(fit_curve_1 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_2 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_3 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_4 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            fit_curve_tot = fit_curve_1+fit_curve_2+fit_curve_3+fit_curve_4;
            plot(fit_curve_tot ,'color',[0.0 1.0 1.0],'linestyle','-','linewidth',2);
            hold off;
        end;

		x_peak = [];
        x_peak(1) = G{1}.b;
        x_peak(2) = G{2}.b;
        x_peak(3) = G{3}.b;
        x_peak(4) = G{4}.b;

        if (MakeFigure == true)
            hold on;
            for i = 1:4
                plot([x_peak(i) x_peak(i)],[0 max(raw)],'color',colorPeak,'linestyle','-','linewidth',2);
            end;
            for i = 1:length(actual_cutoffs)
                plot([actual_cutoffs(i) actual_cutoffs(i)],[0 max(raw)],'color',colorCutoff,'linestyle','-','linewidth',2);
            end;
            hold off;
        end;
    else   % procede as if approximates tetrasomy.
        if (MakeFigure == true)
            figure(fig);
            hold on;
            time1_1 = 1:floor(G{1}.b);
            time1_2 = ceil(G{1}.b):200;
            if (time1_1(end) == time1_2(1));    time1_1(end) = [];  end;
            time2_1 = 1:floor(G{2}.b);
            time2_2 = ceil(G{2}.b):200;
            if (time2_1(end) == time2_2(1));    time2_1(end) = [];  end;
			time3   = 1:200;
            time4_1 = 1:floor(G{4}.b);
            time4_2 = ceil(G{4}.b):200;
            if (time4_1(end) == time4_2(1));    time4_2(1) = [];    end;
            time5_1 = 1:floor(G{5}.b);
            time5_2 = ceil(G{5}.b):200;
            if (time5_1(end) == time5_2(1));    time5_2(1) = [];    end;

            fit_curve_1 = [G{1}.a*exp(-0.5*((time1_1-G{1}.b)./G{1}.c).^2) ...
                           G{1}.a*exp(-0.5*((time1_2-G{1}.b)./G{1}.c/(skew_factor/G{1}.b)).^2)];
            fit_curve_2 = [G{2}.a*exp(-0.5*((time2_1-G{2}.b)./G{2}.c).^2) ...
                           G{2}.a*exp(-0.5*((time2_2-G{2}.b)./G{2}.c/(skew_factor/G{2}.b)).^2)];
            fit_curve_3 =  G{3}.a*exp(-0.5*((time3-G{3}.b)./G{3}.c).^2);
            fit_curve_4 = [G{4}.a*exp(-0.5*((time4_1-G{4}.b)./G{4}.c/(skew_factor/(200-G{4}.b))).^2) ...
                           G{4}.a*exp(-0.5*((time4_2-G{4}.b)./G{4}.c).^2)];
            fit_curve_5 = [G{5}.a*exp(-0.5*((time5_1-G{5}.b)./G{5}.c/(skew_factor/(200-G{5}.b))).^2) ...
                           G{5}.a*exp(-0.5*((time5_2-G{5}.b)./G{5}.c).^2)];

            area(           max([1 round(G{1}.b-G{1}.c)]):min([200 round(G{1}.b+G{1}.c*(skew_factor/G{1}.b))]), ...
                fit_curve_1(max([1 round(G{1}.b-G{1}.c)]):min([200 round(G{1}.b+G{1}.c*(skew_factor/G{1}.b))])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{2}.b-G{2}.c)]):min([200 round(G{2}.b+G{2}.c*(skew_factor/G{2}.b))]), ...
                fit_curve_2(max([1 round(G{2}.b-G{2}.c)]):min([200 round(G{2}.b+G{2}.c*(skew_factor/G{2}.b))])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{3}.b-G{3}.c)]):min([200 round(G{3}.b+G{3}.c)]), ...
                fit_curve_3(max([1 round(G{3}.b-G{3}.c)]):min([200 round(G{3}.b+G{3}.c)])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{4}.b-G{4}.c*(skew_factor/(200-G{4}.b)))]):min([200 round(G{4}.b+G{4}.c)]), ...
                fit_curve_4(max([1 round(G{4}.b-G{4}.c*(skew_factor/(200-G{4}.b)))]):min([200 round(G{4}.b+G{4}.c)])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{5}.b-G{5}.c*(skew_factor/(200-G{5}.b)))]):min([200 round(G{5}.b+G{5}.c)]), ...
                fit_curve_5(max([1 round(G{5}.b-G{5}.c*(skew_factor/(200-G{5}.b)))]):min([200 round(G{5}.b+G{5}.c)])),'FaceColor',[0.5000 1.0 0.5000]);

			plot(fit_curve_1 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_2 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_3 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_4 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_5 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            fit_curve_tot = fit_curve_1+fit_curve_2+fit_curve_3+fit_curve_4+fit_curve_5;
            plot(fit_curve_tot ,'color',[0.0 1.0 1.0],'linestyle','-','linewidth',2);
            hold off;
        end;

        x_peak = [];
        x_peak(1) = G{1}.b;
        x_peak(2) = G{2}.b;
        x_peak(3) = G{3}.b;
        x_peak(4) = G{4}.b;
        x_peak(5) = G{5}.b;

        if (MakeFigure == true)
            hold on;
            for i = 1:5
                plot([x_peak(i) x_peak(i)],[0 max(raw)],'color',colorPeak,'linestyle','-','linewidth',2);
            end;
            for i = 1:length(actual_cutoffs)
                plot([actual_cutoffs(i) actual_cutoffs(i)],[0 max(raw)],'color',colorCutoff,'linestyle','-','linewidth',2);
            end;
            hold off;
        end;
    end;
elseif (chrCopyNum{chromosome}(segment) <= 5)
    fraction = 5-chrCopyNum{chromosome}(segment);
    peaks = [];
    peaks(1) = pentasomy_peak(1);
    peaks(2) = pentasomy_peak(2)*(1-fraction) + tetrasomy_peak(2)*fraction;
    peaks(3) = pentasomy_peak(3)*(1-fraction) + tetrasomy_peak(3)*fraction;
    peaks(4) = pentasomy_peak(4)*(1-fraction) + tetrasomy_peak(3)*fraction;
    peaks(5) = pentasomy_peak(5)*(1-fraction) + tetrasomy_peak(4)*fraction;
    peaks(6) = pentasomy_peak(6);
    G = [];
    [G{1}.a,G{1}.b,G{1}.c, G{2}.a,G{2}.b,G{2}.c, G{3}.a,G{3}.b,G{3}.c, G{4}.a,G{4}.b,G{4}.c, G{5}.a,G{5}.b,G{5}.c, G{6}.a,G{6}.b,G{6}.c, skew_factor] = ...
        fit_Gaussian_model_pentasomy_2(smoothed,peaks*200,sigma,fraction,skew_factor,ErrorType,show_fitting);
    G{1}.d = skew_factor;
    G{2}.d = skew_factor;
    G{3}.d = skew_factor;
    G{4}.d = skew_factor;
    G{5}.d = skew_factor;
    G{6}.d = skew_factor;
    [list] = FindHighestGaussian_2(G);

    actual_cutoffs = [];
    mostLikelyGaussians = [];
    for i = 1:199
        if (list(i) ~= list(i+1))   % we've found a boundary.
            actual_cutoffs = [actual_cutoffs FindGaussianCrossover_2(G{list(i)},G{list(i+1)},i)];
            mostLikelyGaussians = [mostLikelyGaussians list(i)];
        end;
    end;
	mostLikelyGaussians = [mostLikelyGaussians list(200)];

    % If central two gaussians overlap by their 1st SD, then assume one
    % central peak and recalculate.   The difference from 2n is not
    % significant.
    if (G{3}.b+G{3}.c*OverlapLimit > G{4}.b-G{4}.c*OverlapLimit)   % recalculate as if tetrasomy.
        chrCopyNum{chromosome}(segment) = 4;
        fraction = 0;
        peaks = [];
        peaks(1) = tetrasomy_peak(1);
        peaks(2) = tetrasomy_peak(2)*(1-fraction) + trisomy_peak(2)*fraction;
        peaks(3) = tetrasomy_peak(3);
        peaks(4) = tetrasomy_peak(4)*(1-fraction) + trisomy_peak(3)*fraction;
        peaks(5) = tetrasomy_peak(5);
        G = [];
        [G{1}.a,G{1}.b,G{1}.c, G{2}.a,G{2}.b,G{2}.c, G{3}.a,G{3}.b,G{3}.c, G{4}.a,G{4}.b,G{4}.c, G{5}.a,G{5}.b,G{5}.c, skew_factor] = ...
            fit_Gaussian_model_tetrasomy_2(smoothed,peaks*200,sigma,fraction,skew_factor,ErrorType,show_fitting);
        G{1}.d = skew_factor;
        G{2}.d = skew_factor;
        G{3}.d = skew_factor;
        G{4}.d = skew_factor;
        G{5}.d = skew_factor;
        [list] = FindHighestGaussian_2(G);

        actual_cutoffs = [];
        mostLikelyGaussians = [];
        for i = 1:199
            if (list(i) ~= list(i+1))   % we've found a boundary.
				actual_cutoffs = [actual_cutoffs FindGaussianCrossover_2(G{list(i)},G{list(i+1)},i)];
                mostLikelyGaussians = [mostLikelyGaussians list(i)];
            end;
        end;
        mostLikelyGaussians = [mostLikelyGaussians list(200)];

        if (MakeFigure == true)
            figure(fig);
            hold on;
            time1_1 = 1:floor(G{1}.b);
            time1_2 = ceil(G{1}.b):200;
            if (time1_1(end) == time1_2(1));    time1_1(end) = [];  end;
            time2_1 = 1:floor(G{2}.b);
            time2_2 = ceil(G{2}.b):200;
            if (time2_1(end) == time2_2(1));    time2_1(end) = [];  end;
            time3   = 1:200;
            time4_1 = 1:floor(G{4}.b);
            time4_2 = ceil(G{4}.b):200;
            if (time4_1(end) == time4_2(1));    time4_2(1) = [];    end;
            time5_1 = 1:floor(G{5}.b);
            time5_2 = ceil(G{5}.b):200;
            if (time5_1(end) == time5_2(1));    time5_2(1) = [];    end;

            fit_curve_1 = [G{1}.a*exp(-0.5*((time1_1-G{1}.b)./G{1}.c).^2) ...
                           G{1}.a*exp(-0.5*((time1_2-G{1}.b)./G{1}.c/(skew_factor/G{1}.b)).^2)];
            fit_curve_2 = [G{2}.a*exp(-0.5*((time2_1-G{2}.b)./G{2}.c).^2) ...
                           G{2}.a*exp(-0.5*((time2_2-G{2}.b)./G{2}.c/(skew_factor/G{2}.b)).^2)];
            fit_curve_3 =  G{3}.a*exp(-0.5*((time3-G{3}.b)./G{3}.c).^2);
			fit_curve_4 = [G{4}.a*exp(-0.5*((time4_1-G{4}.b)./G{4}.c/(skew_factor/(200-G{4}.b))).^2) ...
                           G{4}.a*exp(-0.5*((time4_2-G{4}.b)./G{4}.c).^2)];
            fit_curve_5 = [G{5}.a*exp(-0.5*((time5_1-G{5}.b)./G{5}.c/(skew_factor/(200-G{5}.b))).^2) ...
                           G{5}.a*exp(-0.5*((time5_2-G{5}.b)./G{5}.c).^2)];

            area(           max([1 round(G{1}.b-G{1}.c)]):min([200 round(G{1}.b+G{1}.c*(skew_factor/G{1}.b))]), ...
                fit_curve_1(max([1 round(G{1}.b-G{1}.c)]):min([200 round(G{1}.b+G{1}.c*(skew_factor/G{1}.b))])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{2}.b-G{2}.c)]):min([200 round(G{2}.b+G{2}.c*(skew_factor/G{2}.b))]), ...
                fit_curve_2(max([1 round(G{2}.b-G{2}.c)]):min([200 round(G{2}.b+G{2}.c*(skew_factor/G{2}.b))])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{3}.b-G{3}.c)]):min([200 round(G{3}.b+G{3}.c)]), ...
                fit_curve_3(max([1 round(G{3}.b-G{3}.c)]):min([200 round(G{3}.b+G{3}.c)])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{4}.b-G{4}.c*(skew_factor/(200-G{4}.b)))]):min([200 round(G{4}.b+G{4}.c)]), ...
                fit_curve_4(max([1 round(G{4}.b-G{4}.c*(skew_factor/(200-G{4}.b)))]):min([200 round(G{4}.b+G{4}.c)])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{5}.b-G{5}.c*(skew_factor/(200-G{5}.b)))]):min([200 round(G{5}.b+G{5}.c)]), ...
                fit_curve_5(max([1 round(G{5}.b-G{5}.c*(skew_factor/(200-G{5}.b)))]):min([200 round(G{5}.b+G{5}.c)])),'FaceColor',[0.5000 1.0 0.5000]);

            plot(fit_curve_1 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_2 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_3 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_4 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_5 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            fit_curve_tot = fit_curve_1+fit_curve_2+fit_curve_3+fit_curve_4+fit_curve_5;
            plot(fit_curve_tot ,'color',[0.0 1.0 1.0],'linestyle','-','linewidth',2);
            hold off;
        end;

        x_peak = [];
        x_peak(1) = G{1}.b;
		x_peak(2) = G{2}.b;
        x_peak(3) = G{3}.b;
        x_peak(4) = G{4}.b;
        x_peak(5) = G{5}.b;

        if (MakeFigure == true)
            hold on;
            for i = 1:5
                plot([x_peak(i) x_peak(i)],[0 max(raw)],'color',colorPeak,'linestyle','-','linewidth',2);
            end;
            for i = 1:length(actual_cutoffs)
                plot([actual_cutoffs(i) actual_cutoffs(i)],[0 max(raw)],'color',colorCutoff,'linestyle','-','linewidth',2);
            end;
            hold off;
        end;
    else % procede as if approximates pentasomy.
        if (MakeFigure == true)
            figure(fig);
            hold on;
            time1_1 = 1:floor(G{1}.b);
            time1_2 = ceil(G{1}.b):200;
            if (time1_1(end) == time1_2(1));    time1_1(end) = [];  end;
            time2_1 = 1:floor(G{2}.b);
            time2_2 = ceil(G{2}.b):200;
            if (time2_1(end) == time2_2(1));    time2_1(end) = [];  end;
            time3_1 = 1:floor(G{3}.b);
            time3_2 = ceil(G{3}.b):200;
            if (time3_1(end) == time3_2(1));    time3_1(end) = [];  end;
			time4_1 = 1:floor(G{4}.b);
            time4_2 = ceil(G{4}.b):200;
            if (time4_1(end) == time4_2(1));    time4_2(1) = [];    end;
            time5_1 = 1:floor(G{5}.b);
            time5_2 = ceil(G{5}.b):200;
            if (time5_1(end) == time5_2(1));    time5_2(1) = [];    end;
            time6_1 = 1:floor(G{6}.b);
            time6_2 = ceil(G{6}.b):200;
            if (time6_1(end) == time6_2(1));    time6_2(1) = [];    end;

            fit_curve_1 = [G{1}.a*exp(-0.5*((time1_1-G{1}.b)./G{1}.c).^2) ...
                           G{1}.a*exp(-0.5*((time1_2-G{1}.b)./G{1}.c/(skew_factor/G{1}.b)).^2)];
            fit_curve_2 = [G{2}.a*exp(-0.5*((time2_1-G{2}.b)./G{2}.c).^2) ...
                           G{2}.a*exp(-0.5*((time2_2-G{2}.b)./G{2}.c/(skew_factor/G{2}.b)).^2)];
            fit_curve_3 = [G{3}.a*exp(-0.5*((time3_1-G{3}.b)./G{3}.c).^2) ...
                           G{3}.a*exp(-0.5*((time3_2-G{3}.b)./G{3}.c/(skew_factor/G{3}.b)).^2)];
            fit_curve_4 = [G{4}.a*exp(-0.5*((time4_1-G{4}.b)./G{4}.c/(skew_factor/(200-G{4}.b))).^2) ...
                           G{4}.a*exp(-0.5*((time4_2-G{4}.b)./G{4}.c).^2)];
            fit_curve_5 = [G{5}.a*exp(-0.5*((time5_1-G{5}.b)./G{5}.c/(skew_factor/(200-G{5}.b))).^2) ...
                           G{5}.a*exp(-0.5*((time5_2-G{5}.b)./G{5}.c).^2)];
            fit_curve_6 = [G{6}.a*exp(-0.5*((time6_1-G{6}.b)./G{6}.c/(skew_factor/(200-G{6}.b))).^2) ...
                           G{6}.a*exp(-0.5*((time6_2-G{6}.b)./G{6}.c).^2)];

            area(           max([1 round(G{1}.b-G{1}.c)]):min([200 round(G{1}.b+G{1}.c*(skew_factor/G{1}.b))]), ...
                fit_curve_1(max([1 round(G{1}.b-G{1}.c)]):min([200 round(G{1}.b+G{1}.c*(skew_factor/G{1}.b))])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{2}.b-G{2}.c)]):min([200 round(G{2}.b+G{2}.c*(skew_factor/G{2}.b))]), ...
                fit_curve_2(max([1 round(G{2}.b-G{2}.c)]):min([200 round(G{2}.b+G{2}.c*(skew_factor/G{2}.b))])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{3}.b-G{3}.c)]):min([200 round(G{3}.b+G{3}.c*(skew_factor/G{3}.b))]), ...
				fit_curve_3(max([1 round(G{3}.b-G{3}.c)]):min([200 round(G{3}.b+G{3}.c*(skew_factor/G{3}.b))])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{4}.b-G{4}.c*(skew_factor/(200-G{4}.b)))]):min([200 round(G{4}.b+G{4}.c)]), ...
                fit_curve_4(max([1 round(G{4}.b-G{4}.c*(skew_factor/(200-G{4}.b)))]):min([200 round(G{4}.b+G{4}.c)])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{5}.b-G{5}.c*(skew_factor/(200-G{5}.b)))]):min([200 round(G{5}.b+G{5}.c)]), ...
                fit_curve_5(max([1 round(G{5}.b-G{5}.c*(skew_factor/(200-G{5}.b)))]):min([200 round(G{5}.b+G{5}.c)])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{6}.b-G{6}.c*(skew_factor/(200-G{6}.b)))]):min([200 round(G{6}.b+G{6}.c)]), ...
                fit_curve_6(max([1 round(G{6}.b-G{6}.c*(skew_factor/(200-G{6}.b)))]):min([200 round(G{6}.b+G{6}.c)])),'FaceColor',[0.5000 1.0 0.5000]);

            plot(fit_curve_1 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_2 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_3 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_4 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_5 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_6 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            fit_curve_tot = fit_curve_1+fit_curve_2+fit_curve_3+fit_curve_4+fit_curve_5+fit_curve_6;
            plot(fit_curve_tot ,'color',[0.0 1.0 1.0],'linestyle','-','linewidth',2);
            hold off;
        end;

        x_peak = [];
        x_peak(1) = G{1}.b;
        x_peak(2) = G{2}.b;
        x_peak(3) = G{3}.b;
        x_peak(4) = G{4}.b;
        x_peak(5) = G{5}.b;
        x_peak(6) = G{6}.b;

        if (MakeFigure == true)
			hold on;
            for i = 1:6
                plot([x_peak(i) x_peak(i)],[0 max(raw)],'color',colorPeak,'linestyle','-','linewidth',2);
            end;
            for i = 1:length(actual_cutoffs)
                plot([actual_cutoffs(i) actual_cutoffs(i)],[0 max(raw)],'color',colorCutoff,'linestyle','-','linewidth',2);
            end;
            hold off;
        end;
    end;
else % (chrCopyNum{chromosome}(segment) <= 6)
    fraction = 6-chrCopyNum{chromosome}(segment);
    peaks = [];
    peaks(1) = hexasomy_peak(1);
    peaks(2) = hexasomy_peak(2)*(1-fraction) + pentasomy_peak(2)*fraction;
    peaks(3) = hexasomy_peak(3)*(1-fraction) + pentasomy_peak(3)*fraction;
    peaks(4) = hexasomy_peak(4);
    peaks(5) = hexasomy_peak(5)*(1-fraction) + pentasomy_peak(4)*fraction;
    peaks(6) = hexasomy_peak(6)*(1-fraction) + pentasomy_peak(5)*fraction;
    peaks(7) = hexasomy_peak(7);
    G = [];
    [G{1}.a,G{1}.b,G{1}.c, G{2}.a,G{2}.b,G{2}.c, G{3}.a,G{3}.b,G{3}.c, G{4}.a,G{4}.b,G{4}.c, G{5}.a,G{5}.b,G{5}.c, G{6}.a,G{6}.b,G{6}.c, G{7}.a,G{7}.b,G{7}.c, skew_factor] = ...
        fit_Gaussian_model_hexasomy_2(smoothed,peaks*200,sigma,fraction,skew_factor,ErrorType,show_fitting);
    G{1}.d = skew_factor;
    G{2}.d = skew_factor;
    G{3}.d = skew_factor;
    G{4}.d = skew_factor;
    G{5}.d = skew_factor;
	G{6}.d = skew_factor;
    G{7}.d = skew_factor;
    [list] = FindHighestGaussian_2(G);

    actual_cutoffs = [];
    mostLikelyGaussians = [];
    for i = 1:199
        if (list(i) ~= list(i+1))   % we've found a boundary.
            actual_cutoffs = [actual_cutoffs FindGaussianCrossover_2(G{list(i)},G{list(i+1)},i)];
            mostLikelyGaussians = [mostLikelyGaussians list(i)];
        end;
    end;
    mostLikelyGaussians = [mostLikelyGaussians list(200)];

    % If central two gaussians overlap by their 1st SD, then assume one
    % central peak and recalculate.   The difference from 2n is not
    % significant.
    if (G{4}.b+G{4}.c*OverlapLimit > G{5}.b-G{5}.c*OverlapLimit)   % recalculate as if pentasomy.
        chrCopyNum{chromosome}(segment) = 5;
        fraction = 0;
        peaks = [];
        peaks(1) = pentasomy_peak(1);
        peaks(2) = pentasomy_peak(2);
        peaks(3) = pentasomy_peak(3);
        peaks(4) = pentasomy_peak(4);
        peaks(5) = pentasomy_peak(5);
        peaks(6) = pentasomy_peak(6);
        G = [];
		[G{1}.a,G{1}.b,G{1}.c, G{2}.a,G{2}.b,G{2}.c, G{3}.a,G{3}.b,G{3}.c, G{4}.a,G{4}.b,G{4}.c, G{5}.a,G{5}.b,G{5}.c, G{6}.a,G{6}.b,G{6}.c, skew_factor] = ...
            fit_Gaussian_model_pentasomy_5(smoothed,peaks*200,sigma,fraction,skew_factor,ErrorType,show_fitting);
        G{1}.d = skew_factor;
        G{2}.d = skew_factor;
        G{3}.d = skew_factor;
        G{4}.d = skew_factor;
        G{5}.d = skew_factor;
        G{6}.d = skew_factor;
        [list] = FindHighestGaussian_2(G);

        actual_cutoffs = [];
        mostLikelyGaussians = [];
        for i = 1:199
            if (list(i) ~= list(i+1))   % we've found a boundary.
                actual_cutoffs = [actual_cutoffs FindGaussianCrossover_2(G{list(i)},G{list(i+1)},i)];
                mostLikelyGaussians = [mostLikelyGaussians list(i)];
            end;
        end;
        mostLikelyGaussians = [mostLikelyGaussians list(200)];

        if (MakeFigure == true)
            figure(fig);
            hold on;
            time1_1 = 1:floor(G{1}.b);
            time1_2 = ceil(G{1}.b):200;
            if (time1_1(end) == time1_2(1));    time1_1(end) = [];  end;
            time2_1 = 1:floor(G{2}.b);
            time2_2 = ceil(G{2}.b):200;
			if (time2_1(end) == time2_2(1));    time2_1(end) = [];  end;
            time3_1 = 1:floor(G{3}.b);
            time3_2 = ceil(G{3}.b):200;
            if (time3_1(end) == time3_2(1));    time3_1(end) = [];  end;
            time4_1 = 1:floor(G{4}.b);
            time4_2 = ceil(G{4}.b):200;
            if (time4_1(end) == time4_2(1));    time4_2(1) = [];    end;
            time5_1 = 1:floor(G{5}.b);
            time5_2 = ceil(G{5}.b):200;
            if (time5_1(end) == time5_2(1));    time5_2(1) = [];    end;
            time6_1 = 1:floor(G{6}.b);
            time6_2 = ceil(G{6}.b):200;
            if (time6_1(end) == time6_2(1));    time6_2(1) = [];    end;

            fit_curve_1 = [G{1}.a*exp(-0.5*((time1_1-G{1}.b)./G{1}.c).^2) ...
                           G{1}.a*exp(-0.5*((time1_2-G{1}.b)./G{1}.c/(skew_factor/G{1}.b)).^2)];
            fit_curve_2 = [G{2}.a*exp(-0.5*((time2_1-G{2}.b)./G{2}.c).^2) ...
                           G{2}.a*exp(-0.5*((time2_2-G{2}.b)./G{2}.c/(skew_factor/G{2}.b)).^2)];
            fit_curve_3 = [G{3}.a*exp(-0.5*((time3_1-G{3}.b)./G{3}.c).^2) ...
                           G{3}.a*exp(-0.5*((time3_2-G{3}.b)./G{3}.c/(skew_factor/G{3}.b)).^2)];
            fit_curve_4 = [G{4}.a*exp(-0.5*((time4_1-G{4}.b)./G{4}.c/(skew_factor/(200-G{4}.b))).^2) ...
                           G{4}.a*exp(-0.5*((time4_2-G{4}.b)./G{4}.c).^2)];
            fit_curve_5 = [G{5}.a*exp(-0.5*((time5_1-G{5}.b)./G{5}.c/(skew_factor/(200-G{5}.b))).^2) ...
                           G{5}.a*exp(-0.5*((time5_2-G{5}.b)./G{5}.c).^2)];
            fit_curve_6 = [G{6}.a*exp(-0.5*((time6_1-G{6}.b)./G{6}.c/(skew_factor/(200-G{6}.b))).^2) ...
                           G{6}.a*exp(-0.5*((time6_2-G{6}.b)./G{6}.c).^2)];

            area(           max([1 round(G{1}.b-G{1}.c)]):min([200 round(G{1}.b+G{1}.c*(skew_factor/G{1}.b))]), ...
				fit_curve_1(max([1 round(G{1}.b-G{1}.c)]):min([200 round(G{1}.b+G{1}.c*(skew_factor/G{1}.b))])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{2}.b-G{2}.c)]):min([200 round(G{2}.b+G{2}.c*(skew_factor/G{2}.b))]), ...
                fit_curve_2(max([1 round(G{2}.b-G{2}.c)]):min([200 round(G{2}.b+G{2}.c*(skew_factor/G{2}.b))])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{3}.b-G{3}.c)]):min([200 round(G{3}.b+G{3}.c*(skew_factor/G{3}.b))]), ...
                fit_curve_3(max([1 round(G{3}.b-G{3}.c)]):min([200 round(G{3}.b+G{3}.c*(skew_factor/G{3}.b))])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{4}.b-G{4}.c*(skew_factor/(200-G{4}.b)))]):min([200 round(G{4}.b+G{4}.c)]), ...
                fit_curve_4(max([1 round(G{4}.b-G{4}.c*(skew_factor/(200-G{4}.b)))]):min([200 round(G{4}.b+G{4}.c)])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{5}.b-G{5}.c*(skew_factor/(200-G{5}.b)))]):min([200 round(G{5}.b+G{5}.c)]), ...
                fit_curve_5(max([1 round(G{5}.b-G{5}.c*(skew_factor/(200-G{5}.b)))]):min([200 round(G{5}.b+G{5}.c)])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{6}.b-G{6}.c*(skew_factor/(200-G{6}.b)))]):min([200 round(G{6}.b+G{6}.c)]), ...
                fit_curve_6(max([1 round(G{6}.b-G{6}.c*(skew_factor/(200-G{6}.b)))]):min([200 round(G{6}.b+G{6}.c)])),'FaceColor',[0.5000 1.0 0.5000]);

            plot(fit_curve_1 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_2 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_3 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_4 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_5 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_6 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            fit_curve_tot = fit_curve_1+fit_curve_2+fit_curve_3+fit_curve_4+fit_curve_5+fit_curve_6;
            plot(fit_curve_tot ,'color',[0.0 1.0 1.0],'linestyle','-','linewidth',2);
            hold off;
        end;

        x_peak = [];
        x_peak(1) = G{1}.b;
        x_peak(2) = G{2}.b;
        x_peak(3) = G{3}.b;
        x_peak(4) = G{4}.b;
		x_peak(5) = G{5}.b;
        x_peak(6) = G{6}.b;

        if (MakeFigure == true)
            hold on;
            for i = 1:6
                plot([x_peak(i) x_peak(i)],[0 max(raw)],'color',colorPeak,'linestyle','-','linewidth',2);
            end;
            for i = 1:length(actual_cutoffs)
                plot([actual_cutoffs(i) actual_cutoffs(i)],[0 max(raw)],'color',colorCutoff,'linestyle','-','linewidth',2);
            end;
            hold off;
        end;
    else   % procede as if approximate hexasomy.
        if (MakeFigure == true)
            figure(fig);
            hold on;
            time1_1 = 1:floor(G{1}.b);
            time1_2 = ceil(G{1}.b):200;
            if (time1_1(end) == time1_2(1));    time1_1(end) = [];  end;
            time2_1 = 1:floor(G{2}.b);
            time2_2 = ceil(G{2}.b):200;
            if (time2_1(end) == time2_2(1));    time2_1(end) = [];  end;
            time3_1 = 1:floor(G{3}.b);
            time3_2 = ceil(G{3}.b):200;
            if (time3_1(end) == time3_2(1));    time3_1(end) = [];  end;
            time4   = 1:200;
            time5_1 = 1:floor(G{5}.b);
			time5_2 = ceil(G{5}.b):200;
            if (time5_1(end) == time5_2(1));    time5_2(1) = [];    end;
            time6_1 = 1:floor(G{6}.b);
            time6_2 = ceil(G{6}.b):200;
            if (time6_1(end) == time6_2(1));    time6_2(1) = [];    end;
            time7_1 = 1:floor(G{7}.b);
            time7_2 = ceil(G{7}.b):200;
            if (time7_1(end) == time7_2(1));    time7_2(1) = [];    end;

            fit_curve_1 = [G{1}.a*exp(-0.5*((time1_1-G{1}.b)./G{1}.c).^2) ...
                           G{1}.a*exp(-0.5*((time1_2-G{1}.b)./G{1}.c/(skew_factor/G{1}.b)).^2)];
            fit_curve_2 = [G{2}.a*exp(-0.5*((time2_1-G{2}.b)./G{2}.c).^2) ...
                           G{2}.a*exp(-0.5*((time2_2-G{2}.b)./G{2}.c/(skew_factor/G{2}.b)).^2)];
            fit_curve_3 = [G{3}.a*exp(-0.5*((time3_1-G{3}.b)./G{3}.c).^2) ...
                           G{3}.a*exp(-0.5*((time3_2-G{3}.b)./G{3}.c/(skew_factor/G{3}.b)).^2)];
            fit_curve_4 =  G{4}.a*exp(-0.5*((time4-G{4}.b)./G{4}.c).^2);
            fit_curve_5 = [G{5}.a*exp(-0.5*((time5_1-G{5}.b)./G{5}.c/(skew_factor/(200-G{5}.b))).^2) ...
                           G{5}.a*exp(-0.5*((time5_2-G{5}.b)./G{5}.c).^2)];
            fit_curve_6 = [G{6}.a*exp(-0.5*((time6_1-G{6}.b)./G{6}.c/(skew_factor/(200-G{6}.b))).^2) ...
                           G{6}.a*exp(-0.5*((time6_2-G{6}.b)./G{6}.c).^2)];
            fit_curve_7 = [G{7}.a*exp(-0.5*((time7_1-G{7}.b)./G{7}.c/(skew_factor/(200-G{7}.b))).^2) ...
                           G{7}.a*exp(-0.5*((time7_2-G{7}.b)./G{7}.c).^2)];

            area(           max([1 round(G{1}.b-G{1}.c)]):min([200 round(G{1}.b+G{1}.c*(skew_factor/G{1}.b))]), ...
                fit_curve_1(max([1 round(G{1}.b-G{1}.c)]):min([200 round(G{1}.b+G{1}.c*(skew_factor/G{1}.b))])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{2}.b-G{2}.c)]):min([200 round(G{2}.b+G{2}.c*(skew_factor/G{2}.b))]), ...
                fit_curve_2(max([1 round(G{2}.b-G{2}.c)]):min([200 round(G{2}.b+G{2}.c*(skew_factor/G{2}.b))])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{3}.b-G{3}.c)]):min([200 round(G{3}.b+G{3}.c*(skew_factor/G{3}.b))]), ...
				fit_curve_3(max([1 round(G{3}.b-G{3}.c)]):min([200 round(G{3}.b+G{3}.c*(skew_factor/G{3}.b))])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{4}.b-G{4}.c)]):min([200 round(G{4}.b+G{4}.c)]), ...
                fit_curve_4(max([1 round(G{4}.b-G{4}.c)]):min([200 round(G{4}.b+G{4}.c)])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{5}.b-G{5}.c*(skew_factor/(200-G{5}.b)))]):min([200 round(G{5}.b+G{5}.c)]), ...
                fit_curve_5(max([1 round(G{5}.b-G{5}.c*(skew_factor/(200-G{5}.b)))]):min([200 round(G{5}.b+G{5}.c)])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{6}.b-G{6}.c*(skew_factor/(200-G{6}.b)))]):min([200 round(G{6}.b+G{6}.c)]), ...
                fit_curve_6(max([1 round(G{6}.b-G{6}.c*(skew_factor/(200-G{6}.b)))]):min([200 round(G{6}.b+G{6}.c)])),'FaceColor',[0.5000 1.0 0.5000]);
            area(           max([1 round(G{7}.b-G{7}.c*(skew_factor/(200-G{7}.b)))]):min([200 round(G{7}.b+G{7}.c)]), ...
                fit_curve_7(max([1 round(G{7}.b-G{7}.c*(skew_factor/(200-G{7}.b)))]):min([200 round(G{7}.b+G{7}.c)])),'FaceColor',[0.5000 1.0 0.5000]);

            plot(fit_curve_1 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_2 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_3 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_4 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_5 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_6 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            plot(fit_curve_7 ,'color',[1.0 1.0 0.0],'linestyle','-','linewidth',2);
            fit_curve_tot = fit_curve_1+fit_curve_2+fit_curve_3+fit_curve_4+fit_curve_5+fit_curve_6+fit_curve_7;
            plot(fit_curve_tot ,'color',[0.0 1.0 1.0],'linestyle','-','linewidth',2);
            hold off;
        end;

        x_peak = [];
        x_peak(1) = G{1}.b;
        x_peak(2) = G{2}.b;
        x_peak(3) = G{3}.b;
        x_peak(4) = G{4}.b;
        x_peak(5) = G{5}.b;
		x_peak(6) = G{6}.b;
        x_peak(7) = G{7}.b;

        if (MakeFigure == true)
            hold on;
            for i = 1:7
                plot([x_peak(i) x_peak(i)],[0 max(raw)],'color',colorPeak,'linestyle','-','linewidth',2);
            end;
            for i = 1:length(actual_cutoffs)
                plot([actual_cutoffs(i) actual_cutoffs(i)],[0 max(raw)],'color',colorCutoff,'linestyle','-','linewidth',2);
            end;
            hold off;
        end;
    end;
end;
if (MakeFigure == true)
    hold on;
    plot(raw     ,'color',[0.75 0.75 1.0],'linestyle','-','linewidth',1);
    plot(smoothed,'color',[0.50 0.50 1.0],'linestyle','-','linewidth',1);
    title([name '; chr ' num2str(chromosome) '; segment ' num2str(segment)],'HorizontalAlign','center','VerticalAlign','middle');
    hold off;
    xlim([1,200]);
    ylim([0,max(raw)]);
    % save then delete figures.
    if ispc  % Windows
        fig_dir = 'figures\scatterHist_perChr\';
        if (isdir([file_dir 'figures\scatterHist_perChr']) == 0)
            mkdir([file_dir 'figures\scatterHist_perChr']);
		end;
    else     % MacOS
        fig_dir = 'figures/scatterHist_perChr/';
        if (isdir([file_dir 'figures/scatterHist_perChr']) == 0)
            mkdir([file_dir 'figures/scatterHist_perChr']);
        end;
    end;
    %saveas(fig, [file_dir fig_dir strrep(name,' ','_') '_chr-' num2str(chromosome) '_seg-' num2str(segment) '.eps'], 'epsc');

    saveas(fig, [file_dir fig_dir strrep(name,' ','_') '_chr-' num2str(chromosome) '_seg-' num2str(segment) '.png'], 'png');
    delete(fig);
end;
end

