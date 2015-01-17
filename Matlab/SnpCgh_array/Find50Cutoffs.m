function [raw,smoothed,x_peaks,actual_cutoffs,mostLikelyGaussians,chrCopyNum] = Find50Cutoffs(probeset1,chrCopyNum,chr_breaks,chr_size,dataset,chromosome,segment, ...
    monosomy_peak,disomy_peak,trisomy_peak,tetrasomy_peak,pentasomy_peak,hexasomy_peak,name,file_dir,MakeFigure,~,DataTypeToUse, workingDir)

% Overlap of 1st SD before assume identity.

% initialize vectors for scattergram analysis
data1_0 = [];   data2_0 = [];   data1_1 = [];
data2_1 = [];   data1_2 = [];   data2_2 = [];

% Gathers SNP probe data for this chromosome segment.
for i = 1:2:length(probeset1)
    if (probeset1(i).probe_location > chr_breaks{chromosome}(segment)*chr_size(chromosome)) && (probeset1(i).probe_location <= chr_breaks{chromosome}(segment+1)*chr_size(chromosome))
        % determines if probe pair is useful; both probes have data.
        if (isnan(probeset1(i).probe_Ratio(dataset)) == 0) && (isnan(probeset1(i+1).probe_Ratio(dataset)) == 0)
            % swaps probes if polarity switch indicates to do so.
            if (probeset1(i).probe_chromosome == chromosome)
                val1 = probeset1(i).probe_Ratio(dataset);
                val2 = probeset1(i+1).probe_Ratio(dataset);
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

%% FindGaussianCutoffs Finds cutoffs as intersections of Gaussians, fit to
%    the data at each peak location.
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
    for i = 1:length(peaks)-1
        actual_cutoffs(i)      = (peaks(i)+peaks(i+1))/2*200;
        mostLikelyGaussians(i) = i;
    end;
    mostLikelyGaussians = [mostLikelyGaussians length(peaks)];
    
    for i = 1:length(peaks)
        x_peaks(i) = peaks(i)*200;
    end;
    if (MakeFigure == true)
        for i = 1:2
            plot([x_peaks(i) x_peaks(i)],[0 max(raw)],'color',colorPeak,'linestyle','-','linewidth',2);
        end;
    end;
    for i = 1:length(actual_cutoffs)
        if (MakeFigure == true)
            plot([actual_cutoffs(i) actual_cutoffs(i)],[0 max(raw)],'color',colorCutoff,'linestyle','-','linewidth',2);
        end;
    end;
elseif (chrCopyNum{chromosome}(segment) <= 2)
    fraction = 2-chrCopyNum{chromosome}(segment);
    peaks = [];
    peaks(1) = disomy_peak(1);
    peaks(2) = disomy_peak(2);
    peaks(3) = disomy_peak(3);
    for i = 1:length(peaks)-1
        actual_cutoffs(i)      = (peaks(i)+peaks(i+1))/2*200;
        mostLikelyGaussians(i) = i;
    end;
    mostLikelyGaussians = [mostLikelyGaussians length(peaks)];

    for i = 1:length(peaks)
        x_peaks(i) = peaks(i)*200;
    end;
    if (MakeFigure == true)
        for i = 1:3
            plot([x_peaks(i) x_peaks(i)],[0 max(raw)],'color',colorPeak,'linestyle','-','linewidth',2);
        end;
    end;    
    for i = 1:length(actual_cutoffs)
        if (MakeFigure == true)
            plot([actual_cutoffs(i) actual_cutoffs(i)],[0 max(raw)],'color',colorCutoff,'linestyle','-','linewidth',2);
        end;
    end;
elseif (chrCopyNum{chromosome}(segment) <= 3)
    fraction = 3-chrCopyNum{chromosome}(segment);
    peaks = [];
    peaks(1) = trisomy_peak(1);
    peaks(2) = trisomy_peak(2)*(1-fraction) + disomy_peak(2)*fraction;
    peaks(3) = trisomy_peak(3)*(1-fraction) + disomy_peak(2)*fraction;
    peaks(4) = trisomy_peak(4);
    for i = 1:length(peaks)-1
        actual_cutoffs(i)      = (peaks(i)+peaks(i+1))/2*200;
        mostLikelyGaussians(i) = i;
    end;
    mostLikelyGaussians = [mostLikelyGaussians length(peaks)];

    for i = 1:length(peaks)
        x_peaks(i) = peaks(i)*200;
    end;
    if (MakeFigure == true)
        for i = 1:4
            plot([x_peaks(i) x_peaks(i)],[0 max(raw)],'color',colorPeak,'linestyle','-','linewidth',2);
        end;
    end;        
    for i = 1:length(actual_cutoffs)
        if (MakeFigure == true)
            plot([actual_cutoffs(i) actual_cutoffs(i)],[0 max(raw)],'color',colorCutoff,'linestyle','-','linewidth',2);
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
    for i = 1:length(peaks)-1
        actual_cutoffs(i)      = (peaks(i)+peaks(i+1))/2*200;
        mostLikelyGaussians(i) = i;
    end;
    mostLikelyGaussians = [mostLikelyGaussians length(peaks)];
    
    for i = 1:length(peaks)
        x_peaks(i) = peaks(i)*200;
    end;
    if (MakeFigure == true)
        for i = 1:5
            plot([x_peaks(i) x_peaks(i)],[0 max(raw)],'color',colorPeak,'linestyle','-','linewidth',2);
        end;
    end;        
    for i = 1:length(actual_cutoffs)
        if (MakeFigure == true)
            plot([actual_cutoffs(i) actual_cutoffs(i)],[0 max(raw)],'color',colorCutoff,'linestyle','-','linewidth',2);
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
    for i = 1:length(peaks)-1
        actual_cutoffs(i)      = (peaks(i)+peaks(i+1))/2*200;
        mostLikelyGaussians(i) = i;
    end;
    mostLikelyGaussians = [mostLikelyGaussians length(peaks)];

    for i = 1:length(peaks)
        x_peaks(i) = peaks(i)*200;
    end;
    if (MakeFigure == true)
        for i = 1:6
            plot([x_peaks(i) x_peaks(i)],[0 max(raw)],'color',colorPeak,'linestyle','-','linewidth',2);
        end;
    end;        
    for i = 1:length(actual_cutoffs)
        if (MakeFigure == true)
                plot([actual_cutoffs(i) actual_cutoffs(i)],[0 max(raw)],'color',colorCutoff,'linestyle','-','linewidth',2);
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
    for i = 1:length(peaks)-1
        actual_cutoffs(i)      = (peaks(i)+peaks(i+1))/2*200;
        mostLikelyGaussians(i) = i;
    end;
    mostLikelyGaussians = [mostLikelyGaussians length(peaks)];

    for i = 1:length(peaks)
        x_peaks(i) = peaks(i)*200;
    end;
    if (MakeFigure == true)
        for i = 1:7
            plot([x_peaks(i) x_peaks(i)],[0 max(raw)],'color',colorPeak,'linestyle','-','linewidth',2);
        end;
    end;        
    for i = 1:length(actual_cutoffs)
        if (MakeFigure == true)
            plot([actual_cutoffs(i) actual_cutoffs(i)],[0 max(raw)],'color',colorCutoff,'linestyle','-','linewidth',2);
        end;
    end;
end;
if (MakeFigure == true)
    figure(fig);
    plot(raw     ,'color',[0.75 0.75 1.0],'linestyle','-','linewidth',1);
    plot(smoothed,'color',[0.50 0.50 1.0],'linestyle','-','linewidth',1);
    title([name '; chr ' num2str(chromosome) '; segment ' num2str(segment)],'HorizontalAlign','center','VerticalAlign','middle');
    hold off;
    xlim([1,200]);
    ylim([0,max(raw)]);
    % save then delete figures.

    % saveas(fig, [workingDir '50cutoffs_chr-' num2str(chromosome) '_seg-' num2str(segment) '.eps'], 'epsc');
      saveas(fig, [workingDir '50cutoffs_chr-' num2str(chromosome) '_seg-' num2str(segment) '.png'], 'png');
    delete(fig);
end;
end
