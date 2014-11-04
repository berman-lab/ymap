function [] = CNV_SNP_v1(parent_name,child_name,genome,ploidyEstimateString,ploidyBaseString,ave_copy_num, ...
			CNV_verString,SNP_verString,LOH_verString,workingDir,figureDir,displayBREAKS)
%% ========================================================================
% Generate graphs involved in analysis of CNV and SNP information from NGS
% datasets.
%==========================================================================
%    Centromere_format          : Controls how centromeres are depicted.   [0..2]   '2' is pinched cartoon default.
%    bases_per_bin              : Controls bin sizes for CGH fractions of plot.
%    scale_type                 : 'Ratio' or 'Log2Ratio' y-axis scaling of copy number.
%                                 'Log2Ratio' does not properly scale CGH data by ploidy.
%    Chr_max_width              : max width of chrs as fraction of figure width.
Centromere_format_default   = 0;
Yscale_nearest_even_ploidy  = true;
HistPlot                    = true;
ChrNum                      = true;
Chr_max_width               = 0.8;
show_annotations            = true;
analyze_rDNA                = true;
colorBars                   = true;
blendColorBars              = false;
Linear_display              = true;
Low_quality_ploidy_estimate = true;

%%=========================================================================
% Control variables.
%--------------------------------------------------------------------------
% Defines chr sizes in bp. (diploid total=28,567,7888)
% Defines centromere locations in bp.
% Defines annotation locations in bp.
[centromeres, chr_sizes, figure_details, annotations, ploidy_default] = Load_genome_information_1(workingDir,figureDir,genome);            
[Aneuploidy]                                                          = Load_dataset_information_1(child_name,workingDir);

for i = 1:length(chr_sizes)
    chr_size(chr_sizes(i).chr)    = chr_sizes(i).size;
end;
for i = 1:length(centromeres)
    cen_start(centromeres(i).chr) = centromeres(i).start;
    cen_end(centromeres(i).chr)   = centromeres(i).end;
end;
if (length(annotations) > 0)
    fprintf(['\nAnnotations for ' genome '.\n']);
    for i = 1:length(annotations)
        annotation_chr(i)       = annotations(i).chr;
        annotation_type{i}      = annotations(i).type;
        annotation_start(i)     = annotations(i).start;
        annotation_end(i)       = annotations(i).end;
        annotation_fillcolor{i} = annotations(i).fillcolor;
        annotation_edgecolor{i} = annotations(i).edgecolor;
        annotation_size(i)      = annotations(i).size;
        fprintf(['\t[' num2str(annotations(i).chr) ':' annotations(i).type ':' num2str(annotations(i).start) ':' num2str(annotations(i).end) ':' annotations(i).fillcolor ':' annotations(i).edgecolor ':' num2str(annotations(i).size) ']\n']);
    end;
end;
for i = 1:length(figure_details)
    if (figure_details(i).chr == 0)
        key_posX   = figure_details(i).posX;
        key_posY   = figure_details(i).posY;
        key_width  = figure_details(i).width;
        key_height = figure_details(i).height;
    else
        chr_id    (figure_details(i).chr) = figure_details(i).chr;
        chr_label {figure_details(i).chr} = figure_details(i).label;
        chr_name  {figure_details(i).chr} = figure_details(i).name;
        chr_posX  (figure_details(i).chr) = figure_details(i).posX;
        chr_posY  (figure_details(i).chr) = figure_details(i).posY;
        chr_width (figure_details(i).chr) = figure_details(i).width;
        chr_height(figure_details(i).chr) = figure_details(i).height;
    end;
end;
num_chr = length(chr_size);
bases_per_bin			= max(chr_size)/700;
chr_length_scale_multiplier	= 1/bases_per_bin;
CGD_bases_per_bin		= 1000;
CGD_chr_length_scale_multiplier	= 1/CGD_bases_per_bin;

%%=========================================================================
%%= No further control variables below. ===================================
%%=========================================================================

% Sanitize user input of euploid state.
ploidyBase = round(str2num(ploidyBaseString));
if (ploidyBase > 4);   ploidyBase = 4;   end;
if (ploidyBase < 1);   ploidyBase = 1;   end;
fprintf(['\nEuploid base = "' num2str(ploidyBase) '"\n']);

if (strcmp(parent_name,child_name) == 1)
    fprintf(['\nGenerating CNV-SNP map figure from ' parent_name ' genome data.\n']);
else
    fprintf(['\nGenerating CNV-LOH map figure from ' parent_name '(parent) and ' child_name '(child) genome data.\n']);
end;

%% Load CNV figure data from earlier analysis.
if (exist([workingDir 'matlab_dir/' child_name '.CNV_' CNV_verString '.ploidy_' ploidyEstimateString '.mat'],'file') == 0)
    fprintf('\nMAT(CNV) file not found, exiting.\n');
    exit;
else
    fprintf('\nMAT(CNV) file found, loading.\n');
    load([workingDir 'matlab_dir/' child_name '.CNV_' CNV_verString '.ploidy_' ploidyEstimateString '.mat']);
    % 'chr_CNVdata','CGD_chr_CNVdata'.
end;

%% Load SNP/LOH figure data from earlier analysis.
% file names differ if the SNPs of one or two strains are being examined.
if (strcmp(parent_name,child_name) == 1)
    fileName_SNP = [child_name '.SNP_' SNP_verString '.mat'];
else
    fileName_SNP = [parent_name '->' child_name '.LOH_' LOH_verString '.mat'];
end;
if (exist([workingDir 'matlab_dir/' fileName_SNP],'file') == 0)
    fprintf('\nMAT(SNP) file not found, exiting.\n');
    exit;
else
    fprintf('\nMAT(SNP) file found, loading.\n');
    load([workingDir 'matlab_dir/' fileName_SNP]);
    % 'chr_SNPdata','ave_copy_number'.
end;

%% generate plots.
% basic plot parameters not defined per genome.
TickSize        = -0.005;  %negative for outside, percentage of longest chr figure.
maxY            = ploidyBase*2;

%define colors for colorBars plot
colorNoData = [1.0   1.0   1.0  ]; %used when no data is available for the bin.
colorInit   = [0.5   0.5   0.5  ]; %external; used in blending at ends of chr.
colorHET    = [0.0   0.0   0.0  ]; % near 1:1 ratio SNPs
colorOddHET = [0.0   1.0   0.0  ]; % Het, but not near 1:1 ratio SNPs.
colorHOM    = [1.0   0.0   0.0  ]; % Hom SNPs;

%% -----------------------------------------------------------------------------------------
% Setup for main figure generation.
%-------------------------------------------------------------------------------------------
fig = figure(1);
set(gcf, 'Position', [0 70 1024 600]);

%% CNV pre-figure calculations.
for chr = 1:num_chr
    CNVplot{chr} = chr_CNVdata{chr};
    chr_max(chr) = max(CNVplot{chr});
    chr_med(chr) = median(CNVplot{chr});
end;
for chr = 1:num_chr
    max_count     = max(chr_max);
    median_count  = sum(chr_med)/length(chr_med);
    CNVplot2{chr} = CNVplot{chr}/median_count;
end;
ploidy = str2num(ploidyEstimateString);
[chr_breaks, chrCopyNum, ploidyAdjust] = FindChrSizes_2(Aneuploidy,CNVplot2,ploidy,num_chr);
fprintf('\n');

%% SNP pre-figure calculations.
data_mode = 3;
for chr = 1:num_chr
    if (data_mode == 1)
        % Regenerate chr plot data if the save file does not exist.
        TOTplot{chr}                           = chr_SNPdata{chr,1}+chr_SNPdata{chr,2}+chr_SNPdata{chr,3};  % TOT data
        TOTave{chr}                            = sum(TOTplot{chr})/length(TOTplot{chr});
        TOTplot2{chr}                          = TOTplot{chr}/ave_copy_num;
        TOTplot2{chr}(TOTplot2{chr} > 1)       = 1;
        TOTave2{chr}                           = sum(TOTplot2{chr})/length(TOTplot2{chr});
                
        HETplot{chr}                           = chr_SNPdata{chr,1};  % HET data
        HETave{chr}                            = sum(HETplot{chr})/length(HETplot{chr});
        HETplot2{chr}                          = HETplot{chr}/ave_copy_num;
        HETplot2{chr}(HETplot2{chr} > 1)       = 1;
        HETave2{chr}                           = sum(HETplot2{chr})/length(HETplot2{chr});
        
        oddHETplot{chr}                        = chr_SNPdata{chr,2};  % oddHET data
        oddHETave{chr}                         = sum(oddHETplot{chr})/length(oddHETplot{chr});
        oddHETplot2{chr}                       = oddHETplot{chr}/ave_copy_num;
        oddHETplot2{chr}(oddHETplot2{chr} > 1) = 1;
        oddHETave2{chr}                        = sum(oddHETplot2{chr})/length(oddHETplot2{chr});
    
        HOMplot{chr}                           = chr_SNPdata{chr,3};  % HOM data
        HOMave{chr}                            = sum(HOMplot{chr})/length(HOMplot{chr});
        HOMplot2{chr}                          = HOMplot{chr}/ave_copy_num;
        HOMplot2{chr}(HOMplot2{chr} > 1)       = 1;
        HOMave2{chr}                           = sum(HOMplot2{chr})/length(HOMplot2{chr});
    elseif (data_mode == 2)
        %% Details from LOH_v2a.m :
        % Regenerate chr plot data if the save file does not exist.
        TOTplot{chr}                                  = chr_SNPdata{chr,1}+chr_SNPdata{chr,2}+chr_SNPdata{chr,3};  % TOT data
        HETplot{chr}                                  = chr_SNPdata{chr,1};  % HET data
        oddHETplot{chr}                               = chr_SNPdata{chr,2};  % oddHET data
        HOMplot{chr}                                  = chr_SNPdata{chr,3};  % HOM data 
        
        TOTave{chr}                                   = sum(TOTplot{chr})/length(TOTplot{chr});
        TOTplot2{chr}                                 = TOTplot{chr}/ave_copy_num;
        TOTplot2{chr}(TOTplot2{chr} > 1)              = 1;
        TOTave2{chr}                                  = sum(TOTplot2{chr})/length(TOTplot2{chr});
        
        HETave{chr}                                   = sum(HETplot{chr})/length(HETplot{chr});
        HETplot2{chr}                                 = HETplot{chr}/ave_copy_num;
        HETplot2{chr}(HETplot2{chr} > 1)              = 1;
        HETave2{chr}                                  = sum(HETplot2{chr})/length(HETplot2{chr});
        oddHETave{chr}                                = sum(oddHETplot{chr})/length(oddHETplot{chr});
        oddHETplot2{chr}                              = oddHETplot{chr}/ave_copy_num;   
        oddHETplot2{chr}(oddHETplot2{chr} > 1)        = 1;
        oddHETave2{chr}                               = sum(oddHETplot2{chr})/length(oddHETplot2{chr});
        HOMave{chr}                                   = sum(HOMplot{chr})/length(HOMplot{chr});
        HOMplot2{chr}                                 = HOMplot{chr}/ave_copy_num;
        HOMplot2{chr}(HOMplot2{chr} > 1)              = 1;
        HOMave2{chr}                                  = sum(HOMplot2{chr})/length(HOMplot2{chr});
    elseif (data_mode == 3)
        %% Details from LOH_v3a.m :
        % Regenerate chr plot data if the save file does not exist.
        TOTplot{chr}                                  = chr_SNPdata{chr,1}+chr_SNPdata{chr,2}+chr_SNPdata{chr,3};  % TOT data
        TOTave{chr}                                   = sum(TOTplot{chr})/length(TOTplot{chr});
        TOTplot2{chr}                                 = TOTplot{chr}/ave_copy_num;
        TOTplot2{chr}(TOTplot2{chr} > 1)              = 1;
        TOTave2{chr}                                  = sum(TOTplot2{chr})/length(TOTplot2{chr});
        
        HETplot{chr}                                  = chr_SNPdata{chr,1};  % HET data
        HETave{chr}                                   = sum(HETplot{chr})/length(HETplot{chr});
        HETplot2{chr}                                 = HETplot{chr}/ave_copy_num;
        HETplot2{chr}(HETplot2{chr} > 1)              = 1;
        HETave2{chr}                                  = sum(HETplot2{chr})/length(HETplot2{chr});
        
        oddHETplot{chr}                               = chr_SNPdata{chr,2};  % oddHET data
        oddHETave{chr}                                = sum(oddHETplot{chr})/length(oddHETplot{chr});
        oddHETplot2{chr}                              = oddHETplot{chr}/ave_copy_num;
        oddHETplot2{chr}(oddHETplot2{chr} > 1)        = 1;
        
        HOMplot{chr}                                  = chr_SNPdata{chr,3};  % HOM data
        HOMave{chr}                                   = sum(HOMplot{chr})/length(HOMplot{chr});  
        HOMplot2{chr}                                 = HOMplot{chr}/ave_copy_num;
        HOMplot2{chr}(HOMplot2{chr} > 1)              = 1;
    end;
end;
fprintf('\n');
largestChr = find(chr_width == max(chr_width));

%% -----------------------------------------------------------------------------------------
% Setup for linear-view figure generation.
%-------------------------------------------------------------------------------------------
if (Linear_display == true)
    Linear_fig = figure(2);
    Linear_genome_size     = sum(chr_size);
    Linear_Chr_max_width   = 0.85;
    Linear_left_start      = 0.07;
    Linear_left_chr_gap    = 0.01;
    Linear_height          = 0.6;
    Linear_base            = 0.1;
    Linear_TickSize        = -0.01;  %negative for outside, percentage of longest chr figure.
    Linear_maxY            = ploidyBase*2;

    Linear_left = Linear_left_start;
end;

%% -----------------------------------------------------------------------------------------
% Make figures
%-------------------------------------------------------------------------------------------
for chr = 1:num_chr
    figure(fig);
    % make standard chr cartoons.
    left   = chr_posX(chr);
    bottom = chr_posY(chr);
    width  = chr_width(chr);
    height = chr_height(chr);
    subplot('Position',[left bottom width height]);
    fprintf(['figposition = [' num2str(left) ' | ' num2str(bottom) ' | ' num2str(width) ' | ' num2str(height) ']\t']);
    hold on;

    %% SNP(LOH) plot section.
    c_prev = colorInit;
    c_post = colorInit;
    c_     = c_prev;

    infill = zeros(1,length(HETplot2{chr}));
    colors = [];
    
    % determines the color of each bin.
    for i = 1:length(TOTplot2{chr})+1;
        if (i-1 < length(TOTplot2{chr}))
            c_tot_post = TOTplot2{chr}(i)+TOTplot2{chr}(i);
            if (c_tot_post == 0)
                c_post = colorNoData;
            else
                %c_post = colorHET*HETplot2{chr}(i) + ...
                %         colorHOM*HOMplot2{chr}(i) + ...
                %         colorNoData*(1-min([HETplot2{chr}(i)+HOMplot2{chr}(i) 1]));
                %colorMix = colorHET   *HETplot2   {chr}(i)/TOTplot2{chr}(i) + ...
                %           colorOddHET*oddHETplot2{chr}(i)/TOTplot2{chr}(i) + ...
                %           colorHOM   *HOMplot2   {chr}(i)/TOTplot2{chr}(i);
                colorMix = colorHET   *   HETplot2{chr}(i)/TOTplot2{chr}(i) + ...
                           colorOddHET*oddHETplot2{chr}(i)/TOTplot2{chr}(i) + ...
                           colorHOM   *   HOMplot2{chr}(i)/TOTplot2{chr}(i);
                c_post =   colorMix   *   min(1,TOTplot2{chr}(i)) + ...
                           colorNoData*(1-min(1,TOTplot2{chr}(i)));
                %colorNoData*(1-min([HETplot2{chr}(i)+oddHETplot2{chr}(i)+HOMplot2{chr}(i) 1]));
            end;
        else
            c_post = colorInit;
        end;
        colors(i,1) = c_post(1);
        colors(i,2) = c_post(2);
        colors(i,3) = c_post(3);
    end;
    % draw colorbars.
    for i = 1:length(HETplot2{chr})+1;
        x_ = [i i i-1 i-1];
        y_ = [0 maxY maxY 0];
        c_post(1) = colors(i,1);
        c_post(2) = colors(i,2);
        c_post(3) = colors(i,3);
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
    
    %% cgh plot section.
    c_ = [0 0 0];
    fprintf(['main-plot : chr' num2str(chr) ':' num2str(length(CNVplot2{chr})) '\n']);
    for i = 1:length(CNVplot2{chr});
        x_ = [i i i-1 i-1];
        if (CNVplot2{chr}(i) == 0)
            CNVhistValue = 1;
        else
            CNVhistValue = CNVplot2{chr}(i);
        end;

        % The CNV-histogram values were normalized to a median value of 1.
        % The ratio of 'ploidy' to 'ploidyBase' determines where the data is displayed relative to the median line.
        startY = maxY/2;
        if (Low_quality_ploidy_estimate == true)
            endY = CNVhistValue*ploidy*ploidyAdjust;
        else
            endY = CNVhistValue*ploidy;
        end;
        y_ = [startY endY endY startY];

        % makes a blackbar for each bin.
        f = fill(x_,y_,c_);
        set(f,'linestyle','none');
    end;
    x2 = chr_size(chr)*chr_length_scale_multiplier;
    plot([0; x2], [maxY/2; maxY/2],'color',[0 0 0]);  % 2n line.
    
    %% draw lines across plots for easier interpretation of CNV regions.
    switch ploidyBase
        case 1
        case 2
            line([0 x2], [maxY/4*1 maxY/4*1],'Color',[0.85 0.85 0.85]);
            line([0 x2], [maxY/4*3 maxY/4*3],'Color',[0.85 0.85 0.85]);
        case 3
            line([0 x2], [maxY/6*1 maxY/6*1],'Color',[0.85 0.85 0.85]);
            line([0 x2], [maxY/6*2 maxY/6*2],'Color',[0.85 0.85 0.85]);
            line([0 x2], [maxY/6*4 maxY/6*4],'Color',[0.85 0.85 0.85]);
            line([0 x2], [maxY/6*5 maxY/6*5],'Color',[0.85 0.85 0.85]);
        case 4
            line([0 x2], [maxY/8*1 maxY/8*1],'Color',[0.85 0.85 0.85]);
            line([0 x2], [maxY/8*2 maxY/8*2],'Color',[0.85 0.85 0.85]);
            line([0 x2], [maxY/8*3 maxY/8*3],'Color',[0.85 0.85 0.85]);
            line([0 x2], [maxY/8*5 maxY/8*5],'Color',[0.85 0.85 0.85]);
            line([0 x2], [maxY/8*6 maxY/8*6],'Color',[0.85 0.85 0.85]);
            line([0 x2], [maxY/8*7 maxY/8*7],'Color',[0.85 0.85 0.85]);
    end;
    %% end cgh plot section.
    
    %axes labels etc.
    hold off;
    % limit x-axis to range of chromosome.
    xlim([0,chr_size(chr)*chr_length_scale_multiplier]);
    
    %% modify y axis limits to show annotation locations if any are provided.
    if (length(annotations) > 0)
        ylim([-maxY/10*1.5,maxY]);
    else
        ylim([0,maxY]);
    end;
    
    set(gca,'YTick',[]);
    set(gca,'TickLength',[(TickSize*chr_size(largestChr)/chr_size(chr)) 0]); %ensures same tick size on all subfigs.
    ylabel(chr_label{chr}, 'Rotation', 90, 'HorizontalAlign', 'center', 'VerticalAlign', 'bottom');
    set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
    set(gca,'XTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8','2.0','2.2','2.4','2.6','2.8','3.0','3.2'});

    % This section sets the Y-axis labelling.
    switch ploidyBase
        case 1
            set(gca,'YTick',[0 maxY/2 maxY]);
            set(gca,'YTickLabel',{'','1','2'});
        case 2
            set(gca,'YTick',[0 maxY/4 maxY/2 maxY/4*3 maxY]);
            set(gca,'YTickLabel',{'','1','2','3','4'});
        case 3
            set(gca,'YTick',[0 maxY/6 maxY/3 maxY/2 maxY/3*2 maxY/6*5 maxY]);
            set(gca,'YTickLabel',{'','','','3','','','6'});
        case 4
            set(gca,'YTick',[0 maxY/8 maxY/4 maxY/8*3 maxY/2 maxY/8*5 maxY/4*3 maxY/8*7 maxY]);
            set(gca,'YTickLabel',{'','','2','','4','','6','','8'});
    end;

    set(gca,'FontSize',6);
    if (chr == find(chr_posY == max(chr_posY)))
        title([ child_name ' CNV map'],'Interpreter','none','FontSize',12);
    end;
    
    hold on;
    %end axes labels etc.
    
    %show segmental anueploidy breakpoints.
    if (displayBREAKS == true)
        for segment = 2:length(chr_breaks{chr})-1
            bP = chr_breaks{chr}(segment)*length(CNVplot2{chr});
            c_ = [0 0 1];
            x_ = [bP bP bP-1 bP-1];
            y_ = [0 maxY maxY 0];
            f = fill(x_,y_,c_);
            set(f,'linestyle','none');
        end;
    end;
    
    %show centromere.
    if (chr_size(chr) < 100000)
        Centromere_format = 1;
    else
        Centromere_format = Centromere_format_default;
    end;
    x1 = cen_start(chr)*chr_length_scale_multiplier;
    x2 = cen_end(chr)*chr_length_scale_multiplier;
    leftEnd  = 0.5*(5000/bases_per_bin);
    rightEnd = chr_size(chr)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
    if (Centromere_format == 0)
        % standard chromosome cartoons in a way which will not cause segfaults when running via commandline.
        dx     = 5*(5000/bases_per_bin);
        dy     = maxY/5;
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
    elseif (Centromere_format == 1)
        leftEnd  = 0;
        rightEnd = chr_size(chr)*chr_length_scale_multiplier;
        
        % Minimal outline for examining very small sequence regions, such as C.albicans MTL locus.
        plot([leftEnd   leftEnd   rightEnd   rightEnd   leftEnd], [0   maxY   maxY   0   0], 'Color',[0 0 0]);
    end;
    %end show centromere.
    
    %show annotation locations
    if (show_annotations) && (length(annotations) > 0)
        plot([leftEnd rightEnd], [-maxY/10*1.5 -maxY/10*1.5],'color',[0 0 0]);
        hold on;
        annotation_location = (annotation_start+annotation_end)./2;
        for i = 1:length(annotation_location)
            if (annotation_chr(i) == chr)
                annotationloc = annotation_location(i)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
                plot(annotationloc,-maxY/10*1.5,'k:o','MarkerEdgeColor',annotation_edgecolor{i}, ...
                                                      'MarkerFaceColor',annotation_fillcolor{i}, ...
                                                      'MarkerSize',     annotation_size(i));
            end;
        end;
        hold off;
    end;
    %end show annotation locations.
    
    % make CGH histograms to the right of the main chr cartoons.
    if (HistPlot == true)
        width     = 0.020;
        height    = chr_height(chr);
        bottom    = chr_posY(chr);
        chr_CNVdata;
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
            histAll{segment}(histAll{segment}<=0)             = [];
            histAll{segment}(histAll{segment}>ploidyBase*2+2) = ploidyBase*2+2;
            histAll{segment}(length(histAll{segment})+1)      = 0;   % endpoints added to ensure histogram bounds.
            histAll{segment}(length(histAll{segment})+1)      = ploidyBase*2+2;
            smoothed{segment}    = smooth_gaussian(hist(histAll{segment},(ploidyBase*2+2)*50),5,20);
            % make a smoothed version of just the endpoints used to ensure histogram bounds.
            histAll2{segment}(1) = 0;
            histAll2{segment}(2) = ploidyBase*2+2;
            smoothed2{segment}   = smooth_gaussian(hist(histAll2{segment},(ploidyBase*2+2)*50),5,20)*4;
            % subtract the smoothed endpoints from the histogram to remove the influence of the added endpoints.
            smoothed{segment}    = (smoothed{segment}-smoothed2{segment});
            smoothed{segment}    = smoothed{segment}/max(smoothed{segment});

            hold on;
            for i = 1:(ploidyBase*2-1)
                plot([0; 1],[i*50; i*50],'color',[0.75 0.75 0.75]);
            end;
            area(smoothed{segment},1:length(smoothed{segment}),'FaceColor',[0 0 0]);
            hold off;
            set(gca,'YTick',[]);
            set(gca,'XTick',[]);
            xlim([0,1]);
	    ylim([0,(ploidyBase*2)*50]);
        end;
    end;
    
    % places chr copy number to the right of the main chr cartoons.
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
        if (length(chrCopyNum{chr}) == 1)
            chr_string = num2str(chrCopyNum{chr}(1));
        else%		chr_string = num2str(chrCopyNum{chr}(1));
            for i = 2:length(chrCopyNum{chr})
                chr_string = [chr_string ',' num2str(chrCopyNum{chr}(i))];
            end;
        end;
        text(0.1,0.5, chr_string,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',12);
    end;

    %% Linear figure draw section
    if (Linear_display == true)
        figure(Linear_fig);  
        Linear_width = Linear_Chr_max_width*chr_size(chr)/Linear_genome_size;
        subplot('Position',[Linear_left Linear_base Linear_width Linear_height]);
        Linear_left = Linear_left + Linear_width + Linear_left_chr_gap;
        hold on;
        title(chr_label{chr},'Interpreter','none','FontSize',10);

        % draw colorbars.
        for i = 1:length(HETplot2{chr})+1;
            x_ = [i i i-1 i-1];
            y_ = [0 maxY maxY 0];
            c_post(1) = colors(i,1);
            c_post(2) = colors(i,2);
            c_post(3) = colors(i,3);
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

        %% cgh plot section.
        c_ = [0 0 0];
        fprintf(['linear-plot : chr' num2str(chr) ':' num2str(length(CNVplot2{chr})) '\n']);
        for i = 1:length(CNVplot2{chr}); 
            x_ = [i i i-1 i-1];
            if (CNVplot2{chr}(i) == 0)
                CNVhistValue = 1;
            else
                CNVhistValue = CNVplot2{chr}(i);
            end;
    
            % The CNV-histogram values were normalized to a median value of 1.
            % The ratio of 'ploidy' to 'ploidyBase' determines where the data is displayed relative to the median line.
            startY = maxY/2;
            if (Low_quality_ploidy_estimate == true)
                endY = CNVhistValue*ploidy*ploidyAdjust;
            else
                endY = CNVhistValue*ploidy;
            end;
            y_ = [startY endY endY startY];

            % makes a blackbar for each bin.
            f = fill(x_,y_,c_);
            set(f,'linestyle','none');
        end;
        x2 = chr_size(chr)*chr_length_scale_multiplier;
        plot([0; x2], [maxY/2; maxY/2],'color',[0 0 0]);  % 2n line.

        %% draw lines across plots for easier interpretation of CNV regions.
        switch ploidyBase
            case 1
            case 2
                line([0 x2], [maxY/4*1 maxY/4*1],'Color',[0.85 0.85 0.85]);
                line([0 x2], [maxY/4*3 maxY/4*3],'Color',[0.85 0.85 0.85]);
            case 3
                line([0 x2], [maxY/6*1 maxY/6*1],'Color',[0.85 0.85 0.85]);
                line([0 x2], [maxY/6*2 maxY/6*2],'Color',[0.85 0.85 0.85]);
                line([0 x2], [maxY/6*4 maxY/6*4],'Color',[0.85 0.85 0.85]);
                line([0 x2], [maxY/6*5 maxY/6*5],'Color',[0.85 0.85 0.85]);
            case 4
                line([0 x2], [maxY/8*1 maxY/8*1],'Color',[0.85 0.85 0.85]);
                line([0 x2], [maxY/8*2 maxY/8*2],'Color',[0.85 0.85 0.85]);
                line([0 x2], [maxY/8*3 maxY/8*3],'Color',[0.85 0.85 0.85]);
                line([0 x2], [maxY/8*5 maxY/8*5],'Color',[0.85 0.85 0.85]);
                line([0 x2], [maxY/8*6 maxY/8*6],'Color',[0.85 0.85 0.85]);
                line([0 x2], [maxY/8*7 maxY/8*7],'Color',[0.85 0.85 0.85]);
        end;
        %% end cgh plot section.
                    
        %show segmental anueploidy breakpoints.
        if (displayBREAKS == true)
            for segment = 2:length(chr_breaks{chr})-1
                bP = chr_breaks{chr}(segment)*length(CNVplot2{chr});
                c_ = [0 0 1];
                x_ = [bP bP bP-1 bP-1];
                y_ = [0 maxY maxY 0];
                f = fill(x_,y_,c_);   
                set(f,'linestyle','none');
            end;
        end;

        %show centromere.
        if (chr_size(chr) < 100000)
            Centromere_format = 1;
        else
            Centromere_format = Centromere_format_default;
        end;
        x1 = cen_start(chr)*chr_length_scale_multiplier;
        x2 = cen_end(chr)*chr_length_scale_multiplier;
        leftEnd  = 0.5*(5000/bases_per_bin);
        rightEnd = chr_size(chr)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
        if (Centromere_format == 0)
            % standard chromosome cartoons in a way which will not cause segfaults when running via commandline.
            dx     = 5*(5000/bases_per_bin);
            dy     = maxY/10;
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
                 [dy        maxY-dy   maxY         maxY    maxY-dy   maxY-dy   maxY    maxY          maxY-dy    dy         0             0       dy   dy   0       0            dy],...
                 'Color',[0 0 0]);
        elseif (Centromere_format == 1)
            leftEnd  = 0;
            rightEnd = chr_size(chr)*chr_length_scale_multiplier;
            
            % Minimal outline for examining very small sequence regions, such as C.albicans MTL locus.
            plot([leftEnd   leftEnd   rightEnd   rightEnd   leftEnd], [0   maxY   maxY   0   0], 'Color',[0 0 0]);
        end;
        %end show centromere.

        %show annotation locations
        if (show_annotations) && (length(annotations) > 0)
            plot([leftEnd rightEnd], [-maxY/10*1.5 -maxY/10*1.5],'color',[0 0 0]);
            hold on;
            annotation_location = (annotation_start+annotation_end)./2;
            for i = 1:length(annotation_location)
                if (annotation_chr(i) == chr)
                    annotationloc = annotation_location(i)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
                    plot(annotationloc,-maxY/10*1.5,'k:o','MarkerEdgeColor',annotation_edgecolor{i}, ...
                                                          'MarkerFaceColor',annotation_fillcolor{i}, ...
                                                          'MarkerSize',annotation_size(i));
                end;
            end;
            hold off;
        end;
        %end show annotation locations.

        %% Final formatting stuff.
        xlim([0,chr_size(chr)*chr_length_scale_multiplier]);
        % modify y axis limits to show annotation locations if any are provided.
        if (length(annotations) > 0)
            ylim([-maxY/10*1.5,maxY]);
        else
            ylim([0,maxY]);
        end;
        set(gca,'TickLength',[(Linear_TickSize*chr_size(1)/chr_size(chr)) 0]); %ensures same tick size on all subfigs.
        set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
        set(gca,'XTickLabel',[]);
        if (chr == 1)
            if (strcmp(parent_name,child_name) == 1)
                ylabel(parent_name,                 ...
                       'Rotation',        0,        ...
                       'HorizontalAlign', 'right',  ...
                       'VerticalAlign',   'bottom', ...
                       'Interpreter',     'none',   ...
                       'FontSize',        5);
            else
                ylabel({parent_name;'vs.';child_name}, ...
                       'Rotation',        0,           ...
                       'HorizontalAlign', 'right',     ...
                       'VerticalAlign',   'bottom',    ...
                       'Interpreter',     'none',      ...
                       'FontSize',        5);
            end;

            % This section sets the Y-axis labelling.
            switch ploidyBase
                case 1
                    set(gca,'YTick',[0 maxY/2 maxY]);
                    set(gca,'YTickLabel',{'','1','2'});
                case 2
                    set(gca,'YTick',[0 maxY/4 maxY/2 maxY/4*3 maxY]);
                    set(gca,'YTickLabel',{'','1','2','3','4'});
                case 3
                    set(gca,'YTick',[0 maxY/6 maxY/3 maxY/2 maxY/3*2 maxY/6*5 maxY]);
                    set(gca,'YTickLabel',{'','','','3','','','6'});
                case 4
                    set(gca,'YTick',[0 maxY/8 maxY/4 maxY/8*3 maxY/2 maxY/8*5 maxY/4*3 maxY/8*7 maxY]);
                    set(gca,'YTickLabel',{'','','2','','4','','6','','8'});
            end;
        else
            set(gca,'YTick',[]);
            set(gca,'YTickLabel',[]);
        end;
        set(gca,'FontSize',6);
        %end final reformatting.
        
        % shift back to main figure generation.
        figure(fig);
        hold on;
    end;
end;

%	% Main figure colors key.
%	left   = key_posX;
%	bottom = key_posY;
%	width  = str2num(key_width);
%	height = key_height;
%	colorNoData = [1.0   1.0   1.0  ]; %used when no data is available for the bin.
%	colorCNV    = [0.0   0.0   0.0  ]; %external; used in blending at ends of chr.

%	subplot('Position',[left bottom width height]);
%	axis off square;
%	xlim([-0.1,1]);
%	ylim([-0.1,1.6]);
%	set(gca,'XTick',[]);
%	set(gca,'YTick',[]);
%	patch([0 0.2 0.2 0], [1.4 1.4 1.5 1.5], colorCNV);      text(0.3,1.05,'Copy number variation (CNV).');
%	patch([0 0.2 0.2 0], [1.0 1.0 1.1 1.1], colorNoData);   text(0.3,1.45,'No CNV.');

if (strcmp(parent_name,child_name) == 1)
    saveas(fig, [figureDir parent_name '.CNV-SNP-map.1.eps'], 'epsc');
    set(Linear_fig,'PaperPosition',[0 0 8 0.62222222]);
    saveas(Linear_fig, [figureDir parent_name '.CNV-SNP-map.2.eps'], 'epsc');
else
    saveas(fig, [figureDir parent_name '->' child_name '.CNV-LOH-map.1.eps'], 'epsc');
    set(Linear_fig,'PaperPosition',[0 0 8 0.62222222]);
    saveas(fig, [figureDir parent_name '->' child_name '.CNV-LOH-map.2.eps'], 'epsc');
end;
delete(fig);
delete(Linear_fig);
end
