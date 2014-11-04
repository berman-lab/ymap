function [] = LOH_v2_1(parent_name,child_name,genome,ave_copy_num,SNP_verString,LOH_verString,workingDir,figureDir,displayBREAKS)
%% ========================================================================
%    Centromere_format          : Controls how centromeres are depicted.   [0..2]   '2' is pinched cartoon default.
%    bases_per_bin              : Controls bin sizes for SNP/CGH fractions of plot.
%    scale_type                 : 'Ratio' or 'Log2Ratio' y-axis scaling of copy number.
%                                 'Log2Ratio' does not properly scale CGH data by ploidy.
%    Chr_max_width              : max width of chrs as fraction of figure width.
Centromere_format           = 0;
Chr_max_width               = 0.8;
colorBars                   = true;
blendColorBars              = false;
show_annotations            = true;
Yscale_nearest_even_ploidy  = true;
Linear_display              = true;

%%=========================================================================
% Control variables for Candida albicans SC5314.
%--------------------------------------------------------------------------
% Defines chr sizes in bp. (diploid total=28,567,7888)
% Defines centromere locations in bp.
% Defines MRS locations in bp.
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

%% This block is normally calculated in FindChrSizes_2 in CNV analysis.
for usedChr = 1:num_chr
    % determine where the endpoints of ploidy segments are.
    chr_breaks{usedChr}(1) = 0.0;
    break_count = 1;
    if (length(Aneuploidy) > 0)
        for i = 1:length(Aneuploidy)
            if (Aneuploidy(i).dataset == dataset) && (Aneuploidy(i).chr == usedChr)
                break_count = break_count+1;
                chr_broken = true;
                chr_breaks{usedChr}(break_count) = Aneuploidy(i).break;
            end;
        end;
    end;
    chr_breaks{usedChr}(length(chr_breaks{usedChr})+1) = 1;
end;

% bases_per_bin			= 5000;
bases_per_bin			= max(chr_size)/700;
chr_length_scale_multiplier	= 1/bases_per_bin;
CGD_bases_per_bin		= 1000;
CGD_chr_length_scale_multiplier	= 1/CGD_bases_per_bin;

%%=========================================================================
%%= No further control variables below. ===================================
%%=========================================================================
if (strcmp(parent_name,child_name) == 1)
    fprintf(['\nGenerating SNP-map figure from ' parent_name ' genome data.\n']);
else
    fprintf(['\nGenerating LOH-map figure from ' parent_name ' and ' child_name ' genome data.\n']);
end;


%%================================================================================================
% Load SNP/LOH data.
%-------------------------------------------------------------------------------------------------
if (strcmp(parent_name,child_name) == 1)
    % if only one strain is examined, then file name should indicate SNP analysis.
    fileName = [parent_name '.SNP_' SNP_verString '.mat'];
else
    % if two strains are examined, then file name should indicate LOH analysis.
    fileName = [parent_name '->' child_name '.LOH_' LOH_verString '.mat'];
end;

if (exist([workingDir 'matlab_dir/' fileName],'file') == 0)
    % Preprocess "[*****]_putatitve_SNPs.txt" files into separate data structures for
    %    homozygous (aa), heterozygous (ab), & homozygous (bb) potential SNPs.

    % Preprocess 'parent' dataset for SNP categories if they have not already been processed.
    if (exist([workingDir 'matlab_dir/' parent_name '.SNP_' SNP_verString '.het.mat'],'file') == 0)
        preprocess_SNP_data_files_1(parent_name,SNP_verString,workingDir);
    end;

    % Preprocess 'child' dataset for SNP categories if they have not already been processed.
    if (exist([workingDir 'matlab_dir/' child_name '.SNP_' SNP_verString '.het.mat'],'file') == 0) && ...
        (exist([workingDir 'matlab_dir/' child_name '.SNP_' SNP_verString '.hom1.mat'],'file') == 0) && ...
            (exist([workingDir 'matlab_dir/' child_name '.SNP_' SNP_verString '.hom2.mat'],'file') == 0)
        preprocess_SNP_data_files_1(child_name,SNP_verString,workingDir);
    end;

    % If 'parent' and 'child' are the same, load the 'child' dataset first.
    % If 'parent' and 'child' are different, load the 'parent' dataset first.
    if (strcmp(parent_name,child_name) == 1)
        % Load 'child' homozygous (aa), heterozygous (ab), & homozygous (bb) datasets.
        if (strcmp(parent_name,child_name) == 1)
            fprintf('Loading data for dataset :\n');
        else
            fprintf('Loading data for (child) dataset :\n');
        end;
        fprintf('\tLoading het SNP data structure.\n');
        load([workingDir 'matlab_dir/' child_name '.SNP_' SNP_verString '.het.mat']);
        fprintf('\tLoading hom SNP data structure. [1/2]\n');
        load([workingDir 'matlab_dir/' child_name '.SNP_' SNP_verString '.hom1.mat']);
        fprintf('\tLoading hom SNP data structure. [2/2]\n');
        load([workingDir 'matlab_dir/' child_name '.SNP_' SNP_verString '.hom2.mat']);

        % Assign 'parent' heterozygous data equal to 'child'.
        parent_het = dataset_het;
    else
        % Load 'parent' heterozygous (ab) dataset.
        fprintf('Loading data for (parent) dataset :\n');
        fprintf(['\tLoading het SNP data for ' parent_name '.\n']);
        load([workingDir 'matlab_dir/' parent_name '.SNP_' SNP_verString '.het.mat']);
        parent_het = dataset_het;
        
        % Load 'child' homozygous (aa), heterozygous (ab), & homozygous (bb) datasets.
        if (strcmp(parent_name,child_name) == 1)
            fprintf('Loading data for dataset :\n');
        else
            fprintf('Loading data for (child) dataset :\n');
        end;
        fprintf('\tLoading het SNP data structure.\n');
        load([workingDir 'matlab_dir/' child_name '.SNP_' SNP_verString '.het.mat']);
        fprintf('\tLoading hom SNP data structure. [1/2]\n');
        load([workingDir 'matlab_dir/' child_name '.SNP_' SNP_verString '.hom1.mat']);
        fprintf('\tLoading hom SNP data structure. [2/2]\n');
        load([workingDir 'matlab_dir/' child_name '.SNP_' SNP_verString '.hom2.mat']);
    end;
    fprintf('\tSNP data structures loaded.\n');
    
    % Initializes vectors used to hold number of SNPs in each interpretation catagory
    %     for each chromosome region.
    for chr = 1:length(chr_sizes)
        % 3 SNP interpretation catagories tracked.
        for j = 1:3
            chr_SNPdata{chr,j} = zeros(1,ceil(chr_size(chr)/bases_per_bin));
            % 1 = heterozygous, 1:1.
            % 2 = heterozygous, not 1:1.
            % 3 = homozygous.
        end;
    end;
    
    % tally het data into bins for main figure histogram.
    HET_parent_length   = length(parent_het);
    HET_dataset_length  = length(dataset_het);
    HOM1_dataset_length = length(dataset_hom1);
    HOM2_dataset_length = length(dataset_hom2);
    HOM_dataset_length  = length(dataset_hom1)+length(dataset_hom2);

    % Tracking counters.
    counter = 0;
    het_place = 1;
    hom_place = 1;

    fprintf('\n\tProcessing het data.');
    fprintf('\n\t\t[');
    for i = 1:HET_parent_length
        counter = counter+1;
        % locus is the genomic location of the SNP.
        locus  = ceil(str2double((parent_het{i}.position))/bases_per_bin);
        
        % Determines distribution of data.
        parent_chr = 0;
        for ii = 1:length(figure_details)
            if (figure_details(ii).chr ~= 0)
                if (strcmp(parent_het{i}.chr,chr_name{ii}) == 1)
                    parent_chr = ii;
                end;
            end;
        end;

        if (parent_chr > 0)
            % Initially define child SNP state as 'null'.
            child = 'null';
            % If parent SNP is found in het dataset for child, define child SNP state as 'HET'.
            for j1 = het_place:HET_dataset_length
                dataset_chr = 0;
                for ii = 1:length(figure_details)
                    if (figure_details(ii).chr ~= 0)
                        if (strcmp(dataset_het{j1}.chr,chr_name{ii}) == 1)
                            dataset_chr = ii;
                        end;
                    end;
                end;
                if (parent_chr == dataset_chr) && (str2double(parent_het{i}.position) == str2double(dataset_het{j1}.position))
                    child = 'HET';
                    break;
                end;
                if (parent_chr == dataset_chr) && (str2double(parent_het{i}.position) < str2double(dataset_het{j1}.position))
                    j1 = j1-1;
                    if (j1 == 0);   j1 = 1;   end;
                    break;
                end;
                if (parent_chr < dataset_chr)
                    break;
                end;
            end;
            het_place = j1;
            % If child state has not been identified as 'HET', but is found in the hom SNP datasets, define child SNP as 'oddHET'.
            if (strcmp(child,'null') == 1)
                for j2 = hom_place:(HOM1_dataset_length+HOM2_dataset_length)
                    if (j2 <= HOM1_dataset_length)
                        dataset_chr = 0;
                        for ii = 1:length(figure_details)
                            if (figure_details(ii).chr ~= 0)
                                if (strcmp(dataset_hom1{j2}.chr,chr_name{ii}) == 1)
                                    dataset_chr = ii;
                                end;
                            end;
                        end;
                        if (parent_chr == dataset_chr) && (str2num(parent_het{i}.position) == str2num(dataset_hom1{j2}.position))
                            child = 'oddHET';
                            break;
                        end;
                        if (parent_chr == dataset_chr) && (str2num(parent_het{i}.position) < str2num(dataset_hom1{j2}.position))
                            j2 = j2-1;
                            if (j2 == 0);   j2 = 1;   end;
                            break;
                        end;
                        if (parent_chr < dataset_chr)
                            break;
                        end;
                    else
                        dataset_chr = 0;
                        for ii = 1:length(figure_details)
                            if (figure_details(ii).chr ~= 0)
                                if (strcmp(dataset_hom2{j2-HOM1_dataset_length}.chr,chr_name{ii}) == 1)
                                    dataset_chr = ii;
                                end;
                            end;
                        end;
                        if (parent_chr == dataset_chr) && (str2num(parent_het{i}.position) == str2num(dataset_hom2{j2-HOM1_dataset_length}.position))
                            child = 'oddHET';
                            break;
                        end;
                        if (parent_chr == dataset_chr) && (str2num(parent_het{i}.position) < str2num(dataset_hom2{j2-HOM1_dataset_length}.position))
                            j2 = j2-1;
                            if (j2 == 0);   j2 = 1;   end;
                            break;
                        end;
                        if (parent_chr < dataset_chr)
                            break;
                        end;
                    end;
                end;
                hom_place = j2;
            end;
            % If child SNP has not yet been identified, it is not found in child dataset and is defined as 'HOM'.
            if (strcmp(child,'null') == 1)
                child = 'HOM';
            end;
            
            if (strcmp(child,'HET') == 1)
                chr_SNPdata{parent_chr,1}(locus) = chr_SNPdata{parent_chr,1}(locus)+1;
            elseif (strcmp(child,'oddHET') == 1)
                chr_SNPdata{parent_chr,2}(locus) = chr_SNPdata{parent_chr,2}(locus)+1;
            elseif (strcmp(child,'HOM') == 1)
                chr_SNPdata{parent_chr,3}(locus) = chr_SNPdata{parent_chr,3}(locus)+1;
            else   % It was heterozygous and now is not, thus it is now homozygous.
                chr_SNPdata{parent_chr,3}(locus) = chr_SNPdata{parent_chr,3}(locus)+1;
            end;
        end;
        if (mod(i,10000) == 0)
            fprintf('.');
        end;
        if (mod(i,6000000) == 0)
            fprintf('\n\t\t');
        end;
    end;
    fprintf(']');
    if exist('h','var'); delete(h); clear h; end;
    save([workingDir 'matlab_dir/' fileName],'chr_SNPdata');
else
    % Load chr plot data if the save file does exist.
    load([workingDir 'matlab_dir/' fileName]);
end;

% basic plot parameters not defined per genome.
TickSize        = -0.005;  %negative for outside, percentage of longest chr figure.
maxY            = 10;

%define colors for colorBars plot
colorNoData = [1.0   1.0   1.0  ]; %used when no data is available for the bin.
colorInit   = [0.5   0.5   0.5  ]; %external; used in blending at ends of chr.
colorHET    = [0.0   0.0   0.0  ]; % near 1:1 ratio SNPs
colorOddHET = [0.0   1.0   0.0  ]; % Het, but not near 1:1 ratio SNPs.
colorHOM    = [1.0   0.0   0.0  ]; % Hom SNPs;

%% -----------------------------------------------------------------------------------------
% Setup for main figure generation.
%------------------------------------------------------------------------------------------
fig = figure(1);
set(gcf, 'Position', [0 70 1024 600]);
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
    Linear_maxY            = 10;   

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
    fprintf(['\tfigposition = [' num2str(left) ' | ' num2str(bottom) ' | ' num2str(width) ' | ' num2str(height) ']\n']);
    hold on;

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
                %c_post =   colorHET*HETplot2{chr}(i) + ...
                %           colorHOM*HOMplot2{chr}(i) + ...
                %           colorNoData*(1-min([HETplot2{chr}(i)+HOMplot2{chr}(i) 1]));
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

    % axes labels etc.
    hold off;
    xlim([0,chr_size(chr)*chr_length_scale_multiplier]);
    
    %% modify y axis limits to show annotation locations if any are provided.
    if (length(annotations) > 0)
        ylim([-1.5,maxY]);
    else
        ylim([0,maxY]);
    end;
    set(gca,'YTick',[]);
    set(gca,'TickLength',[(TickSize*chr_size(largestChr)/chr_size(chr)) 0]); %ensures same tick size on all subfigs.
    ylabel(chr_label{chr}, 'Rotation', 90, 'HorizontalAlign', 'center', 'VerticalAlign', 'bottom');
    set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
    set(gca,'XTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8','2.0','2.2','2.4','2.6','2.8','3.0','3.2'});
    set(gca,'YTick',[maxY/4 maxY/2 maxY/4*3 maxY]);
    set(gca,'YTickLabel',{'1','2','3','4'});
    set(gca,'FontSize',6);
    if (chr == find(chr_posY == max(chr_posY)))
        if (strcmp(parent_name,child_name) == 1)
            title([ parent_name ' SNP map'],'Interpreter','none','FontSize',12);
        else
            title([ parent_name ' vs. ' child_name ' SNP/LOH map'],'Interpreter','none','FontSize',12);
        end;
    end;
    hold on;
    %end axes labels etc.
    
    %show centromere outlines and horizontal marks.
    x1 = cen_start(chr)*chr_length_scale_multiplier;
    x2 = cen_end(chr)*chr_length_scale_multiplier;
    leftEnd  = 0.5*(5000/bases_per_bin);
    rightEnd = chr_size(chr)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
    if (Centromere_format == 0)
        % standard chromosome cartoons in a way which will not cause segfaults when running via commandline.
        dx     = 5*(5000/bases_per_bin);
        dy     = 2;
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
    end;
    %end show centromere.
    
    %show annotation locations
    if (show_annotations) && (length(annotations) > 0)
        plot([leftEnd rightEnd], [-1.5 -1.5],'color',[0 0 0]);
        hold on;
        annotation_location = (annotation_start+annotation_end)./2;
        for i = 1:length(annotation_location)
            if (annotation_chr(i) == chr)
                annotationloc = annotation_location(i)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
                plot(annotationloc,-1.5,'k:o','MarkerEdgeColor',annotation_edgecolor{i},'MarkerFaceColor',annotation_fillcolor{i},'MarkerSize',annotation_size(i));
                % plot(annotationloc,-1.5,'k:o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);
            end;
        end;
        hold off;
    end;
    %end show annotation locations.

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

%{
        %% cgh plot section.
        c_ = [0 0 0];
        fprintf(['chr' num2str(chr) ':' num2str(length(CNVplot2{chr})) '\n']);
        for i = 1:length(CNVplot2{chr}); 
            x_ = [i i i-1 i-1];
            % ploidy_factor : is used to adjust default ploidy in figure generation.  Haploid default ploidy position is set to 1/4 of the
            % scale of 1..4 while diploid position is set to 1/2 of the scale.
            if (ploidy_default == 1)
                ploidy_factor = 4;
            elseif (ploidy_default == 2)
                ploidy_factor = 2;
            end;
            if (ploidy == 0) % no value; no scaling.
                y_ = [maxY/ploidy_factor CNVplot2{chr}(i)*maxY/2 CNVplot2{chr}(i)*maxY/2 maxY/ploidy_factor];
            else
                if (Yscale_nearest_even_ploidy == true)
                    if (abs(2-ploidy) < abs(4-ploidy)) % nearer to diploid   
                        y_ = [maxY/ploidy_factor CNVplot2{chr}(i)*ploidy/2*maxY/2 CNVplot2{chr}(i)*ploidy/2*maxY/2 maxY/ploidy_factor];
                    else
                        y_ = [maxY/ploidy_factor CNVplot2{chr}(i)*ploidy/4*maxY/2 CNVplot2{chr}(i)*ploidy/4*maxY/2 maxY/ploidy_factor];
                    end;
                else
                    y_ = [maxY/ploidy_factor CNVplot2{chr}(i)*ploidy/2*maxY/2 CNVplot2{chr}(i)*ploidy/2*maxY/2 maxY/ploidy_factor];
                end; 
            end;
        
            % makes a blackbar for each bin.
            f = fill(x_,y_,c_);
            set(f,'linestyle','none');
        end;
        x2 = chr_size(chr)*chr_length_scale_multiplier;
        plot([0; x2], [maxY/2; maxY/2],'color',[0 0 0]);  % 2n line.
        % draw lines across plots for easier interpretation of CNV regions.
        if (ploidy_default == 2)
            line([0 x2], [maxY/4*3 maxY/4*3],'Color',[0.85 0.85 0.85]);
            line([0 x2], [maxY/4   maxY/4],  'Color',[0.85 0.85 0.85]);
        elseif (ploidy_default == 1)
            line([0 x2], [maxY/4*3 maxY/4*3],'Color',[0.85 0.85 0.85]);      
            line([0 x2], [maxY/2   maxY/2],  'Color',[0.85 0.85 0.85]);
        end;
        %% end cgh plot section.
%}
                    
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
        x1 = cen_start(chr)*chr_length_scale_multiplier;
        x2 = cen_end(chr)*chr_length_scale_multiplier;
        leftEnd  = 0.5*(5000/bases_per_bin);
        rightEnd = chr_size(chr)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);

        if (Centromere_format == 0)
            % standard chromosome cartoons in a way which will not cause segfaults when running via commandline.
            dx     = 5*(5000/bases_per_bin);
            dy     = 2;
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
        end;
        %end show centromere.

        %show annotation locations
        if (show_annotations) && (length(annotations) > 0)
            plot([leftEnd rightEnd], [-1.5 -1.5],'color',[0 0 0]);
            hold on;
            annotation_location = (annotation_start+annotation_end)./2;
            for i = 1:length(annotation_location)
                if (annotation_chr(i) == chr)
                    annotationloc = annotation_location(i)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
                    plot(annotationloc,-1.5,'k:o','MarkerEdgeColor',annotation_edgecolor{i},'MarkerFaceColor',annotation_fillcolor{i},'MarkerSize',annotation_size(i));
                    % plot(annotationloc,-1.5,'k:o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);
                end;
            end;
            hold off;
        end;
        %end show annotation locations.

        %% Final formatting stuff.
        xlim([0,chr_size(chr)*chr_length_scale_multiplier]);
        % modify y axis limits to show annotation locations if any are provided.
        if (length(annotations) > 0)
            ylim([-1.5,maxY]);
        else
            ylim([0,maxY]);
        end;
        set(gca,'TickLength',[(Linear_TickSize*chr_size(1)/chr_size(chr)) 0]); %ensures same tick size on all subfigs.
        set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
        set(gca,'XTickLabel',[]);
        if (chr == 1)
            if (strcmp(parent_name,child_name) == 1)
                ylabel(parent_name, 'Rotation', 0, 'HorizontalAlign', 'right', 'VerticalAlign', 'bottom','Interpreter','none','FontSize',5);
            else
                ylabel({parent_name;'vs.';child_name}, 'Rotation', 0, 'HorizontalAlign', 'right', 'VerticalAlign', 'bottom','Interpreter','none','FontSize',5);
            end;
            set(gca,'YTick',[Linear_maxY/4 Linear_maxY/2 Linear_maxY/4*3 Linear_maxY]);
            set(gca,'YTickLabel',{'1','2','3','4'});
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

%   % Main figure colors key.
%   subplot('Position',[0.65 0.2 0.2 0.4]);
%   axis off square;
%   xlim([-0.1,1]);
%   ylim([-0.1,1.6]);
%   set(gca,'XTick',[]);
%   set(gca,'YTick',[]);
%   patch([0 0.2 0.2 0], [1.4 1.4 1.5 1.5], colorNoData);   text(0.3,1.45,'Low SNP density');
%   patch([0 0.2 0.2 0], [1.2 1.2 1.3 1.3], colorHET);      text(0.3,1.25,'Heterozygous SNP density');
%   patch([0 0.2 0.2 0], [1.0 1.0 1.1 1.1], colorOddHET);   text(0.3,1.05,'non-1:1 ratio Heterozygous SNP density');
%   patch([0 0.2 0.2 0], [0.8 0.8 0.9 0.9], colorHOM);      text(0.3,0.85,'Homozygous SNP density');
if (strcmp(parent_name,child_name) == 1)
    saveas(fig,        [figureDir parent_name '.SNP-map.1.eps'], 'epsc');
    set(Linear_fig,'PaperPosition',[0 0 8 0.62222222]);
    saveas(Linear_fig, [figureDir parent_name '.SNP-map.2.eps'], 'epsc');
else
    saveas(fig,        [figureDir parent_name '->' child_name '.SNP-map.1.eps'], 'epsc');
    set(Linear_fig,'PaperPosition',[0 0 8 0.62222222]);
    saveas(Linear_fig, [figureDir parent_name '->' child_name '.SNP-map.2.eps'], 'epsc');
end;
delete(fig);
delete(Linear_fig);

%% ========================================================================
% end stuff
%==========================================================================
end
