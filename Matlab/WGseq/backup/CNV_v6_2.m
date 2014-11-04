function [] = CNV_v6_2(projectName,genome,ploidyString,ploidyBaseString, CNV_verString,rDNA_verString,workingDir,figureDir,displayBREAKS,  referenceCHR)
%% ========================================================================
Centromere_format_default   = 0;
Yscale_nearest_even_ploidy  = true;
HistPlot                    = true;
ChrNum                      = true;
show_annotations            = true;
analyze_rDNA                = true;
Linear_display              = true;
Low_quality_ploidy_estimate = true;

% determine which GBrowse output format to use.
if (strcmp(genome,'Ca_a') == 1) || (strcmp(genome,'Ca_g') == 1)
    Output_CNVs_as_CGD_annotations = true;
    Output_CNVs_as_SGD_annotations = false;
elseif (strcmp(genome,'Sa_c') == 1)
    % Saccharomyces cerevisiae genome fasta file had chromosomes renamed from what was recieved from SGD.
    % This was done so fasta chr names would be a single string which could be used for matching preprocessed data with chrs in figure generation.
    % The chromosome names in SGD Genome Browser are yet again something else, so output to SGD annotations has to be handled a little differently than CGD.
    Output_CNVs_as_CGD_annotations = false;
    Output_CNVs_as_SGD_annotations = true;
else
    Output_CNVs_as_CGD_annotations = false;
    Output_CNVs_as_SGD_annotations = false;
end;
%%=========================================================================
% Control variables.
%--------------------------------------------------------------------------
% Defines chr sizes in bp. (diploid total=28,567,7888)
% Defines centromere locations in bp.
% Defines annotation locations in bp.
fprintf(['\n$$ 1 : ' workingDir '\n$$ 2 : ' figureDir '\n$$ 3 : ' genome '\n']);
[centromeres, chr_sizes, figure_details, annotations, ploidy_default] = Load_genome_information_1(workingDir,figureDir,genome);
[Aneuploidy]                                                          = Load_dataset_information_1(projectName,workingDir);

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
        fprintf(['\t[' num2str(annotations(i).chr) ':' annotations(i).type ':' num2str(annotations(i).start) ':' num2str(annotations(i).end) ':' annotations(i).fillcolor ':' ...
            annotations(i).edgecolor ':' num2str(annotations(i).size) ']\n']);
    end;
end;

for i = 1:length(figure_details)
    if (figure_details(i).chr == 0)
	if (strcmp(figure_details(i).label,'Key') == 1)
	    key_posX   = figure_details(i).posX;
	    key_posY   = figure_details(i).posY;
	    key_width  = figure_details(i).width;
	    key_height = figure_details(i).height;
	end;
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
% bases_per_bin			= 5000;
bases_per_bin			= max(chr_size)/700;
chr_length_scale_multiplier	= 1/bases_per_bin;
CGD_bases_per_bin		= 1000;
CGD_chr_length_scale_multiplier	= 1/CGD_bases_per_bin;

% Defines rDNA region.
if (analyze_rDNA == true)
    if (strcmp(genome,'Ca_a') == 1)
        chrRDNA       = 8;
        positionRDNA1 = 1884387;
        positionRDNA2 = 1897582;
        sizeRDNA      = positionRDNA2-positionRDNA1+1;
        sizeREF       = chr_size(referenceCHR);
    else
        fprintf(['\nWARNING:[analyze_CNVs.m]: rDNA copy number estimate not yet implemented for genome: ''' genome '''\n']);
	analyze_rDNA = false;
    end;
end;

%%=========================================================================
%%= No further control variables below. ===================================
%%=========================================================================

% Sanitize user input of euploid state.
ploidyBase = round(str2num(ploidyBaseString));
if (ploidyBase > 4);   ploidyBase = 4;   end;
if (ploidyBase < 1);   ploidyBase = 1;   end; 
fprintf(['\nEuploid base = "' num2str(ploidyBase) '"\n']);

fprintf(['\nGenerating CNV figure from ''' projectName ''' sequence data.\n']);

% Initializes vectors used to hold copy number data.
for chr = 1:num_chr   % number of chrs.
    for j = 1:2       % 2 categories tracked : total read counts in bin and size of bin in base-pairs.
        chr_CNVdata{chr,j}     = zeros(1,ceil(chr_size(chr)/bases_per_bin));
        CGD_chr_CNVdata{chr,j} = zeros(1,ceil(chr_size(chr)/CGD_bases_per_bin));
    end;
end;
% Initializes vectors used to hold other site data, for normalization purposes.
for j = 1:2
    nonchr_CNVdata{1,j}        = zeros(1,ceil(chr_size(chr)/bases_per_bin));
end;

if (exist([workingDir 'matlab_dir/' projectName '.CNV_' CNV_verString '.ploidy_' ploidyString '.mat'],'file') == 0) || ...
	(exist([workingDir 'matlab_dir/' projectName '.rDNA-CNV_' rDNA_verString '.ploidy_' ploidyString '.mat'],'file') == 0)
    fprintf('\nMAT file not found, regenerating.\n');
    datafile = [workingDir 'pileup_dir/' projectName '_putative_CNVs_' CNV_verString '.txt'];
    data     = fopen(datafile);
    k = 0;
    j = 0;
    
    if (analyze_rDNA == 1)
        countRDNA = 0;
        countREF  = 0;
    end;
    for i = 1:length(chr_name)
        fprintf(['\nchr' num2str(i) ' = ''' chr_name{i} '''.\tCGHlength = ' num2str(length(chr_CNVdata{i,1}))]);
    end;
    fprintf('\n');
    error_counter_1 = 0;
    error_counter_2 = 0;
    old_chr = -1;
    while ~feof(data)
        j = j+1;
        
        chrID    = fscanf(data,'%s',1);	% chr ID of bp.
        position = fscanf(data,'%s',1);	% chr position of bp.
        count    = fscanf(data,'%s',1);	% read count at bp.
        quality  = fscanf(data,'%s',1);	% quality of bp, not used, but needed to ensure proper data loading.
        
        % Identify which chromosome this data point corresponds to.
        chr = 0;
        for i = 1:length(chr_name)
            if (strcmp(chrID,chr_name(i)) == 1)
                chr = chr_id(i);
            end;
        end;
        
        % tabulate rDNA copies compared to reference chromosome.
        if (analyze_rDNA == 1)
            if (chr == 8)
                if (str2num(position) >= positionRDNA1) && (str2num(position) <= positionRDNA2)
                    countRDNA = countRDNA+str2num(count);
                end;
            end;
            if (chr == referenceCHR)
                countREF = countREF+str2num(count);
            end;
        end;
        
        % val is the genomic location bin of the bp.
        val     = ceil(str2double(position)/bases_per_bin);
        CGD_val	= ceil(str2double(position)/CGD_bases_per_bin);
        % count of data points at locus.
        
        if (chr > 0)
	    % Current data point belongs to a used chromosome.
            if (val <= length(chr_CNVdata{chr,1}))
                chr_CNVdata{chr,1}(val) = chr_CNVdata{chr,1}(val)+str2num(count);
                chr_CNVdata{chr,2}(val) = chr_CNVdata{chr,2}(val)+1;
	    else
		error_counter_2 = error_counter_2+1;
            end;
            if (CGD_val <= length(CGD_chr_CNVdata{chr,1}))
                CGD_chr_CNVdata{chr,1}(CGD_val) = CGD_chr_CNVdata{chr,1}(CGD_val)+str2num(count);
                CGD_chr_CNVdata{chr,2}(CGD_val) = CGD_chr_CNVdata{chr,2}(CGD_val)+1;
            end;
        elseif (chr == 0)
	    % Current data point does not belong to a used chromosome.
            error_counter_1 = error_counter_1+1;
        end;

	% Log file outputs of status.
	if (old_chr ~= chr)
	    if (chr == 0)
		fprintf(['\n\tData not used in figure.\n\t\t']);
	    else
		fprintf(['\n\tchr' chrID '\n\t\t']);
	    end;
	end;
	old_chr = chr;
	if (mod(j,10000) == 0)
	    fprintf('.');
	end;
	if (mod(j,10000*80) == 0)
	    fprintf('\n\t\t');
	end;
    end;
    fprintf(['\n# lines not matching chrs (mito, etc) : ' num2str(error_counter_1) ]);
    fprintf(['\n# lines with positions outside chr    : ' num2str(error_counter_2) ]);
    
    save([workingDir 'matlab_dir/' projectName '.CNV_' CNV_verString '.ploidy_' ploidyString '.mat'],'chr_CNVdata','CGD_chr_CNVdata');
    if (analyze_rDNA == 1)
        save([workingDir 'matlab_dir/' projectName '.rDNA-CNV_' rDNA_verString '.ploidy_' ploidyString '.mat'],'countRDNA','sizeRDNA','countREF','sizeREF');
    end;
else
    fprintf('\nMAT file found, loading.\n');
    load([workingDir 'matlab_dir/' projectName '.CNV_' CNV_verString '.ploidy_' ploidyString '.mat']);
end;

% basic plot parameters not defined per genome.
TickSize        = -0.005;  %negative for outside, percentage of longest chr figure.
maxY            = ploidyBase*2;

%% -----------------------------------------------------------------------------------------
% Setup for CGD output.
%-------------------------------------------------------------------------------------------
if (Output_CNVs_as_CGD_annotations == true) || (Output_CNVs_as_SGD_annotations == true)
    for chr = 1:num_chr
        CGD_CNVplot{chr} = CGD_chr_CNVdata{chr};
        CGD_chr_max(chr) = max(CGD_CNVplot{chr});
        CGD_chr_med(chr) = median(CGD_CNVplot{chr});
    end;
    for chr = 1:num_chr
        CGD_max_count     = max(CGD_chr_max);
        CGD_median_count  = sum(CGD_chr_med)/length(CGD_chr_med);
        CGD_CNVplot2{chr} = CGD_CNVplot{chr}/CGD_median_count;
    end;
    ploidy = str2num(ploidyString);
    fprintf(['Ploidy string = "' ploidyString '"\n']);

    [chr_breaks, chrCopyNum, ploidyAdjust] = FindChrSizes_2(Aneuploidy,CGD_CNVplot2,ploidy,num_chr);
    
    CGDid = fopen([figureDir projectName '.CNV-map.fff'],'w');
    fprintf(CGDid,'[CNV_estimates]\n');
    fprintf(CGDid,'glyph=xyplot\n');
    fprintf(CGDid,'graph_type=boxes\n');
    fprintf(CGDid,'scale=left\n');
    fprintf(CGDid,'fgcolor=black\n');
    fprintf(CGDid,'bgcolor=black\n');
    fprintf(CGDid,'height=40\n');
    if (ploidy_default == 1)
        fprintf(CGDid,'min_score=-1\n');
        fprintf(CGDid,'max_score=3\n');
        CGD_maxY = 4;
        CGD_min  = -1;
        CGD_max  = 3;
    else % (ploidy_default == 1)
        if (abs(2-ploidy) < abs(4-ploidy)) % nearer to diploid
            fprintf(CGDid,'min_score=-2\n');
            fprintf(CGDid,'max_score=2\n');
            CGD_maxY = 4;
            CGD_min  = -2;
            CGD_max  = 2;
        else
            fprintf(CGDid,'min_score=-4\n');
            fprintf(CGDid,'max_score=4\n');
            CGD_maxY = 8;
            CGD_min  = -4;
            CGD_max  = 4;
        end;
    end;
    fprintf(CGDid,'label=1\n');
    fprintf(CGDid,'key=CNV estimate\n');
    fprintf(CGDid,'\n');
    
    for chr = 1:num_chr
        if (Output_CNVs_as_CGD_annotations == true)
            annotation_chr_string = chr_name{chr};
        elseif (Output_CNVs_as_SGD_annotations == true)
            annotation_chr_string = ['Chr' chr_label{chr}];
        end;
        fprintf(CGDid,['reference="' annotation_chr_string '"\n']);
        
        for i = 1:length(CGD_CNVplot2{chr});
            if (ploidy == 0) % no value; no scaling.
                fprintf(CGDid,['CNV_estimates "Differences from expected" ' annotation_chr_string ':' num2str(1+(i-1)*CGD_bases_per_bin) '..' num2str(i*CGD_bases_per_bin) ' score=' num2str(CGD_CNVplot2{chr}(i)*CGD_maxY/2 + CGD_min) '\n']);
                % fprintf(CGDid,['CNV_estimates\t' annotation_chr_string '\t' num2str(1+(i-1)*CGD_bases_per_bin) '..' num2str(i*CGD_bases_per_bin) '\tscore=' num2str(CGD_CNVplot2{chr}(i)*CGD_maxY/2 + CGD_min) '\n']);
                % fprintf(CGDid,[annotation_chr_string '\t' num2str((i-1)*CGD_bases_per_bin) '\t' num2str(i*CGD_bases_per_bin-1) '\tscore=' num2str(CGD_CNVplot2{chr}(i)*CGD_maxY/2 + CGD_min) '\n']);
                % y_ = [maxY/ploidy_factor CNVplot2{chr}(i)*maxY/2 CNVplot2{chr}(i)*maxY/2 maxY/ploidy_factor];
            else
                if (Yscale_nearest_even_ploidy == true)
                    if (abs(2-ploidy) < abs(4-ploidy)) % nearer to diploid
                        fprintf(CGDid,['CNV_estimates "Differences from expected" ' annotation_chr_string ':' num2str(1+(i-1)*CGD_bases_per_bin) '..' num2str(i*CGD_bases_per_bin) ' score=' num2str(CGD_CNVplot2{chr}(i)*ploidy/2*CGD_maxY/2 + CGD_min) '\n']);
                        % fprintf(CGDid,['CNV_estimates\t' annotation_chr_string '\t' num2str(1+(i-1)*CGD_bases_per_bin) '..' num2str(i*CGD_bases_per_bin) '\tscore=' num2str(CGD_CNVplot2{chr}(i)*ploidy/2*CGD_maxY/2 + CGD_min) '\n']);
                        % fprintf(CGDid,[annotation_chr_string '\t' num2str((i-1)*CGD_bases_per_bin) '\t' num2str(i*CGD_bases_per_bin-1) '\tscore=' num2str(CGD_CNVplot2{chr}(i)*ploidy/2*CGD_maxY/2 + CGD_min) '\n']);
                        % y_ = [maxY/ploidy_factor CNVplot2{chr}(i)*ploidy/2*maxY/2 CNVplot2{chr}(i)*ploidy/2*maxY/2 maxY/ploidy_factor];
                    else
                        fprintf(CGDid,['CNV_estimates "Differences from expected" ' annotation_chr_string ':' num2str(1+(i-1)*CGD_bases_per_bin) '..' num2str(i*CGD_bases_per_bin) ' score=' num2str(CGD_CNVplot2{chr}(i)*ploidy/4*CGD_maxY/2 + CGD_min) '\n']);
                        % fprintf(CGDid,['CNV_estimates\t' annotation_chr_string '\t' num2str(1+(i-1)*CGD_bases_per_bin) '..' num2str(i*CGD_bases_per_bin) '\tscore=' num2str(CGD_CNVplot2{chr}(i)*ploidy/4*CGD_maxY/2 + CGD_min) '\n']);
                        % fprintf(CGDid,[annotation_chr_string '\t' num2str((i-1)*CGD_bases_per_bin) '\t' num2str(i*CGD_bases_per_bin-1) '\tscore=' num2str(CGD_CNVplot2{chr}(i)*ploidy/4*CGD_maxY/2 + CGD_min) '\n']);
                        % y_ = [maxY/ploidy_factor CNVplot2{chr}(i)*ploidy/4*maxY/2 CNVplot2{chr}(i)*ploidy/4*maxY/2 maxY/ploidy_factor];
                    end;
                else
                    fprintf(CGDid,['CNV_estimates "Differences from expected" ' annotation_chr_string ':' num2str(1+(i-1)*CGD_bases_per_bin) '..' num2str(i*CGD_bases_per_bin) ' score=' num2str(CGD_CNVplot2{chr}(i)*ploidy/2*CGD_maxY/2 + CGD_min) '\n']);
                    % fprintf(CGDid,['CNV_estimates\t' annotation_chr_string '\t' num2str(1+(i-1)*CGD_bases_per_bin) '..' num2str(i*CGD_bases_per_bin) '\tscore=' num2str(CGD_CNVplot2{chr}(i)*ploidy/2*CGD_maxY/2 + CGD_min) '\n']);
                    % fprintf(CGDid,[annotation_chr_string '\t' num2str((i-1)*CGD_bases_per_bin) '\t' num2str(i*CGD_bases_per_bin-1) '\tscore=' num2str(CGD_CNVplot2{chr}(i)*ploidy/2*CGD_maxY/2 + CGD_min) '\n']);
                    % y_ = [maxY/ploidy_factor CNVplot2{chr}(i)*ploidy/2*maxY/2 CNVplot2{chr}(i)*ploidy/2*maxY/2 maxY/ploidy_factor];
                end;
            end;
        end;
        fprintf(CGDid,'\n');
    end;
    fclose(CGDid);
end;

%% -----------------------------------------------------------------------------------------
% Setup for main figure generation.
%-------------------------------------------------------------------------------------------
fig = figure(1);
set(gcf, 'Position', [0 70 1024 600]);
for chr = 1:num_chr
    for pos = 1:length(chr_CNVdata{chr,1})
        % Add to the plotting array : the sum of the data in each region, divided by the number of data points in each region.
        if (chr_CNVdata{chr,2}(pos) == 0)
            % If there are no data elements in the region, assign a null value.
            CNVplot{chr}(pos) = 0;
        else
            % if there are data elements in the region, divide the sum by the number of data elements.
            CNVplot{chr}(pos) = chr_CNVdata{chr,1}(pos)/chr_CNVdata{chr,2}(pos);
        end;
    end;
    % CNVplot{chr} = chr_CNVdata{chr};

    chr_max(chr) = max(CNVplot{chr});
    chr_med(chr) = median(CNVplot{chr});
end;
max_count     = max(chr_max);
median_count  = sum(chr_med)/length(chr_med);
for chr = 1:num_chr
    CNVplot2{chr} = CNVplot{chr}/median_count;
end;
ploidy = str2num(ploidyString);
[chr_breaks, chrCopyNum, ploidyAdjust] = FindChrSizes_2(Aneuploidy,CNVplot2,ploidy,num_chr);
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

%% Initialize copy numbers string.
stringChrCNVs = '';

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
    
    %% cgh plot section.
    c_ = [0 0 0];
    fprintf(['chr' num2str(chr) ':' num2str(length(CNVplot2{chr})) '\n']);
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

    % modify y axis limits to show annotation locations if any are provided.
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
        title([ projectName ' CNV map'],'Interpreter','none','FontSize',12);
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

	stringChrCNVs = [stringChrCNVs ';' chr_string];
    end;


    %% Linear figure draw section
    if (Linear_display == true)
        figure(Linear_fig);
        Linear_width = Linear_Chr_max_width*chr_size(chr)/Linear_genome_size;
        subplot('Position',[Linear_left Linear_base Linear_width Linear_height]);
        Linear_left = Linear_left + Linear_width + Linear_left_chr_gap;
        hold on;
        title(chr_label{chr},'Interpreter','none','FontSize',10);

        %% cgh plot section.
        c_ = [0 0 0];
        fprintf(['chr' num2str(chr) ':' num2str(length(CNVplot2{chr})) '\n']);
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
            ylabel(projectName, 'Rotation', 0, 'HorizontalAlign', 'right', 'VerticalAlign', 'bottom','Interpreter','none','FontSize',5);

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

% Save primary genome figure.
saveas(fig,        [figureDir projectName '.CNV-map.1.eps'], 'epsc');
delete(fig);

% Save horizontal aligned genome figure.
set(Linear_fig,'PaperPosition',[0 0 8 0.62222222]);
saveas(Linear_fig, [figureDir projectName '.CNV-map.2.eps'], 'epsc');
delete(Linear_fig);

% Output chromosome copy number estimates.
textFileName = [figureDir projectName '.CNV-map.3.txt'];
fprintf(['Text output of CNVs : "' textFileName '"\n']);
textFileID = fopen(textFileName,'w');
fprintf(textFileID,stringChrCNVs);
fclose(textFileID);

end
