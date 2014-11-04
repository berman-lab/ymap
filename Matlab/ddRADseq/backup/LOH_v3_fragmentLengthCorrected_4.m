function [] = LOH_v3_fragmentLengthCorrected_3(parent_name,child_name,genome,ave_copy_num, ploidyEstimateString,ploidyBaseString, ...
						SNP_verString,LOH_verString,workingDir,figureDir,displayBREAKS)
%% ========================================================================
%    Centromere_format          : Controls how centromeres are depicted.   [0..2]   '2' is pinched cartoon default.
%    scale_type                 : 'Ratio' or 'Log2Ratio' y-axis scaling of copy number.
%                                 'Log2Ratio' does not properly scale CGH data by ploidy.
%    Chr_max_width              : max width of chrs as fraction of figure width.

%
% Restriction fragments are called as heterozygous if they contain at least one good heterozygous locus, rather than an average value across length of fragment.
%

Centromere_format_default   = 0;
Chr_max_width               = 0.8;
colorBars                   = true;
show_annotations            = true;
Yscale_nearest_even_ploidy  = true;
Linear_display              = true;
Low_quality_ploidy_estimate = false;
HistPlot                    = false;
ChrNum                      = false;

%%=========================================================================
% Control variables for Candida albicans SC5314.
%--------------------------------------------------------------------------
% Defines chr sizes in bp. (diploid total=28,567,7888)
% Defines centromere locations in bp.
% Defines MRS locations in bp.
[centromeres, chr_sizes, figure_details, annotations, ploidy_default] = Load_genome_information_1(workingDir,figureDir,genome);
[Aneuploidy]                                                          = Load_dataset_information_1(parent_name,workingDir);

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
% bases_per_bin                 = 5000;
bases_per_bin                   = max(chr_size)/700;
chr_length_scale_multiplier     = 1/bases_per_bin;

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

%%=========================================================================
%%= No further control variables below. ===================================
%%=========================================================================
if (strcmp(parent_name,child_name) == 1)
    fprintf(['\nGenerating SNP-map figure from ' parent_name ' genome data.\n']);
else
    fprintf(['\nGenerating LOH-map figure from ' parent_name ' and ' child_name ' genome data.\n']);
end;

%% Determine reference genome FASTA file in use.
%  Read in and parse : "links_dir/main_script_dir/genome_specific/[genome]/reference.txt"
reference_file   = [workingDir 'main_script_dir/genomeSpecific/' genome '/reference.txt'];
refernce_fid     = fopen(reference_file, 'r');
refFASTA         = fgetl(refernce_fid);
fclose(refernce_fid);

%% load pre-processed CNV data : project
if (exist([workingDir 'matlab_dir/' child_name '.fragment_CNV_data.mat'],'file') == 0)
    fprintf('## Process project for CNV first.');
    error('## Process project for CNV first.');
else
    load([workingDir 'matlab_dir/' child_name '.fragment_CNV_data.mat']);
end;
project_fragments_CNV = fragments_CNV;

if (strcmp(parent_name,child_name) == 0)
    %% load pre-processed CNV data : reference
    if (exist([workingDir 'matlab_dir/' parent_name '.fragment_CNV_data.mat'],'file') == 0)
	fprintf('## Process reference for CNV first.');
	error('## Process reference for CNV first.');
    else
	load([workingDir 'matlab_dir/' parent_name '.fragment_CNV_data.mat']);
    end;
    reference_fragments_CNV = fragments_CNV;
else
    reference_fragments_CNV = project_fragments_CNV;
end;

numFragments = length(fragments_CNV);

%% Standardize 'fragment_data' data structure to hold results.
fragment_data = project_fragments_CNV;
numFragments  = length(fragment_data);
for frag = 1:numFragments
    fragment_data(frag).ref_read_count = reference_fragments_CNV(frag).read_count;  
    fragment_data(frag).ref_read_max   = reference_fragments_CNV(frag).read_max;
    fragment_data(frag).ref_read_ave   = reference_fragments_CNV(frag).read_ave;
%    fragment_data(frag).repet_count    = fragments_repet(frag).repet_count;
%    fragment_data(frag).repet_max      = fragments_repet(frag).repet_max;
%    fragment_data(frag).repet_ave      = fragments_repet(frag).repet_ave;
end;


%%================================================================================================
% Load SNP/LOH data.
%-------------------------------------------------------------------------------------------------

%% Load pre-processed ddRADseq fragment SNP data for project.
if (exist([workingDir 'matlab_dir/' child_name '.fragment_SNP_data_raw.mat'],'file') == 0)
    % ###	Coordinate annotation.  "[ ]" not usable fragment.   "[@]" odd usable fragment.   "[*]" even usable fragment.
    % ###		[chrID]:[bp]	A	T	G	C	Allelic fraction : calculated as 1-max(max(A,T,G,C)/(A+T+G+C),0.5)
    % ###											[copyEstimate]:[fragmentID]
    % ##	[@]	chr1:14021	1	51	0	0	0.01923076923077	3:[1]
    % ##	[@]	chr1:14248	0	2	45	0	0.04255319148936	3:[1]
    % ### chr	bp_Start	bp_end	copyNum	binsData
    % 1       14001   14291   2       20:0:0
    % 1       15900   16123   3       18:0:0
    % 1       21632   21834   2       14:0:0
    fprintf('Loading project results from Python script, which pre-processed the putative_SNPs relative to genome restriction fragments and corrected fragment copy number estimates.\n');
    datafile_RADseq  = [workingDir 'pileup_dir/' child_name '_RADseq_digest_analysis_SNP.txt'];
    data_RADseq      = fopen(datafile_RADseq);
    count            = 0;
    fragments_SNP    = [];
    while ~feof(data_RADseq)
        % Load fragment data from pre-processed text file, single line.
        tline = fgetl(data_RADseq);
        
        % check if line is a comment.
        col_1 = sscanf(tline, '%s',1);
        if (strcmp(col_1,'##') == 1)	% If the first space-delimited string is '##', then the line contains per-bp information.
	    % [*] or [@] for useful coordinates.   [_] for non-useful coordinates.
            new_string = sscanf(tline, '%s',  2 );
            for i = 1:size(sscanf(tline,'%s', 1 ),2);   new_string(1) = [];   end;
            col_2 = new_string;
	    if (strcmp(col_2,'[*]') == 1) || (strcmp(col_2,'[@]') == 1)
		% The number of valid lines found so far...  the number of coordinates with putatitve SNP data on valid restriction fragments.
		count = count + 1;

		% chromosome and coordinate of bp.
		new_string = sscanf(tline, '%s',  3 );
		for i = 1:size(sscanf(tline,'%s', 2 ),2);   new_string(1) = [];   end;
		col_3 = new_string;
		coordinate = strsplit(col_3,':');
		chr = str2num(strrep(coordinate{1},'chr',''));
		pos = str2num(coordinate{2});
        
		% allelic ratio at coordinate.
		new_string = sscanf(tline, '%s',  8 );
		for i = 1:size(sscanf(tline,'%s', 7 ),2);   new_string(1) = [];   end;
		col_8 = new_string;
		ratio = str2num(col_8);

		% Identity of fragment.
		new_string = sscanf(tline, '%s',  9 );
		for i = 1:size(sscanf(tline,'%s', 8 ),2);   new_string(1) = [];   end;
		col_9 = new_string;
		ID_temp_1      = strsplit(col_9,':');
		ID_temp_2      = ID_temp_1{2};
		ID_temp_3      = ID_temp_2(2:(end-1));
		% ID_temp_3(1)   = [];
		% ID_temp_3(end) = [];
		fragmentID = str2num(ID_temp_3);

		% Add fragment data to data structure.
		data_SNP(count).chr        = chr;
		data_SNP(count).pos        = pos;
		data_SNP(count).ratio      = ratio;
		data_SNP(count).fragmentID = fragmentID;
	    end;
        end;
    end;
    fclose(data_RADseq);
    save([workingDir 'matlab_dir/' child_name '.fragment_SNP_data_raw.mat'], 'data_SNP');
else
    load([workingDir 'matlab_dir/' child_name '.fragment_SNP_data_raw.mat']);
end;
child_data_SNP = data_SNP;
clear data_SNP;

if (strcmp(parent_name,child_name) == 0)
    %% Load pre-processed ddRADseq fragment SNP data for project.
    if (exist([workingDir 'matlab_dir/' parent_name '.fragment_SNP_data_raw.mat'],'file') == 0)
        % ###       Coordinate annotation.  "[ ]" not usable fragment.   "[@]" odd usable fragment.   "[*]" even usable fragment.
        % ###               [chrID]:[bp]    A       T       G       C       Allelic fraction : calculated as 1-max(max(A,T,G,C)/(A+T+G+C),0.5)
        % ###                                                                                       [copyEstimate]:[fragmentID]
        % ##        [@]     chr1:14021      1       51      0       0       0.01923076923077        3:[1]
        % ##        [@]     chr1:14248      0       2       45      0       0.04255319148936        3:[1]
        % ### chr   bp_Start        bp_end  copyNum binsData
        % 1       14001   14291   2       20:0:0
        % 1       15900   16123   3       18:0:0
        % 1       21632   21834   2       14:0:0
        fprintf('Loading project results from Python script, which pre-processed the putative_SNPs relative to genome restriction fragments and corrected fragment copy number estimates.\n');
        datafile_RADseq  = [workingDir 'pileup_dir/' parent_name '_RADseq_digest_analysis_SNP.txt'];
        data_RADseq      = fopen(datafile_RADseq);
        count            = 0;
        fragments_SNP    = [];
        while ~feof(data_RADseq)
            % Load fragment data from pre-processed text file, single line.
            tline = fgetl(data_RADseq);
    
            % check if line is a comment.
            col_1 = sscanf(tline, '%s',1);
            if (strcmp(col_1,'##') == 1)    % If the first space-delimited string is '##', then the line contains per-bp information.
                % [*] or [@] for useful coordinates.   [_] for non-useful coordinates.
                new_string = sscanf(tline, '%s',  2 );
                for i = 1:size(sscanf(tline,'%s', 1 ),2);   new_string(1) = [];   end;
                col_2 = new_string;
                if (strcmp(col_2,'[*]') == 1) || (strcmp(col_2,'[@]') == 1)
                    % The number of valid lines found so far...  the number of coordinates with putatitve SNP data on valid restriction fragments.
                    count = count + 1;

                    % chromosome and coordinate of bp.
                    new_string = sscanf(tline, '%s',  2 );
                    for i = 1:size(sscanf(tline,'%s', 1 ),2);   new_string(1) = [];   end;
                    col_3 = new_string;
                    coordinate = strsplit(col_3,':')
                    chr = str2num(strrep(coordinate{1},'chr',''));
                    pos = str2num(coordinate{2});

		    % allelic ratio at coordinate.
                    new_string = sscanf(tline, '%s',  8 );
                    for i = 1:size(sscanf(tline,'%s', 7 ),2);   new_string(1) = [];   end;
                    col_8 = new_string;
                    ratio = str2num(col_8);

                    % Identity of fragment.
                    new_string = sscanf(tline, '%s',  9 );
                    for i = 1:size(sscanf(tline,'%s', 8 ),2);   new_string(1) = [];   end;
                    col_9 = new_string;
                    ID_temp_1 = strsplit(col_9,':');
                    ID_temp_2 = ID_temp_1{2};
                    ID_temp_2(1) = [];
                    ID_temp_2(end) =[];
                    ID = str2num(ID_temp_2);
                    fragmentID = str2num(ID);   

                    % Add fragment data to data structure.
                    data_SNP(count).chr        = chr;
                    data_SNP(count).pos        = pos;  
                    data_SNP(count).ratio      = ratio;
                    data_SNP(count).fragmentID = fragmentID;
                end;
            end;
        end;
        fclose(data_RADseq);
        save([workingDir 'matlab_dir/' parent_name '.fragment_SNP_data_raw.mat'], 'data_SNP');
    else
        load([workingDir 'matlab_dir/' parent_name '.fragment_SNP_data_raw.mat']);
    end;
    parent_data_SNP = data_SNP;
    clear data_SNP;
end;


%
%%
%%%
%%%%
%%%%
%%%%
%%%
%%
%

%%=========================================================================
%%= Analyze corrected SNP data for RADseq datasets. =======================
%%=========================================================================

numRatios = length(child_data_SNP);
fprintf(['Number of all ratios       = ' num2str(numRatios)    '\n']);
fprintf(['Number of all fragments    = ' num2str(numFragments) '\n']);
numUsableFragments = 0;
for frag = 1:numFragments;
    if (fragment_data(frag).usable == 1)
	numUsableFragments = numUsableFragments + 1;
    end;
end;
fprintf(['Number of usable fragments = ' num2str(numUsableFragments) '\n']);

% Initializes vectors used to hold SNP ratio data.
for chr = 1:num_chr   % number of chrs.
    chr_SNPdata_RADseq{chr} = zeros(1,chr_size(chr));
end;
%% Gather SNP data for map display.
if (exist([workingDir 'matlab_dir/' child_name '.raw_SNP_' SNP_verString '.mat'],'file') == 0)
    fprintf('\nMAT file containing project SNP information not found, regenerating from prior data files.\n');
    fprintf('\n# Loading SNP ratio values.');
    for ratioID = 1:numRatios
        chr    = child_data_SNP(ratioID).chr;
        pos    = child_data_SNP(ratioID).pos;
        ratio  = child_data_SNP(ratioID).ratio;
        fragID = child_data_SNP(ratioID).fragmentID;
	usable = fragment_data(fragID).usable;
        if (usable == 1)
            chr_SNPdata_RADseq{chr}(pos) = ratio;
        end;
    end;
    fprintf('\n# Fragment corrected SNP ratio values have been added to map bins.');
                
    save([workingDir 'matlab_dir/' child_name '.raw_SNP_' SNP_verString '.mat'], 'chr_SNPdata_RADseq');
else
    fprintf('\nMAT file containing project SNP information found, loading.\n');
    load([workingDir 'matlab_dir/' child_name '.raw_SNP_' SNP_verString '.mat']);
end;
chr_SNPdata_RADseq_child = chr_SNPdata_RADseq;
clear chr_SNPdata_RADseq;

if (strcmp(parent_name,child_name) == 0)
    % Initializes vectors used to hold SNP ratio data.
    for chr = 1:num_chr   % number of chrs.
        chr_SNPdata_RADseq{chr} = zeros(1,chr_size(chr)); 
    end;
    %% Gather SNP data for map display.
    if (exist([workingDir 'matlab_dir/' child_name '.raw_SNP_' SNP_verString '.mat'],'file') == 0)
        fprintf('\nMAT file containing project SNP information not found, regenerating from prior data files.\n');
        fprintf('\n# Loading SNP ratio values.');
        for ratioID = 1:numRatios
            if (fragment_data(frag).usable == 1)
                chr    = child_data_SNP(count).chr;
                pos    = child_data_SNP(count).pos;
                ratio  = child_data_SNP(count).ratio;
                fragID = child_data_SNP(count).fragmentID;       
                if (fragment_data(fragID).usable == 1)
                    chr_SNPdata_RADseq{chr}(pos) = ratio;
                end;
            end;
        end;
        fprintf('\n# Fragment corrected SNP ratio values have been added to map bins.');

        save([workingDir 'matlab_dir/' child_name '.raw_SNP_' SNP_verString '.mat'], 'chr_SNPdata_RADseq');
    else
        fprintf('\nMAT file containing project SNP information found, loading.\n');
        load([workingDir 'matlab_dir/' child_name '.raw_SNP_' SNP_verString '.mat']);
    end;
    chr_SNPdata_RADseq_parent = chr_SNPdata_RADseq;
    clear chr_SNPdata_RADseq;
end;

% basic plot parameters not defined per genome.
TickSize  = -0.005;  %negative for outside, percentage of longest chr figure.
maxY      = 0.5;

%% -----------------------------------------------------------------------------------------
% Make figures
%-------------------------------------------------------------------------------------------
Standard_fig = figure(1);
set(gcf, 'Position', [0 70 1024 600]);
totalSNPdata = [];
for chr = 1:num_chr
    totalSNPdata = [totalSNPdata chr_SNPdata_RADseq_child{chr}];
end;
max_ratio    = max(totalSNPdata);
min_ratio    = min(totalSNPdata);
mean_ratio   = mean(totalSNPdata);
median_ratio = median(totalSNPdata);
fprintf(['\n max ratio    = ' num2str(max_ratio)]);
fprintf(['\n min ratio    = ' num2str(min_ratio)]);
fprintf(['\n mean ratio   = ' num2str(mean_ratio)]);
fprintf(['\n median ratio = ' num2str(median_ratio) '\n\n']);

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
    Linear_TickSize        = -0.01;   %negative for outside, percentage of longest chr figure.
    maxY                   = 0.5;
    
    Linear_left = Linear_left_start;
end;

%% -----------------------------------------------------------------------------------------
% Make figures
%-------------------------------------------------------------------------------------------
plotData = chr_SNPdata_RADseq_child;
for chr = 1:num_chr
    figure(Standard_fig);
    % make standard chr cartoons.
    left   = chr_posX(chr);
    bottom = chr_posY(chr);
    width  = chr_width(chr);
    height = chr_height(chr);
    subplot('Position',[left bottom width height]);
    fprintf(['\nfigposition = [' num2str(left) ' | ' num2str(bottom) ' | ' num2str(width) ' | ' num2str(height) ']']);
    hold on;

    %% ratio data plot section.   
    c_ = [0 0 0];
    fprintf(['\tchr' num2str(chr) ':' num2str(length(plotData{chr}))]);
    plot_X = 1:length(plotData{chr});
    plot_Y = plotData{chr};
    plot_X(plot_Y == 0) = [];
    plot_Y(plot_Y == 0) = [];
    plot(plot_X*chr_length_scale_multiplier,plot_Y,'o','MarkerSize',0.2);
%dddd


    x2 = chr_size(chr)*chr_length_scale_multiplier;
    %% draw lines across plots for easier interpretation of SNP regions.
    line([0 x2], [maxY/4*0 maxY/4*0],'Color',[0.85 0.85 0.85]);
    line([0 x2], [maxY/4*1 maxY/4*1],'Color',[0.85 0.85 0.85]);
    line([0 x2], [maxY/4*2 maxY/4*2],'Color',[0.85 0.85 0.85]);
    line([0 x2], [maxY/4*3 maxY/4*3],'Color',[0.85 0.85 0.85]);
    line([0 x2], [maxY/4*4 maxY/4*4],'Color',[0.85 0.85 0.85]);
    %% end ratio data plot section.

    %axes labels etc.
    hold off;
    xlim([0,chr_size(chr)*chr_length_scale_multiplier]);
    ylim([0,maxY]);
    
    set(gca,'YTick',[]);
    set(gca,'TickLength',[(TickSize*chr_size(largestChr)/chr_size(chr)) 0]); %ensures same tick size on all subfigs.
    ylabel(chr_label{chr}, 'Rotation', 90, 'HorizontalAlign', 'center', 'VerticalAlign', 'bottom');
    set(gca,'XTick',0:(40*5000*chr_length_scale_multiplier):(650*5000*chr_length_scale_multiplier));
    set(gca,'XTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8','2.0','2.2','2.4','2.6','2.8','3.0','3.2'});

    % This section sets the Y-axis labelling.
    set(gca,'YTick',[0 maxY/2 maxY]);
    set(gca,'YTickLabel',{'','',''});
    text(-50000*chr_length_scale_multiplier, 0,   'hom','HorizontalAlignment','right','Fontsize',5);
    text(-50000*chr_length_scale_multiplier, maxY,'1:1','HorizontalAlignment','right','Fontsize',5);

    set(gca,'FontSize',6);
    if (chr == find(chr_posY == max(chr_posY)))
        title([child_name ' SNP map'],'Interpreter','none','FontSize',12);
    end;
    
    hold on;
    %end axes labels etc.

    %show centromere.
    if (chr_size(chr) < 100000)
        Centromere_format = 1;
    else
        Centromere_format = Centromere_format_default;
    end;
    x1 = cen_start(chr)*chr_length_scale_multiplier;
    x2 = cen_end(chr)*chr_length_scale_multiplier;
    leftEnd  = 0.5*5000*chr_length_scale_multiplier;
    rightEnd = chr_size(chr)*chr_length_scale_multiplier-0.5*5000*chr_length_scale_multiplier;
    if (Centromere_format == 0)
        % standard chromosome cartoons in a way which will not cause segfaults when running via commandline.
        dx     = 5*5000*chr_length_scale_multiplier;
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
        plot([leftEnd   leftEnd   leftEnd+dx   x1-dx   x1        x2        x2+dx   rightEnd-dx   rightEnd   rightEnd   rightEnd-dx x2+dx   x2   x1   x1-dx   leftEnd+dx   leftEnd],...
             [dy        maxY-dy   maxY         maxY    maxY-dy   maxY-dy   maxY    maxY          maxY-dy    dy         0           0       dy   dy   0       0            dy     ],...
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
                annotationloc   = (annotation_location(i) - 0.5*5000)*chr_length_scale_multiplier;
		annotationStart = (annotation_start(i)    - 0.5*5000)*chr_length_scale_multiplier;
                annotationEnd   = (annotation_end(i)      - 0.5*5000)*chr_length_scale_multiplier;
                if (strcmp(annotation_type{i},'dot') == 1)
                    plot(annotationloc,-maxY/10*1.5,'k:o','MarkerEdgeColor',annotation_edgecolor{i}, ...
                                                          'MarkerFaceColor',annotation_fillcolor{i}, ...
                                                          'MarkerSize',     annotation_size(i));
                elseif (strcmp(annotation_type{i},'block') == 1)
                    fill([annotationStart annotationStart annotationEnd annotationEnd], ...
                         [-maxY/10*(1.5+0.75) -maxY/10*(1.5-0.75) -maxY/10*(1.5-0.75) -maxY/10*(1.5+0.75)], ...
                         annotation_fillcolor{i},'EdgeColor',annotation_edgecolor{i});
                end;
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
        histAll   = [];
        histAll2  = [];
        smoothed  = [];
        smoothed2 = [];
        for segment = 1:length(chrCopyNum{chr})
            subplot('Position',[(left+chr_width(chr)+0.005)+width*(segment-1) bottom width height]);
        
            % The SNP-histogram values were normalized to a median value of 1.
            for i = round(1+length(plotData{chr})*chr_breaks{chr}(segment)):round(length(plotData{chr})*chr_breaks{chr}(segment+1))
                if (Low_quality_ploidy_estimate == true)
                    histAll{segment}(i) = plotData{chr}(i)*ploidy*ploidyAdjust;
                else
                    histAll{segment}(i) = plotData{chr}(i)*ploidy;
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
            area(smoothed{segment},(1:length(smoothed{segment}))/ploidyBase*2,'FaceColor',[0 0 0]); 
            hold off;
            set(gca,'YTick',[]);
            set(gca,'XTick',[]);
	    xlim([0,1]);
            ylim([0,ploidyBase*2*50]);
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
        else%           chr_string = num2str(chrCopyNum{chr}(1));
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

	%% ratio data plot section.
	c_ = [0 0 0];
	fprintf(['\tLinear chr' num2str(chr) ':' num2str(length(plotData{chr}))]);
	plot_X = 1:length(plotData{chr});
	plot_Y = plotData{chr};
	plot_X(plot_Y == 0) = [];
	plot_Y(plot_Y == 0) = [];
	plot(plot_X*chr_length_scale_multiplier,plot_Y,'o','MarkerSize',0.2);
%dddd

	%show segmental anueploidy breakpoints.
        if (displayBREAKS == true)
            for segment = 2:length(chr_breaks{chr})-1
                bP = chr_breaks{chr}(segment)*length(plotData{chr});
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
        leftEnd  = 0.5*5000*chr_length_scale_multiplier;
        rightEnd = (chr_size(chr) - 0.5*5000)*chr_length_scale_multiplier;
        if (Centromere_format == 0)
            % standard chromosome cartoons in a way which will not cause segfaults when running via commandline.
            dx     = 5*5000*chr_length_scale_multiplier;
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
            
            % draw outlines of chromosome cartoon.   (drawn after horizontal lines to that cartoon edges are not interrupted by
            % horiz lines.
	    plot([leftEnd   leftEnd   leftEnd+dx   x1-dx   x1        x2        x2+dx   rightEnd-dx   rightEnd   rightEnd ...
                  rightEnd-dx   x2+dx   x2   x1   x1-dx   leftEnd+dx  leftEnd],...
                 [dy        maxY-dy   maxY         maxY    maxY-dy   maxY-dy   maxY    maxY          maxY-dy    dy       ...
                  0             0       dy   dy   0       0           dy  ],...
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
                    annotationloc   = (annotation_location(i)- 0.5*5000)*chr_length_scale_multiplier;
		    annotationStart = (annotation_start(i)   - 0.5*5000)*chr_length_scale_multiplier;
                    annotationEnd   = (annotation_end(i)     - 0.5*5000)*chr_length_scale_multiplier;
                    if (strcmp(annotation_type{i},'dot') == 1)
                        plot(annotationloc,-maxY/10*1.5,'k:o','MarkerEdgeColor',annotation_edgecolor{i}, ...
                                                              'MarkerFaceColor',annotation_fillcolor{i}, ...
                                                              'MarkerSize',     annotation_size(i));
                    elseif (strcmp(annotation_type{i},'block') == 1)
                        fill([annotationStart annotationStart annotationEnd annotationEnd], ...
                             [-maxY/10*(1.5+0.75) -maxY/10*(1.5-0.75) -maxY/10*(1.5-0.75) -maxY/10*(1.5+0.75)], ...
                             annotation_fillcolor{i},'EdgeColor',annotation_edgecolor{i});
                    end;
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
        set(gca,'XTick',0:(40*5000):(650*5000));
        set(gca,'XTickLabel',[]);
        if (chr == 1)
            ylabel(child_name, 'Rotation', 0, 'HorizontalAlign', 'right', 'VerticalAlign', 'bottom','Interpreter','none','FontSize',5);
        
            % This section sets the Y-axis labelling.
	    set(gca,'YTick',[0 maxY/2 maxY]);
	    set(gca,'YTickLabel',{'','',''});
	    text(-50000*chr_length_scale_multiplier, 0,   'hom','HorizontalAlignment','right','Fontsize',5);
	    text(-50000*chr_length_scale_multiplier, maxY,'1:1','HorizontalAlignment','right','Fontsize',5);
        else
	    set(gca,'YTick',[]);  
	    set(gca,'YTickLabel',[]);
        end;
        set(gca,'FontSize',6);
        %end final reformatting.
        
        % shift back to main figure generation.
	figure(Standard_fig);
        hold on;
    end;
end;

% Save original arrangement of chromosomes.
saveas(Standard_fig, [figureDir child_name '.raw-SNP-map.1.eps'], 'epsc');
delete(Standard_fig);
        
% Save horizontally arranged chromosomes.
set(Linear_fig,'PaperPosition',[0 0 8 0.62222222]);
saveas(Linear_fig, [figureDir child_name '.raw-SNP-map.2.eps'], 'epsc');
delete(Linear_fig);
        
fprintf('\n');

%% ========================================================================
% end stuff
%==========================================================================
end
