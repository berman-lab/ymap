function [] = CNV_v6_4(main_dir,user,genomeUser,project,genome,ploidyEstimateString,ploidyBaseString, ...
                       CNV_verString,rDNA_verString,displayBREAKS, referenceCHR);

%% ========================================================================
Centromere_format_default   = 0;
Yscale_nearest_even_ploidy  = true;
HistPlot                    = true;
ChrNum                      = true;
show_annotations            = true;
analyze_rDNA                = true;
Linear_display              = true;
Low_quality_ploidy_estimate = true;


%%=========================================================================
% Load FASTA file name from 'reference.txt' file for project.
%--------------------------------------------------------------------------
userReference    = [main_dir 'users/' user '/genomes/' genome '/reference.txt'];
defaultReference = [main_dir 'users/default/genomes/' genome '/reference.txt'];
if (exist(userReference,'file') == 0)
	FASTA_string = strtrim(fileread(defaultReference));
else
	FASTA_string = strtrim(fileread(userReference));
end;
[FastaPath,FastaName,FastaExt] = fileparts(FASTA_string);


%%=========================================================================
% Control variables.
%--------------------------------------------------------------------------
projectDir = [main_dir 'users/' user '/projects/' project '/'];
genomeDir  = [main_dir 'users/' genomeUser '/genomes/' genome '/'];

fprintf(['\n$$ projectDir : ' projectDir '\n']);
fprintf([  '$$ genomeDir  : ' genomeDir  '\n']);
fprintf([  '$$ genome     : ' genome     '\n']);
fprintf([  '$$ project    : ' project    '\n']);
[centromeres, chr_sizes, figure_details, annotations, ploidy_default] = Load_genome_information_1(projectDir,genomeDir, genome);
Aneuploidy = [];

for i = 1:length(chr_sizes)
	chr_size(i) = 0;
end;
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
		chr_in_use(figure_details(i).chr) = str2num(figure_details(i).useChr);
    end;
end;

num_chrs  = length(chr_size);

% bases_per_bin             = 5000;
bases_per_bin               = max(chr_size)/700;
chr_length_scale_multiplier	= 1/bases_per_bin;

% Defines rDNA region.
if (analyze_rDNA == true)
    if (strcmp(genome,'C_albicans_SC5314_version_A21-s02-m03-r03') == 1)
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

% Sanitize user input of euploid state base for species.
ploidyBase = round(str2num(ploidyBaseString));
if (ploidyBase > 4);   ploidyBase = 4;   end;
if (ploidyBase < 1);   ploidyBase = 1;   end; 
fprintf(['\nEuploid base = "' num2str(ploidyBase) '"\n']);

fprintf(['\nGenerating CNV figure from ''' project ''' sequence data.\n']);

% Initializes vectors used to hold copy number data.
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
	    for j = 1:2       % 2 categories tracked : total read counts in bin and size of bin in base-pairs.
	        chr_CNVdata{chr,j}     = zeros(1,ceil(chr_size(chr)/bases_per_bin));
	    end;
	end;
end;
% Initializes vectors used to hold other site data, for normalization purposes.
for j = 1:2
    nonchr_CNVdata{1,j}        = zeros(1,ceil(chr_size(chr)/bases_per_bin));
end;

if (exist([projectDir 'CNV_' CNV_verString '.mat'],'file') == 0)
    fprintf('\nMAT file not found, regenerating.\n');
    datafile = [projectDir 'putative_CNVs_' CNV_verString '.txt'];
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
        % count of data points at locus.
        
        if (chr > 0)
	    % Current data point belongs to a used chromosome.
        if (val <= length(chr_CNVdata{chr,1}))
            chr_CNVdata{chr,1}(val) = chr_CNVdata{chr,1}(val)+str2num(count);
            chr_CNVdata{chr,2}(val) = chr_CNVdata{chr,2}(val)+1;
	   	else
			error_counter_2 = error_counter_2+1;
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
    
    save([projectDir 'CNV_' CNV_verString '.mat'],'chr_CNVdata');
    if (analyze_rDNA == 1)
        save([projectDir 'rDNA-CNV_' rDNA_verString '.mat'],'countRDNA','sizeRDNA','countREF','sizeREF');
    end;
else
    fprintf('\nMAT file found, loading.\n');
    load([projectDir 'CNV_' CNV_verString '.mat']);
end;

% basic plot parameters not defined per genome.
TickSize        = -0.005;  %negative for outside, percentage of longest chr figure.
maxY            = ploidyBase*2;


%% -----------------------------------------------------------------------------------------
% Setup for main figure generation.
%-------------------------------------------------------------------------------------------
fig = figure(1);
set(gcf, 'Position', [0 70 1024 600]);

% calculate CGH bin values.
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
	    for pos = 1:length(chr_CNVdata{chr,1})											% Add to the plotting array : the average of the data points in the region.
	        if (chr_CNVdata{chr,2}(pos) == 0)
	            CNVplot{chr}(pos) = 0;													% If there are no data elements in the region, assign a null value.
	        else
				CNVplot{chr}(pos) = chr_CNVdata{chr,1}(pos)/chr_CNVdata{chr,2}(pos);	% if there are data elements in the region, divide the sum by the number of data elements.
	        end;
	    end;
	end;
end;

% Gather CGH and GCratio data for LOWESS fitting.
CGHdata_all          = [];
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
	    CGHdata_all      = [CGHdata_all     CNVplot{chr}];
	end;
end;
medianRawY = median(CGHdata_all);

for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		CNVplot{chr} = CNVplot{chr}/medianRawY;
	end;
end;


%% ====================================================================
% Apply GC bias correction to the CGH data.
%----------------------------------------------------------------------
% Load standard bin GC_bias data from : standard_bins.GC_ratios.txt
fprintf(['standard_bins_GC_ratios_file :\n\t' main_dir 'users/' genomeUser '/genomes/' genome '/' FastaName '.GC_ratios.standard_bins.txt\n']);
standard_bins_GC_ratios_fid = fopen([main_dir 'users/' genomeUser '/genomes/' genome '/' FastaName '.GC_ratios.standard_bins.txt'], 'r');
fprintf(['\t' num2str(standard_bins_GC_ratios_fid) '\n']);
lines_analyzed = 0;
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
	    chr_GCratioData{chr} = zeros(1,ceil(chr_size(chr)/bases_per_bin));
	end;
end;
while not (feof(standard_bins_GC_ratios_fid))
    dataLine = fgetl(standard_bins_GC_ratios_fid);
	if (length(dataLine) > 0)
		if (dataLine(1) ~= '#')
			lines_analyzed = lines_analyzed+1;
			chr            = str2num(sscanf(dataLine, '%s',1));
			fragment_start = sscanf(dataLine, '%s',2);  for i = 1:size(sscanf(dataLine,'%s',1),2);      fragment_start(1) = []; end;    fragment_start = str2num(fragment_start);
			fragment_end   = sscanf(dataLine, '%s',3);  for i = 1:size(sscanf(dataLine,'%s',2),2);      fragment_end(1) = [];   end;    fragment_end   = str2num(fragment_end);
			GCratio        = sscanf(dataLine, '%s',4);  for i = 1:size(sscanf(dataLine,'%s',3),2);      GCratio(1) = [];        end;    GCratio        = str2num(GCratio);
			position       = ceil(fragment_start/bases_per_bin);
			chr_GCratioData{chr}(position) = GCratio;
		end;
	end;
end;
fclose(standard_bins_GC_ratios_fid);
% Gather GCratio data for LOWESS fitting.
GCratioData_all = [];
CGHdata_all     = [];
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
	    GCratioData_all  = [GCratioData_all chr_GCratioData{chr}];
		CGHdata_all      = [CGHdata_all     CNVplot{chr}        ];
	end;
end;

% Clean up data by:
%    deleting GC ratio data near zero.
%    deleting CGH data beyond 3* the median value.  (rDNA, etc.)
CGHdata_clean                               = CGHdata_all;
GCratioData_clean                           = GCratioData_all;

%CGHdata_clean(GCratioData_clean < 0.01)     = [];
%GCratioData_clean(GCratioData_clean < 0.01) = [];
%GCratioData_clean(CGHdata_clean > 3)        = [];
%CGHdata_clean(CGHdata_clean > 3)            = [];

% Perform LOWESS fitting.
rawData_X1        = GCratioData_clean;
rawData_Y1        = CGHdata_clean;
fprintf(['Lowess X:Y size : [' num2str(size(rawData_X1,1)) ',' num2str(size(rawData_X1,2)) ']:[' num2str(size(rawData_Y1,1)) ',' num2str(size(rawData_Y1,2)) ']\n']);
numFits           = 10;
[fitX1, fitY1]    = optimize_mylowess(rawData_X1,rawData_Y1, numFits);

% Correct data using normalization to LOWESS fitting
Y_target = 1;
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		fprintf(['chr' num2str(chr) ' : ' num2str(length(chr_GCratioData{chr})) ' ... ' num2str(length(CNVplot{chr})) '\n']);
		fprintf(['\tnumbins = ' num2str(ceil(chr_size(chr)/bases_per_bin)) '\n']);
		rawData_chr_X{chr}        = chr_GCratioData{chr};
		rawData_chr_Y{chr}        = CNVplot{chr};
		fitData_chr_Y{chr}        = interp1(fitX1,fitY1,rawData_chr_X{chr},'spline');
		normalizedData_chr_Y{chr} = rawData_chr_Y{chr}./fitData_chr_Y{chr}*Y_target;
	end;
end;

%% Generate figure showing subplots of LOWESS fittings.
GCfig = figure();
subplot(1,2,1);
	hold on;
	for chr = 1:num_chrs
		if (chr_in_use(chr) == 1)
			plot(rawData_chr_X{chr},rawData_chr_Y{chr},'k.','markersize',1);		% raw data
		end;
	end;
	plot(fitX1,fitY1,'r','LineWidth',2);						% LOWESS fit curve.
    hold off;
    xlabel('GC ratio');
    ylabel('CGH data');
	xlim([0.0 1.0]);
	ylim([0 2]);
    axis square;
subplot(1,2,2);
    hold on;
	for chr = 1:num_chrs
		if (chr_in_use(chr) == 1)
			plot(rawData_chr_X{chr},normalizedData_chr_Y{chr},'k.','markersize',1);	% corrected data.
		end;
	end;
	plot([0.0 1.0],[Y_target Y_target],'r','LineWidth',2);			% normalization line.
    hold off;
    xlabel('GC ratio');
    ylabel('corrected CGH data');
	xlim([0.0 1.0]);
	ylim([0 2]);
    axis square;
saveas(GCfig, [projectDir 'fig.GCratio_vs_CGH.png'], 'png');

% Move LOWESS-normalizd CGH data into display pipeline.
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		% GC-corrected.
		CNVplot{chr} = normalizedData_chr_Y{chr};
	end;
end;

%% ====================================================================
% Save presented CNV data in a file format common across pipeline modules.
%----------------------------------------------------------------------
fprintf('\nSaving "Common_CNV" data file.');
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		CNVplot2{chr} = CNVplot{chr};
	end;
end;
genome_CNV = genome;
save([projectDir 'Common_CNV.mat'], 'CNVplot2','genome_CNV');

ploidy = str2num(ploidyEstimateString);
[chr_breaks, chrCopyNum, ploidyAdjust] = FindChrSizes_2(Aneuploidy,CNVplot2,ploidy,num_chrs);
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
for chr = 1:num_chrs
    if (chr_in_use(chr) == 1)
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
	    % The ratio of 'ploidy' to 'ploidyBase' determines where the data is displayed relative to the 
	    % median line.
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
				set(gca,'YTickLabel',{'','',''});
				text(-50000*chr_length_scale_multiplier, maxY/2,   '1','HorizontalAlignment','right','Fontsize',5);
				text(-50000*chr_length_scale_multiplier, maxY,     '2','HorizontalAlignment','right','Fontsize',5);
		    case 2
				set(gca,'YTick',[0 maxY/4 maxY/2 maxY/4*3 maxY]);
				set(gca,'YTickLabel',{'','','','',''});
				text(-50000*chr_length_scale_multiplier, maxY/4,   '1','HorizontalAlignment','right','Fontsize',5);
				text(-50000*chr_length_scale_multiplier, maxY/2,   '2','HorizontalAlignment','right','Fontsize',5);
				text(-50000*chr_length_scale_multiplier, maxY/4*3, '3','HorizontalAlignment','right','Fontsize',5);
				text(-50000*chr_length_scale_multiplier, maxY,     '4','HorizontalAlignment','right','Fontsize',5);
		    case 3
				set(gca,'YTick',[0 maxY/6 maxY/3 maxY/2 maxY/3*2 maxY/6*5 maxY]);
				set(gca,'YTickLabel',{'','','','','','',''});
				text(-50000*chr_length_scale_multiplier, maxY/2,   '3','HorizontalAlignment','right','Fontsize',5);
				text(-50000*chr_length_scale_multiplier, maxY,     '6','HorizontalAlignment','right','Fontsize',5);
		    case 4
				set(gca,'YTick',[0 maxY/8 maxY/4 maxY/8*3 maxY/2 maxY/8*5 maxY/4*3 maxY/8*7 maxY]);
				set(gca,'YTickLabel',{'','','','','','','','',''});
				text(-50000*chr_length_scale_multiplier, maxY/4,   '2','HorizontalAlignment','right','Fontsize',5);
				text(-50000*chr_length_scale_multiplier, maxY/2,   '4','HorizontalAlignment','right','Fontsize',5);
				text(-50000*chr_length_scale_multiplier, maxY/4*3, '6','HorizontalAlignment','right','Fontsize',5);
				text(-50000*chr_length_scale_multiplier, maxY,     '8','HorizontalAlignment','right','Fontsize',5);
		end;

        set(gca,'FontSize',6);
        if (chr == find(chr_posY == max(chr_posY)))
            title([ project ' CNV map'],'Interpreter','none','FontSize',12);
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

            % draw outlines of chromosome cartoon.   (drawn after horizontal lines to that cartoon edges are 
	    % not interrupted by horiz lines.
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
                    annotationLoc   = annotation_location(i)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
		    annotationStart = annotation_start(i)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
		    annotationEnd   = annotation_end(i)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
		    if (strcmp(annotation_type{i},'dot') == 1)
                        plot(annotationLoc,-maxY/10*1.5,'k:o','MarkerEdgeColor',annotation_edgecolor{i}, ... 
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
            chr_CNVdata;
            histAll   = [];
            histAll2  = [];
            smoothed  = [];
            smoothed2 = [];
			fprintf(['chr = ' num2str(chr) '\n']);
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
            else
		chr_string = num2str(chrCopyNum{chr}(1));
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
                        annotationLoc   = annotation_location(i)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
			annotationStart = annotation_start(i)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
                        annotationEnd   = annotation_end(i)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
			if (strcmp(annotation_type{i},'dot') == 1)
                            plot(annotationLoc,-maxY/10*1.5,'k:o','MarkerEdgeColor',annotation_edgecolor{i}, ...
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
            set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
            set(gca,'XTickLabel',[]);
            if (chr == 1)
                ylabel(project, 'Rotation', 0, 'HorizontalAlign', 'right', 'VerticalAlign', 'bottom','Interpreter','none','FontSize',5);

	        % This section sets the Y-axis labelling.
		switch ploidyBase
		case 1
		    set(gca,'YTick',[0 maxY/2 maxY]);
		    set(gca,'YTickLabel',{'','',''});
		    text(-50000*chr_length_scale_multiplier, maxY/2,   '1','HorizontalAlignment','right','Fontsize',5);
		    text(-50000*chr_length_scale_multiplier, maxY,     '2','HorizontalAlignment','right','Fontsize',5);
		case 2
		    set(gca,'YTick',[0 maxY/4 maxY/2 maxY/4*3 maxY]);
		    set(gca,'YTickLabel',{'','','','',''});
		    text(-50000*chr_length_scale_multiplier, maxY/4,   '1','HorizontalAlignment','right','Fontsize',5);
		    text(-50000*chr_length_scale_multiplier, maxY/2,   '2','HorizontalAlignment','right','Fontsize',5);
		    text(-50000*chr_length_scale_multiplier, maxY/4*3, '3','HorizontalAlignment','right','Fontsize',5);
		    text(-50000*chr_length_scale_multiplier, maxY,     '4','HorizontalAlignment','right','Fontsize',5);
		case 3
		    set(gca,'YTick',[0 maxY/6 maxY/3 maxY/2 maxY/3*2 maxY/6*5 maxY]);
		    set(gca,'YTickLabel',{'','','','','','',''});
		    text(-50000*chr_length_scale_multiplier, maxY/2,   '3','HorizontalAlignment','right','Fontsize',5);
		    text(-50000*chr_length_scale_multiplier, maxY,     '6','HorizontalAlignment','right','Fontsize',5);
		case 4
		    set(gca,'YTick',[0 maxY/8 maxY/4 maxY/8*3 maxY/2 maxY/8*5 maxY/4*3 maxY/8*7 maxY]);
		    set(gca,'YTickLabel',{'','','','','','','','',''});
		    text(-50000*chr_length_scale_multiplier, maxY/4,   '2','HorizontalAlignment','right','Fontsize',5);
		    text(-50000*chr_length_scale_multiplier, maxY/2,   '4','HorizontalAlignment','right','Fontsize',5);
		    text(-50000*chr_length_scale_multiplier, maxY/4*3, '6','HorizontalAlignment','right','Fontsize',5);
		    text(-50000*chr_length_scale_multiplier, maxY,     '8','HorizontalAlignment','right','Fontsize',5);
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
end;

% Save primary genome figure.
set(fig,'PaperPosition',[0 0 8 6]*2);
saveas(fig,        [projectDir 'fig.CNV-map.1.eps'], 'epsc');
saveas(fig,        [projectDir 'fig.CNV-map.1.png'], 'png');
delete(fig);

% Save horizontal aligned genome figure.
set(Linear_fig,'PaperPosition',[0 0 8 0.62222222]*2);
saveas(Linear_fig, [projectDir 'fig.CNV-map.2.eps'], 'epsc');
saveas(Linear_fig, [projectDir 'fig.CNV-map.2.png'], 'png');
delete(Linear_fig);

% Output chromosome copy number estimates.
textFileName = [projectDir 'txt.CNV-map.3.txt'];
fprintf(['Text output of CNVs : "' textFileName '"\n']);
textFileID = fopen(textFileName,'w');
fprintf(textFileID,stringChrCNVs);
fclose(textFileID);

end
