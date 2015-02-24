function [centromeres, chrSize, figure_details, annotations, figInfo_ploidy_default] = ...
					Load_genome_information_1(genomeDir,genome)
fprintf(['\nLoad_genome_information_1.m : Genome in use : [' genome ']\n']);

% Load centromere definition file.
%    This is text file containing one header line and two columns.
%    The two columns hold the start and end bp for the centromeres, with
%       respect to each chromosome.
centromeres          = [];
centromere_file_name = [genomeDir '/centromere_locations.txt'];
centromere_fid       = fopen(centromere_file_name, 'r');
discard              = fgetl(centromere_fid);
lines_analyzed       = 0;
fprintf(['\tcentromere_locations.txt\n']);
while not (feof(centromere_fid))
    lineData       = fgetl(centromere_fid);
    lines_analyzed = lines_analyzed+1;
    cen_chr        = sscanf(lineData, '%s',1);
    cen_start      = sscanf(lineData, '%s',2);
    for i = 1:size(sscanf(lineData,'%s',1),2);
        cen_start(1) = [];
    end;
    cen_end   = sscanf(lineData, '%s',3);
    for i = 1:size(sscanf(lineData,'%s',2),2);
        cen_end(1) = [];
    end;
    chr = str2double(cen_chr);
    centromeres(chr).chr   = chr;
    centromeres(chr).start = str2double(cen_start);
    centromeres(chr).end   = str2double(cen_end);
	fprintf(['\t\t|' lineData '\n']);
end;
fclose(centromere_fid);
clear cen_start cen_end line lines_analyzed i ans cen_chr centromere_fid chromosome;
if (length(centromeres) == 0)
	error('[analyze_CNVs]: Centromere definition file is missing.');
end;

% Load chromosome size definition file.
%    This is text file containing one header line and two columns.
%    The two columns hold the start and end bp for the centromeres, with
%       respect to each chromosome.
chrSize        = [];
chrSize_fid    = fopen([genomeDir '/chromosome_sizes.txt'],'r');
discard        = fgetl(chrSize_fid);   clear discard; %discard header line.
lines_analyzed = 0;
fprintf(['\tchromosome_sizes.txt\n']);
while not (feof(chrSize_fid))
    lineData       = fgetl(chrSize_fid);
    lines_analyzed = lines_analyzed+1;
    size_chr       = sscanf(lineData, '%s',1);
    size_size      = sscanf(lineData, '%s',2);
    for i = 1:size(sscanf(lineData,'%s',1),2);
        size_size(1) = [];
    end;
    size_name      = sscanf(lineData, '%s',3);
    for i = 1:size(sscanf(lineData,'%s',2),2);
        size_name(1) = [];
    end;
    chr = str2double(size_chr);
    chrSize(chr).chr  = chr;
    chrSize(chr).size = str2double(size_size);
    chrSize(chr).name = size_name;
	fprintf(['\t\t|' lineData '\n']);
end;
fclose(chrSize_fid);
if (length(chrSize) == 0)
	error('[analyze_CNVs]: Chromosome size definition file is missing.');
end;

% Load additional annotation location definition file.
%    This is text file containing one header line and two columns.
%    The two columns hold the start and end bp for the centromeres, with
%       respect to each chromosome.
annotations       = [];
annotations_fid   = fopen([genomeDir '/annotations.txt'], 'r');
discard           = fgetl(annotations_fid);   clear discard; %discard header line.
lines_analyzed    = 0;
annotations_count = 0;
fprintf(['\tannotations.txt\n']);
while not (feof(annotations_fid))
	lineData = fgetl(annotations_fid);
	if (strcmp(lineData(1),'#') == 0)
		lines_analyzed    = lines_analyzed+1;
		annotations_chr   = sscanf(lineData, '%s',1);
		annotations_type  = sscanf(lineData, '%s',2);
		for i = 1:size(sscanf(lineData,'%s',1),2);
		   annotations_type(1) = [];
		end;
		annotations_start = sscanf(lineData, '%s',3);
		for i = 1:size(sscanf(lineData,'%s',2),2);
		    annotations_start(1) = [];
		end;
		annotations_end   = sscanf(lineData, '%s',4);
		for i = 1:size(sscanf(lineData,'%s',3),2);
		    annotations_end(1) = [];
		end;
		annotations_name  = sscanf(lineData, '%s',5);
		for i = 1:size(sscanf(lineData,'%s',4),2);
		    annotations_name(1) = [];
		end;
		annotations_fillcolor  = sscanf(lineData, '%s',6);
		for i = 1:size(sscanf(lineData,'%s',5),2);
		    annotations_fillcolor(1) = [];
		end;
		annotations_edgecolor  = sscanf(lineData, '%s',7);
		for i = 1:size(sscanf(lineData,'%s',6),2);
		    annotations_edgecolor(1) = [];
		end;
		annotations_size  = sscanf(lineData, '%s',8);
		for i = 1:size(sscanf(lineData,'%s',7),2);
		    annotations_size(1) = [];
		end;

		annotations_count                        = annotations_count+1;
		annotations(annotations_count).chr       = str2double(annotations_chr);
		annotations(annotations_count).type      = annotations_type;
		annotations(annotations_count).start     = str2double(annotations_start);
		annotations(annotations_count).end       = str2double(annotations_end);
		annotations(annotations_count).name      = annotations_name;
		annotations(annotations_count).fillcolor = annotations_fillcolor;
		annotations(annotations_count).edgecolor = annotations_edgecolor;
		annotations(annotations_count).size      = str2double(annotations_size);
		fprintf(['\t\t|' lineData '\n']);
	end;
end;
fclose(annotations_fid);

% Load figure definition file.
% This is text file containing one header line and seven columns.
%    Chr #      : Numerical designations of chromosomes.   (0 is used for line defining figure key.)
%    Chr label  : The label to use for identifying the chromosome in the figure.
%    Chr name   : The full name of the chromosome.
%    Chr posX   : The X-position in % from left to right.
%    Chr posY   : The Y-position in % from bottom to top.
%    Chr width  : The width in %.
%    Chr height : The height in %.
figInfo_ploidy_default = 2.0;
figure_details         = [];
figInfo_fid            = fopen([genomeDir '/figure_definitions.txt'], 'r');
discard                = fgetl(figInfo_fid);   clear discard; %discard header line.
lines_analyzed         = 0;
fprintf(['\tfigure_definitions.txt\n']);
while not (feof(figInfo_fid))
    lineData       = fgetl(figInfo_fid);
    figInfo_chr    = sscanf(lineData, '%s',1);
    figInfo_useChr = sscanf(lineData, '%s',2);
    for i = 1:size(sscanf(lineData,'%s',1),2);
        figInfo_useChr(1) = [];
    end;
    figInfo_label = sscanf(lineData, '%s',3);
    for i = 1:size(sscanf(lineData,'%s',2),2);
        figInfo_label(1) = [];
    end;
	if (str2num(figInfo_chr) > 0) || ...
	   (str2num(figInfo_useChr) > 0) || ...
       ((str2num(figInfo_chr) == 0) && (strcmp(figInfo_label,'Key') == 1)) || ...
       ((str2num(figInfo_chr) == 0) && (strcmp(figInfo_label,'Mito') == 1))
		lines_analyzed = lines_analyzed+1;
		figInfo_name   = sscanf(lineData, '%s',4);
		for i = 1:size(sscanf(lineData,'%s',3),2);
		    figInfo_name(1) = [];
		end;
		figInfo_posX   = sscanf(lineData, '%s',5);
		for i = 1:size(sscanf(lineData,'%s',4),2);
		    figInfo_posX(1) = [];
		end;
		figInfo_posY   = sscanf(lineData, '%s',6);
		for i = 1:size(sscanf(lineData,'%s',5),2);
		    figInfo_posY(1) = [];
		end;
		figInfo_width   = sscanf(lineData, '%s',7);
		for i = 1:size(sscanf(lineData,'%s',6),2);
		    figInfo_width(1) = [];
		end;
		figInfo_height   = sscanf(lineData, '%s',8);
		for i = 1:size(sscanf(lineData,'%s',7),2);
		    figInfo_height(1) = [];
		end;
		figure_details(lines_analyzed).chr    = str2double(figInfo_chr);
		figure_details(lines_analyzed).label  = figInfo_label;
		figure_details(lines_analyzed).name   = figInfo_name;
		figure_details(lines_analyzed).useChr = figInfo_useChr;
		figure_details(lines_analyzed).posX   = str2double(figInfo_posX);
		figure_details(lines_analyzed).posY   = str2double(figInfo_posY);
		figure_details(lines_analyzed).width  = figInfo_width;
		figure_details(lines_analyzed).height = str2double(figInfo_height);
    	elseif ((str2num(figInfo_chr) == 0) && (strcmp(figInfo_label,'Ploidy') == 1))
		figInfo_ploidy_default = sscanf(lineData, '%s',4);
		for i = 1:size(sscanf(lineData,'%s',3),2);
		    figInfo_ploidy_default(1) = [];
		end;
		figInfo_ploidy_default = str2num(figInfo_ploidy_default);
    end;
	fprintf(['\t\t|' lineData '\n']);
end;
if (length(figure_details) == 0)
    error('[analyze_CNVs]: Figure display definition file is missing.');
end;


%% figure out widths for chromosomes in figure.
maxFigSize = 0;
maxChrSize = 0;
for i = 1:length(figure_details)
    fprintf(['Fig_chr : [' num2str(figure_details(i).chr) '|']);
    fprintf([figure_details(i).label '|']);
    fprintf([figure_details(i).name '|']);
    fprintf([num2str(figure_details(i).posX) '|']);
    fprintf([num2str(figure_details(i).posY) '|']);
    fprintf([figure_details(i).width '|']);
    fprintf([num2str(figure_details(i).height) ']\n']);

    if (figure_details(i).chr > 0)
        if (strcmp(figure_details(i).width(1),'*') == 0)
            maxFigSize = str2num(figure_details(i).width);
            maxChrSize = chrSize(figure_details(i).chr).size;
        end;
    end;
end;
for i = 1:length(figure_details)
    if (figure_details(i).chr > 0)
        currentChrSize          = chrSize(figure_details(i).chr).size;
        figure_details(i).width = currentChrSize/maxChrSize*maxFigSize;
    end;
end;
fclose(figInfo_fid);

end
