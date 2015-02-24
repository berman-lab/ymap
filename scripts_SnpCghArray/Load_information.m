function [centromere, chrSize, MRS, rDNA] = Load_information()
% Load centromere definition file.
%    This is text file containing one header line and two columns.
%    The two columns hold the start and end bp for the centromeres, with
%       respect to each chromosome.
centromere = [];
fprintf('\n\tLoading centromere definition file.');
centromere_fid = fopen('centromere_locations.txt', 'r');
discard = fgetl(centromere_fid);   clear discard; %discard header line.
lines_analyzed = 0;
while not (feof(centromere_fid))
	line = fgetl(centromere_fid);
	lines_analyzed = lines_analyzed+1;
	cen_chr = sscanf(line, '%s',1);
	cen_start = sscanf(line, '%s',2);
	for i = 1:size(sscanf(line,'%s',1),2);
		cen_start(1) = [];
	end;
	cen_end   = sscanf(line, '%s',3);
	for i = 1:size(sscanf(line,'%s',2),2);
		cen_end(1) = [];
	end;
	chromosome = str2double(cen_chr);
	centromere(chromosome).start = str2double(cen_start);
	centromere(chromosome).end   = str2double(cen_end);
end;
fclose(centromere_fid);
clear cen_start cen_end line lines_analyzed i ans cen_chr centromere_fid chromosome;


% Load chromosome size definition file.
%    This is text file containing one header line and two columns.
%    The two columns hold the start and end bp for the centromeres, with
%       respect to each chromosome.
chrSize = [];
fprintf('\n\tLoading chromosome sizes definition file.');
chrSize_fid = fopen('chromosome_sizes.txt', 'r');
discard = fgetl(chrSize_fid);   clear discard; %discard header line.
lines_analyzed = 0;
while not (feof(chrSize_fid))
	line = fgetl(chrSize_fid);
	lines_analyzed = lines_analyzed+1;
	size_chr = sscanf(line, '%s',1);
	size_size = sscanf(line, '%s',2);
	for i = 1:size(sscanf(line,'%s',1),2);
		size_size(1) = [];
	end;
	size_name   = sscanf(line, '%s',3);
	for i = 1:size(sscanf(line,'%s',2),2);
		size_name(1) = [];
	end;
	chromosome = str2double(size_chr);
	chrSize(chromosome).size = str2double(size_size);
	chrSize(chromosome).name = size_name;
end;
fclose(chrSize_fid);
clear chrSize_fid chromosome i line lines_analyzed size_chr size_name size_size;


% Load MRS location definition file.
%    This is text file containing one header line and four columns.
%    The two columns hold the start and end bp for the centromeres, with
%       respect to each chromosome.
MRS = [];
fprintf('\n\tLoading MRS locations definition file.');
MRS_fid = fopen('MRS_locations.txt', 'r');
discard = fgetl(MRS_fid);   clear discard; %discard header line.
lines_analyzed = 0;
MRS_count = 0;
while not (feof(MRS_fid))
	line = fgetl(MRS_fid);
	lines_analyzed = lines_analyzed+1;
	MRS_chr = sscanf(line, '%s',1);
	MRS_start = sscanf(line, '%s',2);
	for i = 1:size(sscanf(line,'%s',1),2);
		MRS_start(1) = [];
	end;
	MRS_end   = sscanf(line, '%s',3);
	for i = 1:size(sscanf(line,'%s',2),2);
		MRS_end(1) = [];
	end;
	MRS_name   = sscanf(line, '%s',4);
	for i = 1:size(sscanf(line,'%s',3),2);
		MRS_name(1) = [];
	end;
	MRS_count = MRS_count+1;
	MRS(MRS_count).chr        = str2double(MRS_chr);
	MRS(MRS_count).start      = str2double(MRS_start);
	MRS(MRS_count).end        = str2double(MRS_end);
	MRS(MRS_count).name       = MRS_name;
	MRS(MRS_count).color_fill = 'k';
	MRS(MRS_count).color_edge = 'k';
end;
fclose(MRS_fid);


% Load rDNA location definition file.
%    This is text file containing one header line and four columns.
%    The four columns hold the chr, start bp, end bp, and label the rDNA.
rDNA = [];
fprintf('\n\tLoading rDNA locations definition file.');
rDNA_fid = fopen('rDNA_location.txt', 'r');
discard = fgetl(rDNA_fid);   clear discard; %discard header line.
lines_analyzed = 0;
line = fgetl(rDNA_fid);
rDNA_chr = sscanf(line, '%s',1);
rDNA_start = sscanf(line, '%s',2);
for i = 1:size(sscanf(line,'%s',1),2);
	rDNA_start(1) = [];
end;
rDNA_end   = sscanf(line, '%s',3);
for i = 1:size(sscanf(line,'%s',2),2);
	rDNA_end(1) = [];
end;
rDNA_name   = sscanf(line, '%s',4);
for i = 1:size(sscanf(line,'%s',3),2);
	rDNA_name(1) = [];
end;
rDNA.chr        = str2double(rDNA_chr);
rDNA.start      = str2double(rDNA_start);
rDNA.end        = str2double(rDNA_end);
rDNA.name       = rDNA_name;
rDNA.color_fill = 'b';
rDNA.color_edge = 'b';
fclose(rDNA_fid);
end

