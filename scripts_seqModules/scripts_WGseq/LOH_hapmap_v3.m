function [] = LOH_hapmap_v3(main_dir,user,genomeUser,project,hapmap,genome,ploidyEstimateString,ploidyBaseString, ...
                            SNP_verString,LOH_verString,CNV_verString,displayBREAKS);
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
% Control variables for Candida albicans SC5314.
%--------------------------------------------------------------------------
projectDir  = [main_dir 'users/' user '/projects/' project '/'];

if (exist([[main_dir 'users/default/hapmaps/' hapmap '/']],'dir') == 7)
	hapmapDir = [main_dir 'users/default/hapmaps/' hapmap '/'];
elseif (exist([[main_dir 'users/' user '/hapmaps/' hapmap '/']],'dir') == 7)
	hapmapDir = [main_dir 'users/' user '/hapmaps/' hapmap '/'];
else
	hapmapDir = [main_dir 'users/' user '/projects/' project '/'];
end;

genomeDir  = [main_dir 'users/' genomeUser '/genomes/' genome '/'];

[centromeres, chr_sizes, figure_details, annotations, ploidy_default] = Load_genome_information_1(genomeDir, genome);
[Aneuploidy]                                                          = Load_dataset_information_1(projectDir, project);

num_chrs = length(chr_sizes);

for i = 1:length(chr_sizes)
	chr_size(i)  = 0;
	cen_start(i) = 0;
	cen_end(i)   = 0;
end;
for i = 1:length(chr_sizes)
	chr_size(chr_sizes(i).chr)    = chr_sizes(i).size;
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
num_chrs = length(chr_size);

%% This block is normally calculated in FindChrSizes_2 in CNV analysis.
for usedChr = 1:num_chrs
	if (chr_in_use(usedChr) == 1)
		% determine where the endpoints of ploidy segments are.
		chr_breaks{usedChr}(1) = 0.0;
		break_count = 1;
		if (length(Aneuploidy) > 0)
			for i = 1:length(Aneuploidy)
				if (Aneuploidy(i).chr == usedChr)
					break_count = break_count+1;
					chr_broken = true;
					chr_breaks{usedChr}(break_count) = Aneuploidy(i).break;
				end;
			end;
		end;
		chr_breaks{usedChr}(length(chr_breaks{usedChr})+1) = 1;
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

% basic plot parameters not defined per genome.
TickSize         = -0.005;  %negative for outside, percentage of longest chr figure.
bases_per_bin    = max(chr_size)/700;
maxY             = ploidyBase*2;
cen_tel_Xindent  = 5;
cen_tel_Yindent  = maxY/5;

%define colors for colorBars plot
colorNoData = [1.0   1.0   1.0  ]; %used when no data is available for the bin.
colorInit   = [0.5   0.5   0.5  ]; %external; used in blending at ends of chr.
colorHET    = [0.0   0.0   0.0  ]; % near 1:1 ratio SNPs
colorOddHET = [0.0   1.0   0.0  ]; % Het, but not near 1:1 ratio SNPs.
colorHOM    = [1.0   0.0   0.0  ]; % Hom SNPs;

colorAB     = [0.667 0.667 0.667]; % heterozygous.
colorA      = [1.0   0.0   1.0  ]; % homozygous a:magenta.
colorB      = [0.0   1.0   1.0  ]; % homozygous b:cyan.
%colorAAB    = [0.667 0.333 1.0  ]; % heterozygous trisomy, aab.
%colorABB    = [0.333 0.667 1.0  ]; % heterozygous trisomy, abb.
%colorAAAB   = [0.75   0.25   1.0]; % heterozygous tetrasomy, aaab.
%colorABBB   = [0.25   0.75   1.0]; % heterozygous tetrasomy, abbb.
%colorAAAAB  = [0.8    0.2    1.0]; % heterozygous pentasomy, aaaab.
%colorAAABB  = [0.6    0.4    1.0]; % heterozygous pentasomy, aaabb.
%colorAABBB  = [0.4    0.6    1.0]; % heterozygous pentasomy, aabbb.
%colorABBBB  = [0.2    0.8    1.0]; % heterozygous pentasomy, abbbb.
%colorAAAAAB = [0.8333 0.1667 1.0]; % heterozygous hexasomy, aaaaab.
%colorABBBBB = [0.1667 0.8333 1.0]; % heterozygous hexasomy, abbbbb.
%colorPeak   = [0.5   0.5   0.5  ]; % peak locations.
%colorCutoff = [1.0   0.0   0.0  ]; % cutoff locations.

fprintf(['\nGenerating LOH-map figure from ''' project ''' vs. (hapmap)''' hapmap ''' data.\n']);

% Initializes vectors used to hold copy number data.
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		% 1 category tracked : average read counts per bin.
		chr_CNVdata{chr} = zeros(1,ceil(chr_size(chr)/bases_per_bin));
		% fprintf(['0|' num2str(chr) ':' num2str(length(chr_CNVdata{chr})) '\n']);
	end;
end;

% Initializes vectors used to hold number of SNPs in each interpretation catagory for each chromosome region.
for chr = 1:length(chr_sizes)
	% 4 SNP interpretation catagories tracked.
	%	1 : phased ratio data.
	%	2 : unphased ratio data.
	%   3 : phased coordinate data.
	%   4 : unphased coordinate data.
	chr_length = ceil(chr_size(chr)/bases_per_bin);
	for j = 1:4
		chr_SNPdata{chr,j} = cell(1,chr_length);
	end;
	% fprintf(['0|' num2str(chr) ':' num2str(length(chr_SNPdata{chr,1})) '\n']);
end;


%%================================================================================================
% Load SNP/LOH data.
%-------------------------------------------------------------------------------------------------
if (exist([projectDir 'SNP_' SNP_verString '.mat'],'file') == 0)
	fprintf('\nMAT file not found, regenerating.\n');
	datafile = [projectDir 'preprocessed_SNPs.txt'];
	data     = fopen(datafile, 'r');
	lines_analyzed = 0;
	while not (feof(data))
		dataLine = fgetl(data);
		if (length(dataLine) > 0)
			if (dataLine(1) ~= '#')
				% process the loaded line into data channels.
				lines_analyzed              = lines_analyzed+1;
				chr_num                     = sscanf(dataLine, '%s',1);
				fragment_start              = sscanf(dataLine, '%s',2);   for i = 1:size(sscanf(dataLine,'%s',1),2);   fragment_start(1)              = [];   end;
				fragment_end                = sscanf(dataLine, '%s',3);   for i = 1:size(sscanf(dataLine,'%s',2),2);   fragment_end(1)                = [];   end;
				phased_data_string          = sscanf(dataLine, '%s',4);   for i = 1:size(sscanf(dataLine,'%s',3),2);   phased_data_string(1)          = [];   end;
				unphased_data_string        = sscanf(dataLine, '%s',5);   for i = 1:size(sscanf(dataLine,'%s',4),2);   unphased_data_string(1)        = [];   end;
				phased_coordinates_string   = sscanf(dataLine, '%s',6);   for i = 1:size(sscanf(dataLine,'%s',5),2);   phased_coordinates_string(1)   = [];   end;
				unphased_coordinates_string = sscanf(dataLine, '%s',7);   for i = 1:size(sscanf(dataLine,'%s',6),2);   unphased_coordinates_string(1) = [];   end;

				% format = simple, one number per column.
				chr_num                     = str2num(chr_num);
				fragment_start              = str2num(fragment_start);
				fragment_end                = str2num(fragment_end);
				chr_length                  = ceil(chr_size(chr_num)/bases_per_bin);
				chr_bin                     = ceil(fragment_start/bases_per_bin);
				% fprintf(['chr_num' num2str(chr_num) '|chr_bin = ' num2str(chr_bin) '\n']);

				% format = '(number1,number2,...,numberN)'
				phased_data_string(1)       = [];
				phased_data_string(end)     = [];
				if (length(phased_data_string) == 0)
					phased_data             = [];
				else
					commaCount              = length(find(phased_data_string==','));
					if (commaCount == 0)
						phased_data         = str2num(phased_data_string);
					else
						phased_data         = strsplit(phased_data_string,',');   % function converts number lists from strings to numbers.
					end;
				end;

				% format = '(number1,number2,...,numberN)'
                phased_coordinates_string(1)       = [];
                phased_coordinates_string(end)     = [];
				if (length(phased_coordinates_string) == 0)
					phased_coordinates             = [];
				else
					commaCount              = length(find(phased_coordinates_string==','));
					if (commaCount == 0)
						phased_coordinates         = str2num(phased_coordinates_string);
					else
						phased_coordinates         = strsplit(phased_coordinates_string,',');   % function converts number lists from strings to numbers.
					end;
				end;

				% format = '(number1,number2,...,numberN)'
				unphased_data_string(1)     = [];
				unphased_data_string(end)   = [];
				if (length(unphased_data_string) == 0)
					unphased_data           = [];
				else
					commaCount              = length(find(unphased_data_string==','));
					if (commaCount == 0)
						unphased_data       = str2num(unphased_data_string);
					else
						unphased_data       = strsplit(unphased_data_string,',');  % function converts number lists from strings to numbers.
					end;
				end;

				% format = '(number1,number2,...,numberN)'
				unphased_coordinates_string(1)       = [];
				unphased_coordinates_string(end)     = [];
				if (length(unphased_coordinates_string) == 0)
					unphased_coordinates             = [];
				else
					commaCount              = length(find(unphased_coordinates_string==','));
					if (commaCount == 0)
						unphased_coordinates         = str2num(unphased_coordinates_string);
					else
						unphased_coordinates         = strsplit(unphased_coordinates_string,',');   % function converts number lists from strings to numbers.
					end;
				end;

				% add phased and unphased data to storage arrays.
				chr_SNPdata{chr_num,1}{chr_bin} = phased_data;
				chr_SNPdata{chr_num,2}{chr_bin} = unphased_data;

				% add phased and unphased data coordinates to storage arrays.
				chr_SNPdata{chr_num,3}{chr_bin} = phased_coordinates;
				chr_SNPdata{chr_num,4}{chr_bin} = unphased_coordinates;

				% fprintf(['chr' num2str(chr_num) '|chr_bin = ' num2str(chr_bin) '|' num2str(length(chr_SNPdata{chr_num,1})) '\n']);
			end;
		end;
	end;
	fclose(data);

	save([projectDir 'SNP_' SNP_verString '.mat'],'chr_SNPdata');
else
	fprintf('\nMAT file found, loading.\n');
	load([projectDir 'SNP_' SNP_verString '.mat']);
end;


%% -----------------------------------------------------------------------------------------
% Setup for figure generation.
%-------------------------------------------------------------------------------------------
fig = figure(1);
set(gcf, 'Position', [0 70 1024 600]);

% calculate SNP bin values.
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		for chr_bin = 1:length(chr_SNPdata{chr,1})
			SNPplot{chr,1}{chr_bin} = chr_SNPdata{chr,1}{chr_bin};
			SNPplot{chr,2}{chr_bin} = chr_SNPdata{chr,2}{chr_bin};
		end;
	end;
end;

SNPdata_all = [];
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		for chr_bin = 1:length(SNPplot{chr,1})
			TOTplot{chr}{chr_bin} = length(chr_SNPdata{chr,1}{chr_bin}) + length(chr_SNPdata{chr,2}{chr_bin});   % TOT = phased+nonphased;
			SNPdata_all           = [SNPdata_all TOTplot{chr}{chr_bin}];
		end;
	end;
end;
medianRawY = median(SNPdata_all);

% Gather median-normalized SNP data for LOWESS fitting.
SNPdata_all = [];
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		for chr_bin = 1:length(SNPplot{chr,1})
			TOTplot{chr}{chr_bin} = length(chr_SNPdata{chr,1}{chr_bin}/medianRawY) + length(chr_SNPdata{chr,2}{chr_bin}/medianRawY);   % TOT = phased+nonphased;
			SNPdata_all           = [SNPdata_all TOTplot{chr}{chr_bin}];
		end;
	end;
end;


%% ====================================================================
% Apply GC bias correction to SNP data.
%   Number of putative SNPs vs. GCbias per standard bin.
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

% Gather SNP and GCratio data for LOWESS fitting.
GCratioData_all = [];
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		GCratioData_all = [GCratioData_all chr_GCratioData{chr}];
	end;
end;
medianRawY = median(SNPdata_all);
fprintf(['medianRawY = ' num2str(medianRawY) '\n']);
fprintf(['SNPdata_all     => ' num2str(length(SNPdata_all)    ) '\n']);
fprintf(['GCratioData_all => ' num2str(length(GCratioData_all)) '\n']);

%% Clean up data by:
%%    deleting GC ratio data near zero.
%%    deleting CGH data beyond 3* the median value.  (rDNA, etc.)
SNPdata_clean                                            = SNPdata_all;
GCratioData_clean                                        = GCratioData_all;
SNPdata_clean(     GCratioData_clean < 0.01            ) = [];
GCratioData_clean( GCratioData_clean < 0.01            ) = [];
GCratioData_clean( SNPdata_clean > max(medianRawY*3,3) ) = [];
SNPdata_clean(     SNPdata_clean > max(medianRawY*3,3) ) = [];
GCratioData_clean( SNPdata_clean == 0                  ) = [];
SNPdata_clean(     SNPdata_clean == 0                  ) = [];

% Perform LOWESS fitting.
rawData_X1     = GCratioData_clean;
rawData_Y1     = SNPdata_clean;
fprintf(['Lowess X:Y size : [' num2str(size(rawData_X1,1)) ',' num2str(size(rawData_X1,2)) ']:[' num2str(size(rawData_Y1,1)) ',' num2str(size(rawData_Y1,2)) ']\n']);
[fitX1, fitY1] = optimize_mylowess_SNP(rawData_X1,rawData_Y1);

% Correct data using normalization to LOWESS fitting
Y_target = 1;
% Initialize data structures, to simplify debugging.   These nested loops can be deleted without problems arising.
for chr = 1:num_chrs
	for j = 1:2
		rawData_chr_Y{chr,j}        = [];
		normalizedData_chr_Y{chr,j} = [];
		rawData_chr_X{chr}          = [];
		rawDataAll_chr_Y{chr}       = [];
		fitDataAll_chr_Y{chr}       = [];
		fitData_chr_Y{chr}          = [];
	end;
end;
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		rawData_chr_X{chr}           = chr_GCratioData{chr};
		rawDataAll_chr_Y{chr}        = cell2mat(TOTplot{chr});
		fitDataAll_chr_Y{chr}        = interp1(fitX1,fitY1,rawData_chr_X{chr},'spline');
		normalizedDataAll_chr_Y{chr} = rawDataAll_chr_Y{chr}./fitDataAll_chr_Y{chr}*Y_target;
		fitData_chr_Y{chr}           = interp1(fitX1,fitY1,rawData_chr_X{chr},'spline');

 		for chr_bin = 1:length(SNPplot{chr,1})
			rawData_chr_Y{chr,1}{chr_bin}        = length(chr_SNPdata{chr,1}{chr_bin});
			rawData_chr_Y{chr,2}{chr_bin}        = length(chr_SNPdata{chr,2}{chr_bin});

			normalizedData_chr_Y{chr,1}{chr_bin} = length(rawData_chr_Y{chr,1}{chr_bin})./fitData_chr_Y{chr}*Y_target;
			normalizedData_chr_Y{chr,2}{chr_bin} = length(rawData_chr_Y{chr,2}{chr_bin})./fitData_chr_Y{chr}*Y_target;
		end;
	end;
end;

% Move LOWESS-normalizd SNP data back into display pipeline.
SNPdata_all    = [];
SNPdata_all_1  = [];
SNPdata_all_2  = [];
cSNPdata_all   = [];
cSNPdata_all_1 = [];
cSNPdata_all_2 = [];
GCdata_all     = [];
GCdata_all_1   = [];
GCdata_all_2   = [];
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		TOTplot{chr}        = rawDataAll_chr_Y{chr};
		cTOTplot{chr}       = normalizedDataAll_chr_Y{chr};
		SNPdata_all         = [SNPdata_all    TOTplot{chr}              ];
		cSNPdata_all        = [cSNPdata_all   cTOTplot{chr}             ];
		GCdata_all          = [GCdata_all     rawData_chr_X{chr}        ];

		for chr_bin = 1:length(SNPplot{chr,1})
% darren : attempt at data normalization is failing, introduces major bug due to data compaction.
%			chr_SNPdata{chr,1}{chr_bin}  = rawData_chr_Y{chr,1}{chr_bin};
%			chr_SNPdata{chr,2}{chr_bin}  = rawData_chr_Y{chr,2}{chr_bin};
			cchr_SNPdata{chr,1}{chr_bin} = normalizedData_chr_Y{chr,1}{chr_bin};
			cchr_SNPdata{chr,2}{chr_bin} = normalizedData_chr_Y{chr,2}{chr_bin};
			SNPdata_all_1                = [SNPdata_all_1  length(chr_SNPdata{chr,1}{chr_bin})  ];
			SNPdata_all_2                = [SNPdata_all_2  length(chr_SNPdata{chr,2}{chr_bin})  ];
			cSNPdata_all_1               = [cSNPdata_all_1 length(cchr_SNPdata{chr,1}{chr_bin}) ];
			cSNPdata_all_2               = [cSNPdata_all_2 length(cchr_SNPdata{chr,2}{chr_bin}) ];
		end;
	end;
end;


%% Generate figure showing subplots of LOWESS fittings.
GCfig = figure(3);
subplot(2,3,1);
    plot(GCratioData_all,SNPdata_all,'k.','markersize',1);
    hold on;	plot(fitX1,fitY1,'r','LineWidth',2);   hold off;
    xlabel('GC ratio');   ylabel('SNP data');
    xlim([0.0 1.0]);      ylim([0 max(medianRawY*5,5)]);   axis square;
subplot(2,3,2);
	plot(GCdata_all,SNPdata_all_1,'r.','markersize',1);
	hold on;    plot(fitX1,fitY1,'k','LineWidth',2);   hold off;
	xlabel('GC ratio');   ylabel('SNP data');
	xlim([0.0 1.0]);      ylim([0 max(medianRawY*5,5)]);   axis square;
subplot(2,3,3);
	plot(GCdata_all,SNPdata_all_2,'g.','markersize',1);
	hold on;    plot(fitX1,fitY1,'k','LineWidth',2);   hold off;
	xlabel('GC ratio');   ylabel('SNP data');
	xlim([0.0 1.0]);      ylim([0 max(medianRawY*5,5)]);   axis square;

subplot(2,3,4);
	plot(GCratioData_all,cSNPdata_all,'k.','markersize',1);
	hold on;   plot([min(GCratioData_all) max(GCratioData_all)],[Y_target Y_target],'r','LineWidth',2);   hold off;
	xlabel('GC ratio');   ylabel('corrected SNP data');
	xlim([0.0 1.0]);      ylim([0 5]);                    axis square;
subplot(2,3,5);
	plot(GCdata_all,cSNPdata_all_1,'r.','markersize',1);
	hold on;   plot([min(GCratioData_all) max(GCratioData_all)],[Y_target Y_target],'k','LineWidth',2);   hold off;
	xlabel('GC ratio');   ylabel('corrected SNP data');
	xlim([0.0 1.0]);      ylim([0 5]);                    axis square;
subplot(2,3,6);
	plot(GCdata_all,cSNPdata_all_1,'g.','markersize',1);
	hold on;   plot([min(GCratioData_all) max(GCratioData_all)],[Y_target Y_target],'k','LineWidth',2);   hold off;
	xlabel('GC ratio');   ylabel('corrected SNP data');
	xlim([0.0 1.0]);      ylim([0 5]);                    axis square;

saveas(GCfig, [projectDir '/fig.GCratio_vs_SNP.png'], 'png');


%% -----------------------------------------------------------------------------------------
% Setup for main figure generation.
%------------------------------------------------------------------------------------------
% threshold for full color saturation in SNP/LOH figure.
% synced to bases_per_bin as below, or defaulted to 50.
full_data_threshold = floor(bases_per_bin/100);

fig = figure(1);
set(gcf, 'Position', [0 70 1024 600]);
data_mode = 1;
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		if (data_mode == 1)
			for chr_bin = 1:length(chr_SNPdata{chr,1})
				% Regenerate chr plot data if the save file does not exist.
				SNPs_count{chr}(chr_bin)                                     = length(chr_SNPdata{chr,1}{chr_bin}) + length(chr_SNPdata{chr,2}{chr_bin});   % the number of data points in this bin.
				SNPs_to_fullData_ratio{chr}(chr_bin)                         = SNPs_count{chr}(chr_bin)/full_data_threshold;                                % divide by the threshold for full color saturation in SNP/LOH figure.
				SNPs_to_fullData_ratio{chr}(SNPs_to_fullData_ratio{chr} > 1) = 1;                                                                           % any bins with more data than the threshold for full color saturation around limited to full saturation.

				% darren
				%fprintf(['chr' num2str(chr) ':chr_bin' num2str(chr_bin) ' => copyNum' num2str(round()) '; data = (' num2str(chr_SNPdata{chr,1}{chr_bin}) '):(' num2str(chr_SNPdata{chr,2}{chr_bin}) ')\n']);
    
				HOMplot{chr}(chr_bin)                  = length(chr_SNPdata{chr,1}{chr_bin});         % phased data.
				HOMplot2{chr}(chr_bin)                 = HOMplot{chr}(chr_bin)/full_data_threshold;   %
				HOMplot2{chr}(HOMplot2{chr} > 1)       = 1;                                           %

				HETplot{chr}(chr_bin)                  = length(chr_SNPdata{chr,2}{chr_bin});         % unphased data.
				HETplot2{chr}(chr_bin)                 = HETplot{chr}(chr_bin)/full_data_threshold;   %
				HETplot2{chr}(HETplot2{chr} > 1)       = 1;                                           %
			end;
		end;
	end;
end;
fprintf('\n');
largestChr = find(chr_width == max(chr_width));

%% -----------------------------------------------------------------------------------------
% Setup for linear-view figure generation.
%-------------------------------------------------------------------------------------------
if (Linear_display == true)
	Linear_fig = figure(2);
	Linear_genome_size   = sum(chr_size);

	Linear_Chr_max_width = 0.91;               % width for all chromosomes across figure.  1.00 - leftMargin - rightMargin - subfigure gaps.
	Linear_left_start    = 0.01;               % left margin (also right margin).
	Linear_left_chr_gap  = 0.07/(num_chrs-1);  % gaps between chr subfigures.

	Linear_height        = 0.6;
	Linear_base          = 0.1;
	Linear_TickSize      = -0.01;  %negative for outside, percentage of longest chr figure.
	maxY                 = ploidyBase*2;
	Linear_left          = Linear_left_start;
end;


%% -----------------------------------------------------------------------------------------
% Make figures
%-------------------------------------------------------------------------------------------
first_chr = true;
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
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
	    for i = 1:length(SNPs_to_fullData_ratio{chr})+1;
	        if (i-1 < length(SNPs_to_fullData_ratio{chr}))
	            c_tot_post = SNPs_to_fullData_ratio{chr}(i)+SNPs_to_fullData_ratio{chr}(i);
	            if (c_tot_post == 0)
	                c_post = colorNoData;
	            else
					% darren
					colorMix = colorHET   *   HETplot2{chr}(i)/SNPs_to_fullData_ratio{chr}(i) + ...
					           colorHOM   *   HOMplot2{chr}(i)/SNPs_to_fullData_ratio{chr}(i);
					c_post =   colorMix   *   min(1,SNPs_to_fullData_ratio{chr}(i)) + ...
					           colorNoData*(1-min(1,SNPs_to_fullData_ratio{chr}(i)));
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
	    xlim([0,chr_size(chr)/bases_per_bin]);
    
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
	    set(gca,'YTick',[0 maxY/4 maxY/2 maxY/4*3 maxY]);
	    set(gca,'YTickLabel',{'','','','',''});
	    text(-50000/bases_per_bin, maxY/4,   '1','HorizontalAlignment','right','Fontsize',10);
	    text(-50000/bases_per_bin, maxY/2,   '2','HorizontalAlignment','right','Fontsize',10);
	    text(-50000/bases_per_bin, maxY/4*3, '3','HorizontalAlignment','right','Fontsize',10);
	    text(-50000/bases_per_bin, maxY,     '4','HorizontalAlignment','right','Fontsize',10);

	    set(gca,'FontSize',12);
	    if (chr == find(chr_posY == max(chr_posY)))
			title([ project ' vs. (hapmap)' hapmap ' SNP/LOH map'],'Interpreter','none','FontSize',24);
	    end;
	    hold on;
	    %end axes labels etc.
    
	    %show centromere outlines and horizontal marks.
	    x1 = cen_start(chr)/bases_per_bin;
	    x2 = cen_end(chr)/bases_per_bin;
	    leftEnd  = 0.5*5000/bases_per_bin;
	    rightEnd = (chr_size(chr) - 0.5*5000)/bases_per_bin;
	    if (Centromere_format == 0)
	        % standard chromosome cartoons in a way which will not cause segfaults when running via commandline.
	        dx = cen_tel_Xindent; %5*5000/bases_per_bin;
	        dy = cen_tel_Yindent; %maxY/10;
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
	        plot([leftEnd rightEnd], [-maxY/10*1.5 -maxY/10*1.5],'color',[0 0 0]);
	        hold on;
	        annotation_location = (annotation_start+annotation_end)./2;
	        for i = 1:length(annotation_location)
	            if (annotation_chr(i) == chr)
	                annotationloc = annotation_location(i)/bases_per_bin-0.5*(5000/bases_per_bin);
	                annotationStart = annotation_start(i)/bases_per_bin-0.5*(5000/bases_per_bin);
	                annotationEnd   = annotation_end(i)/bases_per_bin-0.5*(5000/bases_per_bin);
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

	    %% Linear figure draw section
	    if (Linear_display == true)
	        figure(Linear_fig);
	        Linear_width = Linear_Chr_max_width*chr_size(chr)/Linear_genome_size;
	        subplot('Position',[Linear_left Linear_base Linear_width Linear_height]);
	        Linear_left = Linear_left + Linear_width + Linear_left_chr_gap;
	        hold on;
	        title(chr_label{chr},'Interpreter','none','FontSize',20);

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

	        %show segmental anueploidy breakpoints.
	        if (displayBREAKS == true)
	            for segment = 2:length(chr_breaks{chr})-1
	                bP = chr_breaks{chr}(segment)*length(HETplot2{chr});
	                c_ = [0 0 1];
	                x_ = [bP bP bP-1 bP-1];
	                y_ = [0 maxY maxY 0];
	                f = fill(x_,y_,c_);   
	                set(f,'linestyle','none');
	            end;
	        end;

	        %show centromere.
	        x1 = cen_start(chr)/bases_per_bin;
	        x2 = cen_end(chr)/bases_per_bin;
	        leftEnd  = 0.5*5000/bases_per_bin;
	        rightEnd = (chr_size(chr) - 0.5*5000)/bases_per_bin;

	        if (Centromere_format == 0)
	            % standard chromosome cartoons in a way which will not cause segfaults when running via commandline.
	            dx = cen_tel_Xindent; %5*5000/bases_per_bin;
	            dy = cen_tel_Yindent; %maxY/10;
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
	            plot([leftEnd rightEnd], [-maxY/10*1.5 -maxY/10*1.5],'color',[0 0 0]);
	            hold on;
	            annotation_location = (annotation_start+annotation_end)./2;
	            for i = 1:length(annotation_location)
	                if (annotation_chr(i) == chr)
	                    annotationloc = annotation_location(i)/bases_per_bin-0.5*(5000/bases_per_bin);
	                    annotationStart = annotation_start(i)/bases_per_bin-0.5*(5000/bases_per_bin);
	                    annotationEnd   = annotation_end(i)/bases_per_bin-0.5*(5000/bases_per_bin);
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
	        xlim([0,chr_size(chr)/bases_per_bin]);
	        % modify y axis limits to show annotation locations if any are provided.
	        if (length(annotations) > 0)
	            ylim([-maxY/10*1.5,maxY]);
	        else
	            ylim([0,maxY]);
	        end;
	        set(gca,'TickLength',[(Linear_TickSize*chr_size(1)/chr_size(chr)) 0]); %ensures same tick size on all subfigs.
	        set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
	        set(gca,'XTickLabel',[]);
	        if (first_chr == true)
				ylabel({project;'vs. (hapmap)';hapmap}, 'Rotation', 0, 'HorizontalAlign', 'right', 'VerticalAlign', 'bottom','Interpreter','none','FontSize',10);
	            set(gca,'YTick',[0 maxY/4 maxY/2 maxY/4*3 maxY]);
	            set(gca,'YTickLabel',{'','','','',''});
		    text(-50000/bases_per_bin, maxY/4,   '1','HorizontalAlignment','right','Fontsize',10);
		    text(-50000/bases_per_bin, maxY/2,   '2','HorizontalAlignment','right','Fontsize',10);
		    text(-50000/bases_per_bin, maxY/4*3, '3','HorizontalAlignment','right','Fontsize',10);
		    text(-50000/bases_per_bin, maxY,     '4','HorizontalAlignment','right','Fontsize',10);
	        else
	            set(gca,'YTick',[]);
	            set(gca,'YTickLabel',[]);
	        end;
	        set(gca,'FontSize',12);
	        %end final reformatting.
	        
	        % shift back to main figure generation.
	        figure(fig);
	        hold on;
			first_chr = false;
	    end;
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

%% Save figures.
set(fig,'PaperPosition',[0 0 8 6]*2);
saveas(fig,        [projectDir 'fig.SNP-map.1.eps'], 'epsc');
saveas(fig,        [projectDir 'fig.SNP-map.1.png'], 'png');
set(Linear_fig,'PaperPosition',[0 0 8 0.62222222]*2);
saveas(Linear_fig, [projectDir 'fig.SNP-map.2.eps'], 'epsc');
saveas(Linear_fig, [projectDir 'fig.SNP-map.2.png'], 'png');

%% Delete figures from memory.
delete(fig);
delete(Linear_fig);

%% ========================================================================
% end stuff
%==========================================================================
end
