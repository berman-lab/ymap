function [] = LOH_hapmap_v4(main_dir,user,genomeUser,project,parent_or_hapmap,genome,ploidyEstimateString,ploidyBaseString, SNP_verString,LOH_verString,CNV_verString,displayBREAKS);
addpath('../');
workingDir = [main_dir 'users/' user '/projects/' project '/'];
fprintf('\n\n\t*===============================================================*\n');
fprintf(    '\t| Generate SNP/LOH only plot in script "LOH_hapmap_v4.m".       |\n');
fprintf(    '\t*---------------------------------------------------------------*\n');
tic;


%% ========================================================================
%    Centromere_format          : Controls how centromeres are depicted.   [0..2]   '2' is pinched cartoon default.
%    bases_per_bin              : Controls bin sizes for SNP/CGH fractions of plot.
%    scale_type                 : 'Ratio' or 'Log2Ratio' y-axis scaling of copy number.
%                                 'Log2Ratio' does not properly scale CGH data by ploidy.
%    Chr_max_width              : max width of chrs as fraction of figure width.
fprintf('\t|\tSetup for processing.\n');
Centromere_format              = 0;
Chr_max_width                  = 0.8;
colorBars                      = true;
blendColorBars                 = false;
show_annotations               = true;
Yscale_nearest_even_ploidy     = true;
Linear_display                 = true;
Linear_displayBREAKS           = false;
AnglePlot                      = true;   % Show histogram of alleleic fraction at the left end of standard figure chromosomes.
FillColors                     = true;   %     Fill histogram using colors.
show_uncalibrated              = false;  %     Fill with single color instead of ratio call colors.



%%=========================================================================
% Load FASTA file name from 'reference.txt' file for project.
%--------------------------------------------------------------------------
fprintf('\t|\tLoad FASTA file name for the genome in use.\n');
userReference                  = [main_dir 'users/' user '/genomes/' genome '/reference.txt'];
defaultReference               = [main_dir 'users/default/genomes/' genome '/reference.txt'];
if (exist(userReference,'file') == 0)   
	FASTA_string               = strtrim(fileread(defaultReference));
else                    
	FASTA_string               = strtrim(fileread(userReference));
end;
[FastaPath,FastaName,FastaExt] = fileparts(FASTA_string);


%%=========================================================================
% Control variables for Candida albicans SC5314.
%--------------------------------------------------------------------------
projectDir = [main_dir 'users/' user '/projects/' project '/'];
genomeDir  = [main_dir 'users/' genomeUser '/genomes/' genome '/'];


fprintf('\t|\tDetermine if hapmap is in use.\n');
if (exist([main_dir 'users/default/hapmaps/' parent_or_hapmap '/'], 'dir') == 7)
	useHapmap  = true;
	hapmapDir  = [main_dir 'users/default/hapmaps/' parent_or_hapmap '/'];   % system hapmap.
	hapmap     = parent_or_hapmap;
	hapmapUser = 'default';
elseif (exist([main_dir 'users/' user '/hapmaps/' parent_or_hapmap '/'], 'dir') == 7)
	useHapmap  = true;
	hapmapDir  = [main_dir 'users/' user '/hapmaps/' parent_or_hapmap '/'];  % user hapmap.
	hapmap     = parent_or_hapmap;
	hapmapUser = user;
else
	useHapmap  = false;
	hapmapDir  = '';
	hapmap     = '';
	hapmapUser = '';
end;


fprintf('\t|\tDetermine if parent project is in use.\n');
% The 'parent' will == the 'project' when no 'parent' is selected in setup.
if (strcmp(project,parent_or_hapmap) == 0)   % different
	useParent  = true;
	if (exist([main_dir 'users/default/projects/' parent_or_hapmap '/'], 'dir') == 7)
		parentDir  = [main_dir 'users/default/projects/' parent_or_hapmap '/'];   % system parent.
		parentUser = 'default';
	else
		parentDir  = [main_dir 'users/' user '/projects/' parent_or_hapmap '/'];  % user parent.
		parentUser = user;
	end;
	parent     = parent_or_hapmap;
else
	useParent  = false;
	parentDir  = projectDir;
	parent     = project;
	parentUser = user;
end;


fprintf('\t|\tLoad deteilas of genome in use.\n');
[centromeres, chr_sizes, figure_details, annotations, ploidy_default] = Load_genome_information(genomeDir);
[Aneuploidy]                                                          = Load_dataset_information(projectDir);
%% check if dataset uses original chromosomes names
originalNamePath = [projectDir 'original.txt'];
if (exist(originalNamePath,'file'))
	useOriginal = true;
else 
	useOriginal = false;
end;

num_chrs                          = length(chr_sizes);
for i = 1:length(chr_sizes)
	chr_size(i)                   = 0;
	cen_start(i)                  = 0;
	cen_end(i)                    = 0;
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
		if (useOriginal && length(figure_details(i).name) < 10)
		    chr_label {figure_details(i).chr} = figure_details(i).name;
		else
		   chr_label {figure_details(i).chr} = figure_details(i).label;
		end;
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
bases_per_bin = max(chr_size)/700;


%% =========================================================================================
%% =========================================================================================
%% =========================================================================================
%% -----------------------------------------------------------------------------------------
%% =========================================================================================
%% =========================================================================================
%% =========================================================================================


fprintf('\t|\tProcess input ploidy.\n');
% Process input ploidy.
ploidy = str2num(ploidyEstimateString);
% Sanitize user input of euploid state.
ploidyBase = round(str2num(ploidyBaseString));
if (ploidyBase > 4);   ploidyBase = 4;   end;
if (ploidyBase < 1);   ploidyBase = 1;   end;
fprintf(['\nEuploid base = "' num2str(ploidyBase) '"\n']);


%% =========================================================================================
% Define colors for figure generation.
%-------------------------------------------------------------------------------------------
fprintf('\t|\tLoad color definitions.\n');
phased_and_unphased_color_definitions;


%%================================================================================================
% Setup for SNP/LOH data calculations.
%-------------------------------------------------------------------------------------------------
fprintf('\t|\tInitialize data vectors for tracking data presentation.\n');
% Initializes vectors used to hold allelic ratios for each chromosome segment.
for chr = 1:length(chr_sizes)
	% Build data structure for SNP information:  chr_SNPdata{chr,j}{chr_bin} = [];
	%       1 : phased SNP ratio data.
	%       2 : unphased SNP ratio data.
	%       3 : phased SNP position data.
	%       4 : unphased SNP position data.
	%       5 : phased SNP allele strings.   (baseCall:alleleA/alleleB)
	%       6 : unphased SNP allele strings.
	chr_length = ceil(chr_size(chr)/bases_per_bin);
	for j = 1:6
		chr_SNPdata{chr,j} = cell(1,chr_length);
	end;
	% Setup to track RGB values used to present SNP/LOH data for each chromosome bin.
	for j = 1:3
		% Track the RGB value sum per standard bin, then divide by the count to reach the average color per standard genome bin.
		chr_SNPdata_colorsC{chr,j}           = zeros(chr_length,1);
		chr_SNPdata_colorsP{chr,j}           = zeros(chr_length,1);
	end;
	% Track the number of SNP colors per standard bin.
	chr_SNPdata_countC{chr} = zeros(chr_length,1);
	chr_SNPdata_countP{chr} = zeros(chr_length,1);
end;


%%================================================================================================
% Load 'Common_CNV.mat' file containing CNV estimates per standard genome bin.
%-------------------------------------------------------------------------------------------------
fprintf('\t|\tLoading "Common_CNV" data file, to be used in copy number estimation.\n');
load([projectDir 'Common_CNV.mat']);   % 'CNVplot2', 'genome_CNV'
[chr_breaks, chrCopyNum, ploidyAdjust] = FindChrSizes_4(Aneuploidy,CNVplot2,ploidy,num_chrs,chr_in_use);


%% =========================================================================================
% Test adjacent segments for no change in copy number estimate.
%...........................................................................................
% Adjacent pairs of segments with the same copy number will be fused into a single segment.
% Segments with a <= zero copy number will be fused to an adjacent segment.
%-------------------------------------------------------------------------------------------
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		if (length(chrCopyNum{chr}) > 1)  % more than one segment, so lets examine if adjacent segments have different copyNums.
			%% Merge any adjacent segments with the same copy number.
			% add break representing left end of chromosome.
			breakCount_new         = 1;
			chr_breaks_new{chr}    = [];
			chrCopyNum_new{chr}    = [];
			chr_breaks_new{chr}(1) = 0.0;
			fprintf(['\nlength(chrCopyNum{chr}) = ' num2str(length(chrCopyNum{chr})) '\n']);
			if (length(chrCopyNum{chr}) > 0)
				fprintf(['chrCopyNum{chr}(1) = ' num2str(chrCopyNum{chr}(1)) '\n']);
				chrCopyNum_new{chr}(1) = chrCopyNum{chr}(1);
				for segment = 1:(length(chrCopyNum{chr})-1)
					if (round(chrCopyNum{chr}(segment)) == round(chrCopyNum{chr}(segment+1)))
						% two adjacent segments have identical copyNum and should be fused into one; don't add boundry to new list.
					else
						% two adjacent segments have different copyNum; add boundry to new list.
						breakCount_new                      = breakCount_new + 1;
						chr_breaks_new{chr}(breakCount_new) = chr_breaks{chr}(segment+1);
						chrCopyNum_new{chr}(breakCount_new) = chrCopyNum{chr}(segment+1);
					end;
				end;
			end;
			% add break representing right end of chromosome.
			breakCount_new = breakCount_new+1;
			chr_breaks_new{chr}(breakCount_new) = 1.0;
			fprintf(['@@@ chr = ' num2str(chr) '\n']);
			fprintf(['@@@    chr_breaks_old = ' num2str(chr_breaks{chr})     '\n']);
			fprintf(['@@@    chrCopyNum_old = ' num2str(chrCopyNum{chr})     '\n']);
			fprintf(['@@@    chr_breaks_new = ' num2str(chr_breaks_new{chr}) '\n']);
			fprintf(['@@@    chrCopyNum_new = ' num2str(chrCopyNum_new{chr}) '\n']);
			% copy new lists to old.
			chr_breaks{chr} = chr_breaks_new{chr};
			chrCopyNum{chr} = [];
			chrCopyNum{chr} = chrCopyNum_new{chr};
		end;
	end;
end;


%%================================================================================================
% Load SNP/LOH data.
%.................................................................................................
% if (useHapmap)
%       chr_SNPdata{chr,1}{chr_bin} = phased SNP ratio data.
%       chr_SNPdata{chr,2}{chr_bin} = unphased SNP ratio data.
%       chr_SNPdata{chr,3}{chr_bin} = phased SNP position data.
%       chr_SNPdata{chr,4}{chr_bin} = unphased SNP position data.
%       chr_SNPdata{chr,5}{chr_bin} = flipper value for phased SNP.
%       chr_SNPdata{chr,6}{chr_bin} = flipper value for unphased SNP.
% elseif (useParent)
%       chr_SNPdata{chr,1}{chr_bin} = parent SNP ratio data.
%       chr_SNPdata{chr,2}{chr_bin} = child SNP ratio data.
%       chr_SNPdata{chr,3}{chr_bin} = parent SNP position data.
%       chr_SNPdata{chr,4}{chr_bin} = child SNP position data.
% else
%       chr_SNPdata{chr,2}{chr_bin} = child SNP ratio data.
%       chr_SNPdata{chr,4}{chr_bin} = child SNP position data.
% end;
%-------------------------------------------------------------------------------------------------
%	if (useHapmap)
%		phased_ratio_data_string    = '(1.0)';
%		unphased_ratio_data_string  = '()';
%		phased_coordinates_string   = '(19846)';
%		unphased_coordinates_string = '()';
%		phased_alleles_string       = '(G:G/C)';
%		unphased_alleles_string     = '()';
%	elseif (useParent)
%		% coordinates that are heterozygous (0.25 .. 0.75) in the parent.
%		phased_ratio_data_string    = '()';
%		unphased_ratio_data_string  = '(1.0,1.0)';
%		phased_coordinates_string   = '()';
%		unphased_coordinates_string = '(19846,20572)';
%		phased_alleles_string       = '()';
%		unphased_alleles_string     = '(Z:G/C,Z:A/C)';
%	else
%		% coordinate sthat are heterozygous (0.25 .. 0.75) in the experimental dataset.
%		phased_ratio_data_string    = '()';
%		unphased_ratio_data_string  = '(0.529411764706,0.538461538462)';
%		phased_coordinates_string   = '()';
%		unphased_coordinates_string = '(19846,20572)';
%		phased_alleles_string       = '()';
%		unphased_alleles_string     = '(Z:G/C,Z:A/C)';
%	end;

fprintf('\t|\tLoad SNP/LOH data.\n');
if (exist([projectDir 'SNP_' SNP_verString '.mat'],'file') == 0)
	fprintf(['\t|\t\tMAT file "SNP_' SNP_verString '.mat" not found, regenerating.']);
	datafile       = [projectDir 'preprocessed_SNPs.txt'];
	data           = fopen(datafile, 'r');
	count          = 0;
	old_chr        = 0;
	gap_string     = '';	
	while not (feof(data))
		dataLine = fgetl(data);
		if (length(dataLine) > 0)
			if (dataLine(1) ~= '#')
				% process the loaded line into data channels.
				lineVariables = textscan(dataLine, '%f %f %f %s %s %s %s %s %s');
				chr_num = lineVariables{1};
				fragment_start = lineVariables{2};
				fragment_end = lineVariables{3};
				phased_ratio_data_string = lineVariables{4}{1};
				unphased_ratio_data_string = lineVariables{5}{1};
				phased_coordinates_string = lineVariables{6}{1};
				unphased_coordinates_string = lineVariables{7}{1};
				phased_alleles_string = lineVariables{8}{1};
				unphased_alleles_string = lineVariables{9}{1};
				
				
				% format = simple, one number per column.
				chr_length                  = ceil(chr_size(chr_num)/bases_per_bin);
				chr_bin                     = ceil(fragment_start/bases_per_bin);

				% Log file output to indicate progression of this section of code.
				count = count+1;
				if (old_chr ~= chr_num)
					fprintf(['\n\t|\t\t' chr_name{chr_num} '\n\t|\t' gap_string]);
				end;
				if (mod(count,10) == 0)
					fprintf('.');
					gap_string = [gap_string ' '];
				end;
				if (count == 800)
					fprintf('\n\t|\t\t');
					count = 0;
					gap_string = '';
				end;
				old_chr = chr_num;

				% format = '(number1,number2,...,numberN)'
				phased_ratio_data_string(1)              = [];
				phased_ratio_data_string(end)            = [];
				if (length(phased_ratio_data_string)    == 0)
					phased_ratio_data                = [];
				else
					commaCount                       = length(find(phased_ratio_data_string==','));
					if (commaCount == 0)
						phased_ratio_data        = str2num(phased_ratio_data_string);
					else
						phased_ratio_data        = strsplit(phased_ratio_data_string,',');   % function converts number lists from strings to numbers.
					end;
				end;

				% format = '(number1,number2,...,numberN)'
				phased_coordinates_string(1)             = [];
				phased_coordinates_string(end)           = [];
				if (length(phased_coordinates_string)   == 0)
					phased_coordinates               = [];
				else
					commaCount                       = length(find(phased_coordinates_string==','));
					if (commaCount == 0)
						phased_coordinates       = str2num(phased_coordinates_string);
					else
						phased_coordinates       = strsplit(phased_coordinates_string,',');   % function converts number lists from strings to numbers.
					end;
				end;

				% format = '(number1,number2,...,numberN)'
				unphased_ratio_data_string(1)            = [];
				unphased_ratio_data_string(end)          = [];
				if (length(unphased_ratio_data_string)  == 0)
					unphased_ratio_data              = [];
				else
					commaCount                       = length(find(unphased_ratio_data_string==','));
					if (commaCount == 0)
						unphased_ratio_data      = str2num(unphased_ratio_data_string);
					else
						unphased_ratio_data      = strsplit(unphased_ratio_data_string,',');  % function converts number lists from strings to numbers.
					end;
				end;

				% format = '(number1,number2,...,numberN)'
				unphased_coordinates_string(1)           = [];
				unphased_coordinates_string(end)         = [];
				if (length(unphased_coordinates_string) == 0)
					unphased_coordinates             = [];
				else
					commaCount                       = length(find(unphased_coordinates_string==','));
					if (commaCount == 0)
						unphased_coordinates     = str2num(unphased_coordinates_string);
					else
						unphased_coordinates     = strsplit(unphased_coordinates_string,',');   % function converts number lists from strings to numbers.
					end;
				end;

				% format = '(A:A/T,C:G/C,...,T:A/T)'
				phased_alleles_string(1)                 = [];
				phased_alleles_string(end)               = [];
				if (length(phased_alleles_string)       == 0)
					phased_alleles                   = [];
				else
					commaCount                       = length(find(phased_alleles_string==','));
					if (commaCount == 0)
						phased_alleles           = phased_alleles_string;
					else
						phased_alleles           = strsplit(phased_alleles_string,',');   % function converts number lists from strings to numbers.
					end;
				end;

				% format = '(A:A/T,C:G/C,...,T:A/T)'
				unphased_alleles_string(1)               = [];
				unphased_alleles_string(end)             = [];
				if (length(unphased_alleles_string)     == 0)
					unphased_alleles                 = [];
				else
					commaCount                       = length(find(unphased_alleles_string==','));
					if (commaCount == 0)
						unphased_alleles         = unphased_alleles_string;
					else
						unphased_alleles         = strsplit(unphased_alleles_string,',');   % function converts number lists from strings to numbers.
					end;
				end;

				% add phased and unphased data to storage arrays.
				chr_SNPdata{chr_num,1}{chr_bin}          = phased_ratio_data;
				chr_SNPdata{chr_num,2}{chr_bin}          = unphased_ratio_data;

				% add phased and unphased data coordinates to storage arrays.
				chr_SNPdata{chr_num,3}{chr_bin}          = phased_coordinates;
				chr_SNPdata{chr_num,4}{chr_bin}          = unphased_coordinates;

				% add phased and unphased data allele strings to storage arrays.
				chr_SNPdata{chr_num,5}{chr_bin}          = phased_alleles;
				chr_SNPdata{chr_num,6}{chr_bin}          = unphased_alleles;
			end;
		end;
	end;
	fclose(data);

	save([projectDir 'SNP_' SNP_verString '.mat'],'chr_SNPdata');
else
	fprintf('\t|\t\tMAT file found, loading.\n');
	load([projectDir 'SNP_' SNP_verString '.mat']);
end;


%%================================================================================================
% Process SNP/hapmap data to determine colors to be presented for each SNP locus.
%-------------------------------------------------------------------------------------------------
%    Calculate allelic fraction cutoffs for each segment and populate data structure containing
%    SNP phasing information.
%        chr_SNPdata{chr,1}{chr_bin} = phased SNP ratio data.
%        chr_SNPdata{chr,2}{chr_bin} = unphased SNP ratio data.
%        chr_SNPdata{chr,3}{chr_bin} = phased SNP position data.
%        chr_SNPdata{chr,4}{chr_bin} = unphased SNP position data.
%        chr_SNPdata{chr,5}{chr_bin} = phased SNP allele strings.   (baseCall:alleleA/alleleB)
%        chr_SNPdata{chr,6}{chr_bin} = unphased SNP allele strings.
%-------------------------------------------------------------------------------------------

fprintf('\t|\tCalculate allelic ratio cutoffs using Gaussian fitting.\n');
temp_holding = chr_SNPdata;
calculate_allelic_ratio_cutoffs;
chr_SNPdata = temp_holding;


%% =========================================================================================
% Define new colors for SNPs, using Gaussian fitting crossover points as ratio cutoffs.
%-------------------------------------------------------------------------------------------
fprintf('\t|\tDetermine display color for each SNP.\n');
for chr = 1:num_chrs
	% avoid running over chromosomes with empty copy number
	if (chr_in_use(chr) == 1 && ~isempty(chrCopyNum{chr}))
		for chr_bin = 1:ceil(chr_size(chr)/bases_per_bin)
			%
			% Determining colors for each SNP coordinate from calculated cutoffs.
			%
			localCopyEstimate                       = round(CNVplot2{chr}(chr_bin)*ploidy*ploidyAdjust);
			allelic_ratios                          = [chr_SNPdata{chr,1}{chr_bin} chr_SNPdata{chr,2}{chr_bin}];
			coordinates                             = [chr_SNPdata{chr,3}{chr_bin} chr_SNPdata{chr,4}{chr_bin}];
			if (length(chr_SNPdata{chr,1}{chr_bin}) == 1) && (length(chr_SNPdata{chr,2}{chr_bin}) == 1)
				allele_strings                  = {chr_SNPdata{chr,5}{chr_bin} chr_SNPdata{chr,6}{chr_bin}};
			else
				allele_strings                  = [chr_SNPdata{chr,5}{chr_bin} chr_SNPdata{chr,6}{chr_bin}];
			end;

			if (length(allelic_ratios) > 0)
				for SNP = 1:length(allelic_ratios)
					% Load phased SNP data from earlier defined structure.
					allelic_ratio                         = allelic_ratios(SNP);
					coordinate                            = coordinates(SNP);
					if (length(allelic_ratios) > 1)
						allele_string                 = allele_strings{SNP};
					else
						allele_string                 = allele_strings;
					end;
					baseCall                              = allele_string(1);
					homologA                              = allele_string(3);
					homologB                              = allele_string(5);					

					% identify the segment containing the SNP.
					segmentID                             = 0;
					for segment = 1:(length(chrCopyNum{chr}))
						segment_start                 = chr_breaks{chr}(segment  )*chr_size(chr);
						segment_end                   = chr_breaks{chr}(segment+1)*chr_size(chr);
						if (coordinate > segment_start) && (coordinate <= segment_end)
							segmentID             = segment;
						end;
					end;

					% Load cutoffs between Gaussian fits performed earlier.
					segment_copyNum                       = round(chrCopyNum{              chr}(segmentID));
					actual_cutoffs                        = chrSegment_actual_cutoffs{     chr}{segmentID};
					mostLikelyGaussians                   = chrSegment_mostLikelyGaussians{chr}{segmentID};

					% Calculate allelic ratio on range of [1..200].
					SNPratio_int                          = (allelic_ratio)*199+1;

					% Identify the allelic ratio region containing the SNP.
					cutoffs                               = [1 actual_cutoffs 200];
					ratioRegionID                         = 0;
					for GaussianRegionID = 1:length(mostLikelyGaussians)
						cutoff_start                  = cutoffs(GaussianRegionID  );
						cutoff_end                    = cutoffs(GaussianRegionID+1);
						if (GaussianRegionID == 1)
							if (SNPratio_int >= cutoff_start) && (SNPratio_int <= cutoff_end)
								ratioRegionID = mostLikelyGaussians(GaussianRegionID);
							end;
						else
							if (SNPratio_int > cutoff_start) && (SNPratio_int <= cutoff_end)
								ratioRegionID = mostLikelyGaussians(GaussianRegionID);
							end;
						end;
					end;

					if (segment_copyNum <= 0);                          colorList = colorNoData;
					elseif (segment_copyNum == 1)
						% allelic fraction cutoffs: [0.50000] => [A B]
						if ((baseCall == homologA) || (baseCall == homologB))
							if (ratioRegionID == 2);            colorList = colorB;
							else                                colorList = colorA;
							end;
						else
							if (useParent)                      colorList = unphased_color_1of1;
							else                                colorList = colorNoData;
							end;
						end;
					elseif (segment_copyNum == 2)
						%   allelic fraction cutoffs: [0.25000 0.75000] => [AA AB BB]
						if ((baseCall == homologA) || (baseCall == homologB))
							if (ratioRegionID == 3);            colorList = colorBB;
							elseif (ratioRegionID == 2);        colorList = colorAB;
							else                                colorList = colorAA;
							end;
						else
							if (ratioRegionID == 3);            colorList = unphased_color_2of2;
							elseif (ratioRegionID == 2);        colorList = unphased_color_1of2;
							else                                colorList = unphased_color_2of2;
							end;
						end;
					elseif (segment_copyNum == 3)
						% allelic fraction cutoffs: [0.16667 0.50000 0.83333] => [AAA AAB ABB BBB]
						if ((baseCall == homologA) || (baseCall == homologB))
							if (ratioRegionID == 4);            colorList = colorBBB;
							elseif (ratioRegionID == 3);        colorList = colorABB;
							elseif (ratioRegionID == 2);        colorList = colorAAB;
							else                                colorList = colorAAA;
							end;
						else
							if (ratioRegionID == 4);            colorList = unphased_color_3of3;
							elseif (ratioRegionID == 3);        colorList = unphased_color_2of3;
							elseif (ratioRegionID == 2);        colorList = unphased_color_2of3;
							else                                colorList = unphased_color_3of3;
							end;
						end;
					elseif (segment_copyNum == 4)
						% allelic fraction cutoffs: [0.12500 0.37500 0.62500 0.87500] => [AAAA AAAB AABB ABBB BBBB]
						if ((baseCall == homologA) || (baseCall == homologB))
							if (ratioRegionID == 5);            colorList = colorBBBB;
							elseif (ratioRegionID == 4);        colorList = colorABBB;
							elseif (ratioRegionID == 3);        colorList = colorAABB;
							elseif (ratioRegionID == 2);        colorList = colorAAAB;
							else                                colorList = colorAAAA;
							end;
						else
							if (ratioRegionID == 5);            colorList = unphased_color_4of4;
							elseif (ratioRegionID == 4);        colorList = unphased_color_3of4;
							elseif (ratioRegionID == 3);        colorList = unphased_color_2of4;
							elseif (ratioRegionID == 2);        colorList = unphased_color_3of4;
							else                                colorList = unphased_color_4of4;
							end;
						end;
					elseif (segment_copyNum == 5)
						% allelic fraction cutoffs: [0.10000 0.30000 0.50000 0.70000 0.90000] => [AAAAA AAAAB AAABB AABBB ABBBB BBBBB]
						if ((baseCall == homologA) || (baseCall == homologB))
							if (ratioRegionID == 6);            colorList = colorBBBBB;
							elseif (ratioRegionID == 5);        colorList = colorABBBB;
							elseif (ratioRegionID == 4);        colorList = colorAABBB;
							elseif (ratioRegionID == 3);        colorList = colorAAABB;
							elseif (ratioRegionID == 2);        colorList = colorAAAAB;
							else                                colorList = colorAAAAA;
							end;
						else
							if (ratioRegionID == 6);            colorList = unphased_color_5of5;
							elseif (ratioRegionID == 5);        colorList = unphased_color_4of5;
							elseif (ratioRegionID == 4);        colorList = unphased_color_3of5;
							elseif (ratioRegionID == 3);        colorList = unphased_color_3of5;
							elseif (ratioRegionID == 2);        colorList = unphased_color_4of5;
							else                                colorList = unphased_color_5of5;
							end;
						end;
					elseif (segment_copyNum == 6)
						% allelic fraction cutoffs: [0.08333 0.25000 0.41667 0.58333 0.75000 0.91667] => [AAAAAA AAAAAB AAAABB AAABBB AABBBB ABBBBB BBBBBB]
						if ((baseCall == homologA) || (baseCall == homologB))
							if (ratioRegionID == 7);            colorList = colorBBBBBB;
							elseif (ratioRegionID == 6);        colorList = colorABBBBB;
							elseif (ratioRegionID == 5);        colorList = colorAABBBB;
							elseif (ratioRegionID == 4);        colorList = colorAAABBB;
							elseif (ratioRegionID == 3);        colorList = colorAAAABB;
							elseif (ratioRegionID == 2);        colorList = colorAAAAAB;
							else                                colorList = colorAAAAAA;
							end;
						else
							if (ratioRegionID == 7);            colorList = unphased_color_6of6;
							elseif (ratioRegionID == 6);        colorList = unphased_color_5of6;
							elseif (ratioRegionID == 5);        colorList = unphased_color_4of6;
							elseif (ratioRegionID == 4);        colorList = unphased_color_3of6;
							elseif (ratioRegionID == 3);        colorList = unphased_color_4of6;
							elseif (ratioRegionID == 2);        colorList = unphased_color_5of6;
							else                                colorList = unphased_color_6of6;
							end;
						end;
					elseif (segment_copyNum == 7)
						% allelic fraction cutoffs: [0.07143 0.21429 0.35714 0.50000 0.64286 0.78571 0.92857] => [AAAAAAA AAAAAAB AAAAABB AAAABBB AAABBBB AABBBBB ABBBBBB BBBBBBB]
						if ((baseCall == homologA) || (baseCall == homologB))
							if (ratioRegionID == 8);            colorList = colorBBBBBBB;
							elseif (ratioRegionID == 7);        colorList = colorABBBBBB;
							elseif (ratioRegionID == 6);        colorList = colorAABBBBB;
							elseif (ratioRegionID == 5);        colorList = colorAAABBBB;
							elseif (ratioRegionID == 4);        colorList = colorAAAABBB;
							elseif (ratioRegionID == 3);        colorList = colorAAAAABB;
							elseif (ratioRegionID == 2);        colorList = colorAAAAAAB;
							else                                colorList = colorAAAAAAA;
							end;
						else
							if (ratioRegionID == 8);            colorList = unphased_color_7of7;
							elseif (ratioRegionID == 7);        colorList = unphased_color_6of7;
							elseif (ratioRegionID == 6);        colorList = unphased_color_5of7;
							elseif (ratioRegionID == 5);        colorList = unphased_color_4of7;
							elseif (ratioRegionID == 3);        colorList = unphased_color_4of7;
							elseif (ratioRegionID == 3);        colorList = unphased_color_5of7;
							elseif (ratioRegionID == 2);        colorList = unphased_color_6of7;
							else                                colorList = unphased_color_7of7;
							end;
						end;
					elseif (segment_copyNum == 8)
						% allelic fraction cutoffs: [0.06250 0.18750 0.31250 0.43750 0.56250 0.68750 0.81250 0.93750] => [AAAAAAAA AAAAAAAB AAAAAABB AAAAABBB AAAABBBB AAABBBBB AABBBBBB ABBBBBBB BBBBBBBB]
						if ((baseCall == homologA) || (baseCall == homologB))
							if (ratioRegionID == 9);            colorList = colorBBBBBBBB;
							elseif (ratioRegionID == 8);        colorList = colorABBBBBBB;
							elseif (ratioRegionID == 7);        colorList = colorAABBBBBB;
							elseif (ratioRegionID == 6);        colorList = colorAAABBBBB;
							elseif (ratioRegionID == 5);        colorList = colorAAAABBBB;
							elseif (ratioRegionID == 4);        colorList = colorAAAAABBB;
							elseif (ratioRegionID == 3);        colorList = colorAAAAAABB;
							elseif (ratioRegionID == 2);        colorList = colorAAAAAAAB;
							else                                colorList = colorAAAAAAAA;
							end;
						else
							if (ratioRegionID == 9);            colorList = unphased_color_8of8;
							elseif (ratioRegionID == 8);        colorList = unphased_color_7of8;
							elseif (ratioRegionID == 7);        colorList = unphased_color_6of8;
							elseif (ratioRegionID == 6);        colorList = unphased_color_5of8;
							elseif (ratioRegionID == 5);        colorList = unphased_color_4of8;
							elseif (ratioRegionID == 4);        colorList = unphased_color_5of8;
							elseif (ratioRegionID == 3);        colorList = unphased_color_6of8;
							elseif (ratioRegionID == 2);        colorList = unphased_color_7of8;
							else                                colorList = unphased_color_8of8;
							end;
						end;
					elseif (segment_copyNum >= 9)
						% allelic fraction cutoffs: [0.05556 0.16667 0.27778 0.38889 0.50000 0.61111 0.72222 0.83333 0.94444] => [AAAAAAAAA AAAAAAAAB AAAAAAABB AAAAAABBB AAAAABBBB AAAABBBBB AAABBBBBB AABBBBBBB
						%                                                                                                         ABBBBBBBB BBBBBBBBB]
						if ((baseCall == homologA) || (baseCall == homologB))
							if (ratioRegionID == 10);           colorList = colorBBBBBBBBB;
							elseif (ratioRegionID == 9);        colorList = colorABBBBBBBB;
							elseif (ratioRegionID == 8);        colorList = colorAABBBBBBB;
							elseif (ratioRegionID == 7);        colorList = colorAAABBBBBB;
							elseif (ratioRegionID == 6);        colorList = colorAAAABBBBB;
							elseif (ratioRegionID == 5);        colorList = colorAAAAABBBB;
							elseif (ratioRegionID == 4);        colorList = colorAAAAAABBB;
							elseif (ratioRegionID == 3);        colorList = colorAAAAAAABB;
							elseif (ratioRegionID == 2);        colorList = colorAAAAAAAAB;
							else                                colorList = colorAAAAAAAAA;
							end;
						else
							if (ratioRegionID == 10);           colorList = unphased_color_9of9;
							elseif (ratioRegionID == 9);        colorList = unphased_color_8of9;
							elseif (ratioRegionID == 8);        colorList = unphased_color_7of9;
							elseif (ratioRegionID == 7);        colorList = unphased_color_6of9;
							elseif (ratioRegionID == 6);        colorList = unphased_color_5of9;
							elseif (ratioRegionID == 5);        colorList = unphased_color_5of9;
							elseif (ratioRegionID == 4);        colorList = unphased_color_6of9;
							elseif (ratioRegionID == 3);        colorList = unphased_color_7of9;
							elseif (ratioRegionID == 2);        colorList = unphased_color_8of9;
							else                                colorList = unphased_color_9of9;
							end;
						end;
					end;
					chr_SNPdata_colorsC{chr,1}(chr_bin) = chr_SNPdata_colorsC{chr,1}(chr_bin) + colorList(1);
					chr_SNPdata_colorsC{chr,2}(chr_bin) = chr_SNPdata_colorsC{chr,2}(chr_bin) + colorList(2);
					chr_SNPdata_colorsC{chr,3}(chr_bin) = chr_SNPdata_colorsC{chr,3}(chr_bin) + colorList(3);
					chr_SNPdata_countC{ chr  }(chr_bin) = chr_SNPdata_countC{ chr  }(chr_bin) + 1;

					% Troubleshooting output.
					% fprintf(['chr = ' num2str(chr) '; seg = ' num2str(segment) '; bin = ' num2str(chr_bin) '; ratioRegionID = ' num2str(ratioRegionID) '\n']);
				end;
			end;
		end;

		%
		% Average colors of SNPs found in bin.
		%
		fprintf('\t|\tDetermine average color for SNPs in chromosome bin.\n');
		for chr_bin = 1:ceil(chr_size(chr)/bases_per_bin)
			allelic_ratios                                      = [chr_SNPdata{chr,1}{chr_bin} chr_SNPdata{chr,2}{chr_bin}];
			if (length(allelic_ratios) > 0)
				if (chr_SNPdata_countC{chr}(chr_bin) > 0)
					chr_SNPdata_colorsC{chr,1}(chr_bin) = chr_SNPdata_colorsC{chr,1}(chr_bin)/chr_SNPdata_countC{chr}(chr_bin);
					chr_SNPdata_colorsC{chr,2}(chr_bin) = chr_SNPdata_colorsC{chr,2}(chr_bin)/chr_SNPdata_countC{chr}(chr_bin);
					chr_SNPdata_colorsC{chr,3}(chr_bin) = chr_SNPdata_colorsC{chr,3}(chr_bin)/chr_SNPdata_countC{chr}(chr_bin);
				else
					chr_SNPdata_colorsC{chr,1}(chr_bin) = 1.0;
					chr_SNPdata_colorsC{chr,2}(chr_bin) = 1.0;
					chr_SNPdata_colorsC{chr,3}(chr_bin) = 1.0;
				end;
			else
				chr_SNPdata_colorsC{chr,1}(chr_bin)         = 1.0;
				chr_SNPdata_colorsC{chr,2}(chr_bin)         = 1.0;
				chr_SNPdata_colorsC{chr,3}(chr_bin)         = 1.0;
			end;
		end;
	end;
end;


%%================================================================================================
% Setup for main figure generation.
%-------------------------------------------------------------------------------------------------
% threshold for full color saturation in SNP/LOH figure.
% synced to bases_per_bin as below, or defaulted to 50.
fprintf('\t|\tCalculate color intensity for each chromosome bin.\n');
largestChr          = find(chr_width == max(chr_width));
largestChr          = largestChr(1);
full_data_threshold = floor(bases_per_bin/100);
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		for chr_bin = 1:ceil(chr_size(chr)/bases_per_bin)
			% the number of heterozygous data points in this bin.
			SNPs_count{chr}(chr_bin)                                     = length(chr_SNPdata{chr,1}{chr_bin}) + length(chr_SNPdata{chr,2}{chr_bin});

			% divide by the threshold for full color saturation in SNP/LOH figure.
			SNPs_to_fullData_ratio{chr}(chr_bin)                         = SNPs_count{chr}(chr_bin)/full_data_threshold;

			% any bins with more data than the threshold for full color saturation are limited to full saturation.
			SNPs_to_fullData_ratio{chr}(SNPs_to_fullData_ratio{chr} > 1) = 1;

			phased_plot{ chr}(chr_bin)                                   = length(chr_SNPdata{chr,1}{chr_bin});             % phased data.
			phased_plot2{chr}(chr_bin)                                   = phased_plot{chr}(chr_bin)/full_data_threshold;   %
			phased_plot2{chr}(phased_plot2{chr} > 1)                     = 1;                                               %

			unphased_plot{ chr}(chr_bin)                                 = length(chr_SNPdata{chr,2}{chr_bin});             % unphased data.
			unphased_plot2{chr}(chr_bin)                                 = unphased_plot{chr}(chr_bin)/full_data_threshold; %
			unphased_plot2{chr}(unphased_plot2{chr} > 1)                 = 1;                                               %
		end;
	end;
end;
fprintf('\n');

% load size definitions
[linear_fig_height,linear_fig_width,Linear_left_start,Linear_chr_gap,Linear_Chr_max_width,Linear_height...
    ,Linear_base,rotate,linear_chr_font_size,linear_axis_font_size,linear_gca_font_size,stacked_fig_height,...
    stacked_fig_width,stacked_chr_font_size,stacked_title_size,stacked_axis_font_size,...
    gca_stacked_font_size,stacked_copy_font_size,max_chrom_label_size] = Load_size_info(chr_in_use,num_chrs,chr_label,chr_size);

%%================================================================================================
% Setup for main-view figure generation.
%-------------------------------------------------------------------------------------------------
fprintf('\t|\tIntitial setup of main figure.\n');
fig = figure(1);
%set(gcf, 'Position', [0 70 1024 600]);
% basic plot parameters not defined per genome.
TickSize         = -0.005;  %negative for outside, percentage of longest chr figure.
maxY             = ploidyBase*2;
cen_tel_Xindent  = 5;
cen_tel_Yindent  = maxY/5;


%%================================================================================================
% Setup for linear-view figure generation.
%-------------------------------------------------------------------------------------------------
if (Linear_display == true)
	fprintf('\t|\tInitial setup of linear figure.\n');
	Linear_fig = figure(2);
	Linear_genome_size   = sum(chr_size);
	Linear_TickSize      = -0.01;  %negative for outside, percentage of longest chr figure.
	maxY                 = ploidyBase*2;
	Linear_left          = Linear_left_start;
	axisLabelPosition_horiz = 0.01125;
end;
axisLabelPosition_vert = 0.01125;


%%================================================================================================
% Make figures
%-------------------------------------------------------------------------------------------------
fprintf('\t|\tDraw figures.\n');
first_chr = true;
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		figure(fig);
		% make standard chr cartoons.
		left          = chr_posX(chr);
		bottom        = chr_posY(chr);
		width         = chr_width(chr);
		height        = chr_height(chr);
		subPlotHandle = subplot('Position',[left bottom width height]);
		fprintf(['\tfigposition = [' num2str(left) ' | ' num2str(bottom) ' | ' num2str(width) ' | ' num2str(height) ']\n']);
		hold on;

		c_prev = colorInit;
		c_post = colorInit;
		c_     = c_prev;
		infill = zeros(1,length(phased_plot2{chr}));
		colors = [];


		%% standard : determine color of each bin.
		for chr_bin = 1:ceil(chr_size(chr)/bases_per_bin)
			c_tot_post = SNPs_to_fullData_ratio{chr}(chr_bin)+SNPs_to_fullData_ratio{chr}(chr_bin);
			if (c_tot_post == 0)
				c_post = colorNoData;
				fprintf('.');
				if (mod(chr_bin,100) == 0);   fprintf('\n');   end;
			else
				% Average of SNP position colors defined earlier.
				colorMix = [chr_SNPdata_colorsC{chr,1}(chr_bin) chr_SNPdata_colorsC{chr,2}(chr_bin) chr_SNPdata_colorsC{chr,3}(chr_bin)];

				% Determine color to draw bin, accounting for limited data and data saturation.
				c_post =   colorMix   *   min(1,SNPs_to_fullData_ratio{chr}(chr_bin)) + ...
				           colorNoData*(1-min(1,SNPs_to_fullData_ratio{chr}(chr_bin)));
			end;
			colors(chr_bin,1) = c_post(1);
			colors(chr_bin,2) = c_post(2);
			colors(chr_bin,3) = c_post(3);
		end;
		% standard : end determine color of each bin.


		%% standard : draw colorbars.
		for chr_bin = 1:ceil(chr_size(chr)/bases_per_bin)
			x_ = [chr_bin chr_bin chr_bin-1 chr_bin-1];
			y_ = [0 maxY maxY 0];
			c_post(1) = colors(chr_bin,1);
			c_post(2) = colors(chr_bin,2);
			c_post(3) = colors(chr_bin,3);
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

		% standard : axes labels etc.
		hold off;
		xlim([0,chr_size(chr)/bases_per_bin]);

		%% standard : modify y axis limits to show annotation locations if any are provided.
		if (length(annotations) > 0)
			ylim([-maxY/10*1.5,maxY]);
		else
			ylim([0,maxY]);
		end;
		set(gca,'TickLength',[(TickSize*chr_size(largestChr)/chr_size(chr)) 0]); %ensures same tick size on all subfigs.
		set(gca,'YTick',[]);
		set(gca,'YTickLabel',[]);
		set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
		set(gca,'XTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8','2.0','2.2','2.4','2.6','2.8','3.0','3.2'});
		text(-50000/5000/2*3, maxY/2,     chr_label{chr}, 'Rotation',90, 'HorizontalAlignment','center', 'VerticalAlign','bottom', 'Fontsize',stacked_chr_font_size);

		set(gca,'FontSize',gca_stacked_font_size);
		if (chr == find(chr_posY == max(chr_posY)))
			title([ project ' vs. (hapmap)' hapmap ' SNP/LOH map'],'Interpreter','none','FontSize',stacked_title_size);
		end;
		hold on;
		% standard : end axes labels etc.
    
		% standard : show segmental anueploidy breakpoints.
                if (displayBREAKS == true) && (show_annotations == true)
                        chr_length = ceil(chr_size(chr)/bases_per_bin);
                        for segment = 2:length(chr_breaks{chr})-1
                                bP = chr_breaks{chr}(segment)*chr_length;
                                plot([bP bP], [(-maxY/10*2.5) 0],  'Color',[1 0 0],'LineWidth',2);
                        end;
                end;

		% show centromere outlines and horizontal marks.
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
			for annoteID = 1:length(annotation_location)
				if (annotation_chr(annoteID) == chr)
					annotationloc   = annotation_location(annoteID)/bases_per_bin-0.5*(5000/bases_per_bin);
					annotationStart = annotation_start(   annoteID)/bases_per_bin-0.5*(5000/bases_per_bin);
					annotationEnd   = annotation_end(     annoteID)/bases_per_bin-0.5*(5000/bases_per_bin);
					if (strcmp(annotation_type{annoteID},'dot') == 1)
						plot(annotationloc,-maxY/10*1.5,'k:o','MarkerEdgeColor',annotation_edgecolor{annoteID}, ...
						     'MarkerFaceColor',annotation_fillcolor{annoteID}, ...
						     'MarkerSize',     annotation_size(annoteID));
					elseif (strcmp(annotation_type{annoteID},'block') == 1)
						fill([annotationStart annotationStart annotationEnd annotationEnd], ...
						     [-maxY/10*(1.5+0.75) -maxY/10*(1.5-0.75) -maxY/10*(1.5-0.75) -maxY/10*(1.5+0.75)], ...
						     annotation_fillcolor{annoteID},'EdgeColor',annotation_edgecolor{annoteID});
					end;
				end;
			end;
			hold off;
		end;
		%end show annotation locations.


		%% =========================================================================================
		% Draw angleplots to left of main chromosome cartoons.
		%-------------------------------------------------------------------------------------------
		apply_phasing = true;
		angle_plot_subfigures;


%%%%%%%%%%%%%%%% Linear figure draw section


		%% Linear figure draw section
		if (Linear_display == true)
			figure(Linear_fig);
			Linear_width = Linear_Chr_max_width*chr_size(chr)/Linear_genome_size;
			subplot('Position',[Linear_left Linear_base Linear_width Linear_height]);
			Linear_left = Linear_left + Linear_width + Linear_chr_gap;
			hold on;

			%% linear : draw colorbars.
			for chr_bin = 1:ceil(chr_size(chr)/bases_per_bin)
				x_ = [chr_bin chr_bin chr_bin-1 chr_bin-1];
				y_ = [0 maxY maxY 0];
				c_post(1) = colors(chr_bin,1);
				c_post(2) = colors(chr_bin,2);
				c_post(3) = colors(chr_bin,3);
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
			% linear : end draw colorbars.

			%% linear : show segmental anueploidy breakpoints.
			if (Linear_displayBREAKS == true) && (show_annotations == true)
				chr_length = ceil(chr_size(chr)/bases_per_bin);
				for segment = 2:length(chr_breaks{chr})-1
					bP = chr_breaks{chr}(segment)*chr_length;
					plot([bP bP], [(-maxY/10*2.5) 0],  'Color',[1 0 0],'LineWidth',2);
				end;
			end;

			%% linear : show centromere.
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
			% linear : end show centromere.

			%% linear : show annotation locations
			if (show_annotations) && (length(annotations) > 0)
				plot([leftEnd rightEnd], [-maxY/10*1.5 -maxY/10*1.5],'color',[0 0 0]);
				hold on;
				annotation_location = (annotation_start+annotation_end)./2;
				for annoteID = 1:length(annotation_location)
					if (annotation_chr(annoteID) == chr)
						annotationloc   = annotation_location(annoteID)/bases_per_bin-0.5*(5000/bases_per_bin);
						annotationStart = annotation_start(   annoteID)/bases_per_bin-0.5*(5000/bases_per_bin);
						annotationEnd   = annotation_end(     annoteID)/bases_per_bin-0.5*(5000/bases_per_bin);
						if (strcmp(annotation_type{annoteID},'dot') == 1)
							plot(annotationloc,-maxY/10*1.5,'k:o','MarkerEdgeColor',annotation_edgecolor{annoteID}, ...
							     'MarkerFaceColor',annotation_fillcolor{annoteID}, ...
							     'MarkerSize',     annotation_size(annoteID));
						elseif (strcmp(annotation_type{annoteID},'block') == 1)
							fill([annotationStart annotationStart annotationEnd annotationEnd], ...
							     [-maxY/10*(1.5+0.75) -maxY/10*(1.5-0.75) -maxY/10*(1.5-0.75) -maxY/10*(1.5+0.75)], ...
							     annotation_fillcolor{annoteID},'EdgeColor',annotation_edgecolor{annoteID});
						end;
					end;
				end;
				hold off;
			end;
			% linear : end show annotation locations.

			%% linear : Final formatting stuff.
			xlim([0,chr_size(chr)/bases_per_bin]);
			% modify y axis limits to show annotation locations if any are provided.
			if (length(annotations) > 0)
				ylim([-maxY/10*1.5,maxY]);
			else
				ylim([0,maxY]);
			end;
			set(gca,'TickLength',[(Linear_TickSize*chr_size(largestChr)/chr_size(chr)) 0]); %ensures same tick size on all subfigs.
			set(gca,'YTick',[]);
        		set(gca,'YTickLabel',[]);
			set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
			set(gca,'XTickLabel',[]);
			set(gca,'FontSize',linear_gca_font_size);
			% end final reformatting.
			% adding title in the middle of the cartoon
			% note: adding title is done in the end since if placed upper
			% in the code somehow the plot function changes the title position			
			if (rotate == 0 && chr_size(chr) ~= 0 )
				title(chr_label{chr},'Interpreter','none','FontSize',linear_chr_font_size,'Rotation',rotate);
			else
				text((chr_size(chr)/bases_per_bin)/2,maxY+0.5,chr_label{chr},'Interpreter','none','FontSize',linear_chr_font_size,'Rotation',rotate);
			end;
	        
			% shift back to main figure generation.
			figure(fig);
			hold on;
			first_chr = false;
		end;
	end;
end;

%% Save figures.
set(fig,'PaperPosition',[0 0 stacked_fig_width stacked_fig_height]);
saveas(fig,        [projectDir 'fig.SNP-map.1.eps'], 'epsc');
saveas(fig,        [projectDir 'fig.SNP-map.1.png'], 'png');
delete(fig);

set(Linear_fig,'PaperPosition',[0 0 linear_fig_width linear_fig_height]);
saveas(Linear_fig, [projectDir 'fig.SNP-map.2.eps'], 'epsc');
saveas(Linear_fig, [projectDir 'fig.SNP-map.2.png'], 'png');
delete(Linear_fig);

%% ========================================================================
% end stuff
%==========================================================================
end
