function [] = allelic_ratios_ddRADseq_B(main_dir,user,genomeUser,project,parent,hapmap,genome,ploidyEstimateString,ploidyBaseString,SNP_verString,LOH_verString,CNV_verString,displayBREAKS);
addpath('../');
workingDir = [main_dir 'users/' user '/projects/' project '/'];

%% ================================================================================================
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
Linear_displayBREAKS        = false;

AnglePlot                   = true;   % Show histogram of alleleic fraction at the left end of standard figure chromosomes.
FillColors                  = true;   %     Fill histogram using colors.
show_uncalibrated           = false;  %     Fill with single color instead of ratio call colors.


fprintf('\n');
fprintf('#################################\n');
fprintf('## allelic_ratios_ddRADseq_B.m ##\n');
fprintf('#################################\n');


%%================================================================================================
% Load FASTA file name from 'reference.txt' file for project.
%-------------------------------------------------------------------------------------------------
userReference    = [main_dir 'users/' user '/genomes/' genome '/reference.txt'];
defaultReference = [main_dir 'users/default/genomes/' genome '/reference.txt'];
if (exist(userReference,'file') == 0)   
	FASTA_string = strtrim(fileread(defaultReference));
else                    
	FASTA_string = strtrim(fileread(userReference));
end;
[FastaPath,FastaName,FastaExt] = fileparts(FASTA_string);


%%================================================================================================
% Control variables.
%-------------------------------------------------------------------------------------------------
projectDir = [main_dir 'users/' user '/projects/' project '/'];
genomeDir  = [main_dir 'users/' genomeUser '/genomes/' genome '/'];
if (strcmp(hapmap,'') == 1)
	useHapmap = false;
else
	useHapmap = true;
	if (exist([main_dir 'users/default/hapmaps/' hapmap '/'], 'dir') == 7)
		hapmapDir  = [main_dir 'users/default/hapmaps/' hapmap '/'];   % system hapmap.
		hapmapUser = 'default';
	else
		hapmapDir  = [main_dir 'users/' user '/hapmaps/' hapmap '/'];  % user hapmap.
		hapmapUser = user;
	end;
end;
if (strcmp(project,parent) == 1)
	useParent  = false;
	parentDir  = projectDir;
	parentUSer = user;
else
	useParent = true;
	if (exist([main_dir 'users/default/projects/' parent '/'], 'dir') == 7)
		parentDir  = [main_dir 'users/default/projects/' parent '/'];   % system parent.
		parentUser = 'default';
	else
		parentDir  = [main_dir 'users/' user '/projects/' parent '/'];  % user parent.
		parentUser = user;
	end;
end;


genomeDir  = [main_dir 'users/' genomeUser '/genomes/' genome '/'];
[centromeres, chr_sizes, figure_details, annotations, ploidy_default] = Load_genome_information(genomeDir);
[Aneuploidy]                                                          = Load_dataset_information(projectDir);
%% check if dataset uses original chromosomes names
originalNamePath = [projectDir 'original.txt'];
if (exist(originalNamePath,'file'))
	useOriginal = true;
else 
	useOriginal = false;
end;
num_chrs = length(chr_sizes);

for chrID = 1:length(chr_sizes)
	chr_size(chrID)  = 0;
	cen_start(chrID) = 0;
	cen_end(chrID)   = 0;
end;
for chrID = 1:length(chr_sizes)
	chr_size(chr_sizes(chrID).chr)    = chr_sizes(chrID).size;
	cen_start(centromeres(chrID).chr) = centromeres(chrID).start;
	cen_end(centromeres(chrID).chr)   = centromeres(chrID).end;
end;
if (length(annotations) > 0)
	fprintf(['\nAnnotations for ' genome '.\n']);
	for annoteID = 1:length(annotations)
		annotation_chr(annoteID)       = annotations(annoteID).chr;
		annotation_type{annoteID}      = annotations(annoteID).type;
		annotation_start(annoteID)     = annotations(annoteID).start;
		annotation_end(annoteID)       = annotations(annoteID).end;
		annotation_fillcolor{annoteID} = annotations(annoteID).fillcolor;
		annotation_edgecolor{annoteID} = annotations(annoteID).edgecolor;
		annotation_size(annoteID)      = annotations(annoteID).size;
		fprintf(['\t[' num2str(annotations(annoteID).chr) ':' annotations(annoteID).type ':' num2str(annotations(annoteID).start) ':' num2str(annotations(annoteID).end) ':' annotations(annoteID).fillcolor ':' annotations(annoteID).edgecolor ':' num2str(annotations(annoteID).size) ']\n']);
	end;
end;
for detailID = 1:length(figure_details)
	if (figure_details(detailID).chr == 0)
		if (strcmp(figure_details(detailID).label,'Key') == 1)
			key_posX   = figure_details(detailID).posX;
			key_posY   = figure_details(detailID).posY;
			key_width  = figure_details(detailID).width;
			key_height = figure_details(detailID).height;
		end;
	else
		chr_id    (figure_details(detailID).chr) = figure_details(detailID).chr;
		if (useOriginal && length(figure_details(detailID).name) < 10)
		    chr_label {figure_details(detailID).chr} = figure_details(detailID).name;
		else
		    chr_label {figure_details(detailID).chr} = figure_details(detailID).label;
		end;
		chr_name  {figure_details(detailID).chr} = figure_details(detailID).name;
		chr_posX  (figure_details(detailID).chr) = figure_details(detailID).posX;
		chr_posY  (figure_details(detailID).chr) = figure_details(detailID).posY;
		chr_width (figure_details(detailID).chr) = figure_details(detailID).width;
		chr_height(figure_details(detailID).chr) = figure_details(detailID).height;
		chr_in_use(figure_details(detailID).chr) = str2num(figure_details(detailID).useChr);
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
			for chrBreak = 1:length(Aneuploidy)
				if (Aneuploidy(chrBreak).chr == usedChr)
					break_count = break_count+1;
					chr_broken  = true;
					chr_breaks{usedChr}(break_count) = Aneuploidy(chrBreak).break;
				end;
			end;
		end;
		chr_breaks{usedChr}(length(chr_breaks{usedChr})+1) = 1;
	end;
end;


%% =========================================================================================
%% =========================================================================================
%% =========================================================================================
%% = No further control variables below. ===================================================
%% =========================================================================================
%% =========================================================================================
%% =========================================================================================


% Process input ploidy.
ploidy = str2num(ploidyEstimateString);

% Sanitize user input of euploid state.
ploidyBase = round(str2num(ploidyBaseString));
if (ploidyBase > 4);   ploidyBase = 4;   end;
if (ploidyBase < 1);   ploidyBase = 1;   end;
fprintf(['\nEuploid base = "' num2str(ploidyBase) '"\n']);

% basic plot parameters not defined per genome.
TickSize         = -0.005;  %negative for outside, percentage of longest chr figure.
bases_per_bin    = max(chr_size)/700;
maxY             = 1; % ploidyBase*2;
cen_tel_Xindent  = 5;
cen_tel_Yindent  = maxY/5;


%% =========================================================================================
% Define colors for figure generation.
%-------------------------------------------------------------------------------------------
fprintf('\t|\tDefine colors used in figure generation.\n');
phased_and_unphased_color_definitions;


%%================================================================================================
% Load 'Common_CNV.mat' file containing CNV estimates per standard genome bin.
%-------------------------------------------------------------------------------------------------
fprintf('\nLoading "Common_CNV" data file for ddRADseq project.');
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
% Setup for SNP/LOH data calculations.
%-------------------------------------------------------------------------------------------------
fprintf(['\nGenerating LOH-map figure from ''' project ''' vs. (hapmap)''' hapmap ''' data.\n']);
% Initializes vectors used to hold allelic ratios for each chromosome segment.
if (useHapmap)
	new_bases_per_bin = bases_per_bin;
else
	new_bases_per_bin = bases_per_bin/2;
end;
for chr = 1:length(chr_sizes)
	chr_length = ceil(chr_size(chr)/new_bases_per_bin);

	% chr_SNPdata{chr,i} :: i = [1..4]
	%       1 : phased SNP ratio data.
	%       2 : unphased SNP ratio data.
	%       3 : phased SNP position data.
	%       4 : unphased SNP position data.
	%       5 : phased SNP flipper value.
	%       6 : unphased SNP flipper value.
	for j = 1:6
		for i = 1:chr_length
			chr_SNPdata{chr,j}{i} = [];
		end;
	end;

	% Colors used to illustrate SNP/LOH data.
	%       chr_SNPdata_colorsC           : colors scheme defined by hapmap or red for unspecified LOH.
	%       chr_SNPdata_colorsC{chr,i} :: i = [1..3] = [R, G, B]
	%       chr_SNPdata_colorsP{chr,i} :: i = [1..3] = [R, G, B]
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
% Load SNP/LOH data.
%-------------------------------------------------------------------------------------------------
if (useHapmap)
	if (exist([projectDir 'SNP_' SNP_verString '.all3.mat'],'file') == 0)
		fprintf('\nAllelic fraction MAT file 3 not found, generating.\n');
		process_2dataset_hapmap_allelicRatios(projectDir, projectDir, hapmapDir, chr_size, chr_name, chr_in_use, SNP_verString);
	else
		fprintf('\nAllelic fraction MAT file 3 found, loading.\n');
	end;
	load([projectDir 'SNP_' SNP_verString '.all3.mat']);
	% child data:  'C_chr_SNP_data_positions','C_chr_SNP_data_ratios','C_chr_count','C_chr_baseCall','C_chr_SNP_homologA','C_chr_SNP_homologB','C_chr_SNP_flipHomologs'
	%
	% C_chr_SNP_data_positions = coordinate of SNP.
	% C_chr_SNP_data_ratios    = allelic ratio of SNP.
	% C_chr_count              = number of reads at SNP coordinate.
	% C_chr_baseCall           = majority basecall of SNP.
	% C_chr_SNP_homologA       = hapmap homolog a basecall.
	% C_chr_SNP_homologB       = hapmap homolog b basecall.
	% C_chr_SNP_flipHomologs   = does hapmap entry need flipped?   0 = 'correct phase, no', 1 = 'incorrect phase, yes', 10 = 'no phasing info'.
else
	if (exist([projectDir 'SNP_' SNP_verString '.all1.mat'],'file') == 0)
		fprintf('\nAllelic fraction MAT file 1 not found, generating.\n');
		process_2dataset_allelicRatios(projectDir, parentDir, chr_size, chr_name, chr_in_use, SNP_verString);
	else
		fprintf('\nAllelic fraction MAT file 1 found, loading.\n');
	end;
	load([projectDir 'SNP_' SNP_verString '.all1.mat']);
	% child data:  'C_chr_SNP_data_positions','C_chr_SNP_data_ratios','C_chr_count'
	% parent data: 'P_chr_SNP_data_positions','P_chr_SNP_data_ratios','P_chr_count'
	%   
	% C_chr_SNP_data_positions = coordinate of SNP.
	% C_chr_SNP_data_ratios    = allelic ratio of SNP.
	% C_chr_count              = number of reads at SNP coordinate.
end;


%%================================================================================================
% Save workspace variables for use in "allelic_ratios_ddRADseq_D.m"
%-------------------------------------------------------------------------------------------------
save([projectDir 'allelic_ratios_ddRADseq_B.workspace_variables.mat']);


%% =================================================================================================
% Calculate allelic fraction cutoffs for each segment.
%...................................................................................................
% Populate data structure containing SNP phasing information.
% if (useHapmap)
%	chr_SNPdata{chr,1}{chr_bin} = phased SNP ratio data.
%	chr_SNPdata{chr,2}{chr_bin} = unphased SNP ratio data.
%	chr_SNPdata{chr,3}{chr_bin} = phased SNP position data.
%	chr_SNPdata{chr,4}{chr_bin} = unphased SNP position data.
%	chr_SNPdata{chr,5}{chr_bin} = flipper value for phased SNP.
%	chr_SNPdata{chr,6}{chr_bin} = flipper value for unphased SNP.
% elseif (useParent)
%       chr_SNPdata{chr,1}{chr_bin} = child SNP ratio data.
%       chr_SNPdata{chr,2}{chr_bin} = parent SNP ratio data.
%       chr_SNPdata{chr,3}{chr_bin} = child SNP position data.
%       chr_SNPdata{chr,4}{chr_bin} = parent SNP position data.
% else
%       chr_SNPdata{chr,2}{chr_bin} = child SNP ratio data.
%       chr_SNPdata{chr,4}{chr_bin} = child SNP position data.
% end;
%---------------------------------------------------------------------------------------------------
calculate_allelic_ratio_cutoffs;


%%================================================================================================
% Process SNP/LOH data.
%-------------------------------------------------------------------------------------------------
if (useHapmap)
	%% =========================================================================================
	% Define new colors for SNPs, using Gaussian fitting crossover points as ratio cutoffs.
	%-------------------------------------------------------------------------------------------
	for chr = 1:num_chrs
		if (chr_in_use(chr) == 1)
			if (length(C_chr_count{chr}) > 1)
				%
				% Determining colors for each SNP coordinate from calculated cutoffs.
				%
				for SNP = 1:length(C_chr_SNP_data_positions{chr})  % 'length(C_chr_SNP_data_positions)' is the number of SNPs per chromosome.
					coordinate                      = C_chr_SNP_data_positions{chr}(SNP);
					chr_bin                         = ceil(coordinate/new_bases_per_bin);
					localCopyEstimate               = round(CNVplot2{chr}(chr_bin)*ploidy*ploidyAdjust);
					baseCall                        = C_chr_baseCall{        chr}{SNP};
					homologA                        = C_chr_SNP_homologA{    chr}{SNP};
					homologB                        = C_chr_SNP_homologB{    chr}{SNP};
					flipper                         = C_chr_SNP_flipHomologs{chr}(SNP);
					allelic_ratio                   = C_chr_SNP_data_ratios{ chr}(SNP);
					if (flipper == 1)
						temp                    = homologA;
						homologA                = homologB;
						homologB                = temp;
						if (baseCall == homologA)
							allelic_ratio   = 1-allelic_ratio;
						end;
					elseif (flipper == 0)
						if (baseCall == homologA)
							allelic_ratio   = 1-allelic_ratio;
						end;
					else
						% Variable 'flipper' value of '10' indicates no phasing information is available in the hapmap.
						% Variable 'baseCall' value of 'Z' will prevent either hapmap allele from matching and so unphased ratio colors will be used in the following sections.
						baseCall                = 'Z';
					end;


					% identify the segment containing the SNP.
					segmentID                       = 0;
					for segment = 1:(length(chrCopyNum{chr}))
						segment_start           = chr_breaks{chr}(segment  )*chr_size(chr);
						segment_end             = chr_breaks{chr}(segment+1)*chr_size(chr);
						if (coordinate > segment_start) && (coordinate <= segment_end)
							segmentID       = segment;
						end;
					end;


					% Load cutoffs between Gaussian fits performed earlier.
					segment_copyNum                 = round(chrCopyNum{              chr}(segmentID));
					actual_cutoffs                  = chrSegment_actual_cutoffs{     chr}{segmentID};
					mostLikelyGaussians             = chrSegment_mostLikelyGaussians{chr}{segmentID};


					% Calculate allelic ratio on range of [1..200].
					SNPratio_int                    = (allelic_ratio)*199+1;


					% Identify the allelic ratio region containing the SNP.
					cutoffs                         = [1 actual_cutoffs 200];
					ratioRegionID                   = 0;
					for GaussianRegionID = 1:length(mostLikelyGaussians)
						cutoff_start            = cutoffs(GaussianRegionID  );
						cutoff_end              = cutoffs(GaussianRegionID+1);
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
							if (useHapmap || useParent)         colorList = unphased_color_1of1;
							else                                colorList = noparent_color_1of1;
							end;
						end;
					elseif (segment_copyNum == 2)
						%   allelic fraction cutoffs: [0.25000 0.75000] => [AA AB BB]
						%   cutoffs_1 from Gaussian fittings = [0.18342 0.19347 0.18844 0.16332 0.16332 0.15327 0.12814] => 0.1676;
						%   cutoffs_2 from Gaussian fittings = [0.82663 0.82161 0.80151 0.84171 0.83668 0.84673 0.84673] => 0.8317;
						if ((baseCall == homologA) || (baseCall == homologB))
							if (ratioRegionID == 3);            colorList = colorBB;
							elseif (ratioRegionID == 2);        colorList = colorAB;
							else                                colorList = colorAA;
							end;
						else
							if (useHapmap || useParent)
								if (ratioRegionID == 3);        colorList = unphased_color_2of2;
								elseif (ratioRegionID == 2);    colorList = unphased_color_1of2;
								else                            colorList = unphased_color_2of2;
								end;
							else
								if (ratioRegionID == 3);        colorList = noparent_color_2of2;
								elseif (ratioRegionID == 2);    colorList = noparent_color_1of2;
								else                            colorList = noparent_color_2of2;
								end;
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
							if (useHapmap || useParent)
								if (ratioRegionID == 4);        colorList = unphased_color_3of3;
								elseif (ratioRegionID == 3);    colorList = unphased_color_2of3;
								elseif (ratioRegionID == 2);    colorList = unphased_color_2of3;
								else                            colorList = unphased_color_3of3;
								end;
							else
								if (ratioRegionID == 4);        colorList = noparent_color_3of3;
								elseif (ratioRegionID == 3);    colorList = noparent_color_2of3;
								elseif (ratioRegionID == 2);    colorList = noparent_color_2of3;
								else                            colorList = noparent_color_3of3;
								end;
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
							if (useHapmap || useParent)
								if (ratioRegionID == 5);        colorList = unphased_color_4of4;
								elseif (ratioRegionID == 4);    colorList = unphased_color_3of4;
								elseif (ratioRegionID == 3);    colorList = unphased_color_2of4;
								elseif (ratioRegionID == 2);    colorList = unphased_color_3of4;
								else                            colorList = unphased_color_4of4;
								end;
							else
								if (ratioRegionID == 5);        colorList = noparent_color_4of4;
								elseif (ratioRegionID == 4);    colorList = noparent_color_3of4;
								elseif (ratioRegionID == 3);    colorList = noparent_color_2of4;
								elseif (ratioRegionID == 2);    colorList = noparent_color_3of4;
								else                            colorList = noparent_color_4of4;
								end;
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
							if (useHapmap || useParent)
								if (ratioRegionID == 6);        colorList = unphased_color_5of5;
								elseif (ratioRegionID == 5);    colorList = unphased_color_4of5;
								elseif (ratioRegionID == 4);    colorList = unphased_color_3of5;
								elseif (ratioRegionID == 3);    colorList = unphased_color_3of5;
								elseif (ratioRegionID == 2);    colorList = unphased_color_4of5;
								else                            colorList = unphased_color_5of5;
								end;
							else
								if (ratioRegionID == 6);        colorList = noparent_color_5of5;
								elseif (ratioRegionID == 5);    colorList = noparent_color_4of5;
								elseif (ratioRegionID == 4);    colorList = noparent_color_3of5;
								elseif (ratioRegionID == 3);    colorList = noparent_color_3of5;
								elseif (ratioRegionID == 2);    colorList = noparent_color_4of5;
								else                            colorList = noparent_color_5of5;
								end;
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
							if (useHapmap || useParent)
								if (ratioRegionID == 7);        colorList = unphased_color_6of6;
								elseif (ratioRegionID == 6);    colorList = unphased_color_5of6;
								elseif (ratioRegionID == 5);    colorList = unphased_color_4of6;
								elseif (ratioRegionID == 4);    colorList = unphased_color_3of6;
								elseif (ratioRegionID == 3);    colorList = unphased_color_4of6;
								elseif (ratioRegionID == 2);    colorList = unphased_color_5of6;
								else                            colorList = unphased_color_6of6;
								end;
							else
								if (ratioRegionID == 7);        colorList = noparent_color_6of6;
								elseif (ratioRegionID == 6);    colorList = noparent_color_5of6;
								elseif (ratioRegionID == 5);    colorList = nopanret_color_4of6;
								elseif (ratioRegionID == 4);    colorList = nopanret_color_3of6;
								elseif (ratioRegionID == 3);    colorList = noparent_color_4of6;
								elseif (ratioRegionID == 2);    colorList = noparent_color_5of6;
								else                            colorList = noparent_color_6of6;
								end;
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
							if (useHapmap || useParent)
								if (ratioRegionID == 8);        colorList = unphased_color_7of7;
								elseif (ratioRegionID == 7);    colorList = unphased_color_6of7;
								elseif (ratioRegionID == 6);    colorList = unphased_color_5of7;
								elseif (ratioRegionID == 5);    colorList = unphased_color_4of7;
								elseif (ratioRegionID == 3);    colorList = unphased_color_4of7;
								elseif (ratioRegionID == 3);    colorList = unphased_color_5of7;
								elseif (ratioRegionID == 2);    colorList = unphased_color_6of7;
								else                            colorList = unphased_color_7of7;
								end;
							else
								if (ratioRegionID == 8);        colorList = noparent_color_7of7;
								elseif (ratioRegionID == 7);    colorList = noparent_color_6of7;
								elseif (ratioRegionID == 6);    colorList = noparent_color_5of7;
								elseif (ratioRegionID == 5);    colorList = noparent_color_4of7;
								elseif (ratioRegionID == 3);    colorList = noparent_color_4of7;
								elseif (ratioRegionID == 3);    colorList = noparent_color_5of7;
								elseif (ratioRegionID == 2);    colorList = noparent_color_6of7;
								else                            colorList = noparent_color_7of7;
								end;
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
							if (useHapmap || useParent)
								if (ratioRegionID == 9);        colorList = unphased_color_8of8;
								elseif (ratioRegionID == 8);    colorList = unphased_color_7of8;
								elseif (ratioRegionID == 7);    colorList = unphased_color_6of8;
								elseif (ratioRegionID == 6);    colorList = unphased_color_5of8;
								elseif (ratioRegionID == 5);    colorList = unphased_color_4of8;
								elseif (ratioRegionID == 4);    colorList = unphased_color_5of8;
								elseif (ratioRegionID == 3);    colorList = unphased_color_6of8;
								elseif (ratioRegionID == 2);    colorList = unphased_color_7of8;
								else                            colorList = unphased_color_8of8;
								end;
							else
								if (ratioRegionID == 9);        colorList = noparent_color_8of8;
								elseif (ratioRegionID == 8);    colorList = noparent_color_7of8;
								elseif (ratioRegionID == 7);    colorList = noparent_color_6of8;
								elseif (ratioRegionID == 6);    colorList = noparent_color_5of8;
								elseif (ratioRegionID == 5);    colorList = noparent_color_4of8;
								elseif (ratioRegionID == 4);    colorList = noparent_color_5of8;
								elseif (ratioRegionID == 3);    colorList = noparent_color_6of8;
								elseif (ratioRegionID == 2);    colorList = noparent_color_7of8;
								else                            colorList = noparent_color_8of8;
								end;
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
							if (useHapmap || useParent)
								if (ratioRegionID == 10);       colorList = unphased_color_9of9;
								elseif (ratioRegionID == 9);    colorList = unphased_color_8of9;
								elseif (ratioRegionID == 8);    colorList = unphased_color_7of9;
								elseif (ratioRegionID == 7);    colorList = unphased_color_6of9;
								elseif (ratioRegionID == 6);    colorList = unphased_color_5of9;
								elseif (ratioRegionID == 5);    colorList = unphased_color_5of9;
								elseif (ratioRegionID == 4);    colorList = unphased_color_6of9;
								elseif (ratioRegionID == 3);    colorList = unphased_color_7of9;
								elseif (ratioRegionID == 2);    colorList = unphased_color_8of9;
								else                            colorList = unphased_color_9of9;
								end;
							else
								if (ratioRegionID == 10);       colorList = noparent_color_9of9;
								elseif (ratioRegionID == 9);    colorList = noparent_color_8of9;
								elseif (ratioRegionID == 8);    colorList = noparent_color_7of9;
								elseif (ratioRegionID == 7);    colorList = noparent_color_6of9;
								elseif (ratioRegionID == 6);    colorList = noparent_color_5of9;
								elseif (ratioRegionID == 5);    colorList = noparent_color_5of9;
								elseif (ratioRegionID == 4);    colorList = noparent_color_6of9;
								elseif (ratioRegionID == 3);    colorList = noparent_color_7of9;
								elseif (ratioRegionID == 2);    colorList = noparent_color_8of9;
								else                            colorList = noparent_color_9of9;
								end;
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
	end;

	%%
	%% Average colors per standard genome gin.
	%%
	for chr = 1:num_chrs
		if (chr_in_use(chr) == 1)
			if (length(C_chr_count{chr}) > 1)
				for chr_bin = 1:length(chr_SNPdata_countC{chr})
					if (chr_SNPdata_countC{chr}(chr_bin) > 0)
						chr_SNPdata_colorsC{chr,1}(chr_bin) = chr_SNPdata_colorsC{chr,1}(chr_bin)/chr_SNPdata_countC{chr}(chr_bin);
						chr_SNPdata_colorsC{chr,2}(chr_bin) = chr_SNPdata_colorsC{chr,2}(chr_bin)/chr_SNPdata_countC{chr}(chr_bin);
						chr_SNPdata_colorsC{chr,3}(chr_bin) = chr_SNPdata_colorsC{chr,3}(chr_bin)/chr_SNPdata_countC{chr}(chr_bin);
					else
						chr_SNPdata_colorsC{chr,1}(chr_bin) = 1.0;
						chr_SNPdata_colorsC{chr,2}(chr_bin) = 1.0;
						chr_SNPdata_colorsC{chr,3}(chr_bin) = 1.0;
					end;
				end;
			end;
		end;
	end;
elseif (useParent)
	% Initialize values for display.
	for chr = 1:num_chrs
		if (chr_in_use(chr) == 1)
			for chr_bin = 1:length(chr_SNPdata_countC{chr})
				colorList                       = [1.0 1.0 1.0];
				chr_SNPdata_colorsC{chr,1}(chr_bin) = colorList(1);
				chr_SNPdata_colorsC{chr,2}(chr_bin) = colorList(2);
				chr_SNPdata_colorsC{chr,3}(chr_bin) = colorList(3);

				colorList                       = [1.0 1.0 1.0];
				chr_SNPdata_colorsP{chr,1}(chr_bin) = colorList(1);
				chr_SNPdata_colorsP{chr,2}(chr_bin) = colorList(2);
				chr_SNPdata_colorsP{chr,3}(chr_bin) = colorList(3);
			end;
		end;
	end;
else
	% Initialize values for display.
	for chr = 1:num_chrs
		if (chr_in_use(chr) == 1)
			for chr_bin = 1:length(chr_SNPdata_countC{chr})
				colorList                       = [1.0 1.0 1.0];
				chr_SNPdata_colorsC{chr,1}(chr_bin) = colorList(1);
				chr_SNPdata_colorsC{chr,2}(chr_bin) = colorList(2);
				chr_SNPdata_colorsC{chr,3}(chr_bin) = colorList(3);
			end;
		end;
	end;
end;
save([projectDir 'SNP_' SNP_verString '.reduced.mat'],'chr_SNPdata','new_bases_per_bin','chr_SNPdata_colorsC','chr_SNPdata_colorsP');



%%================================================================================================
% Make histogram of data.
%-------------------------------------------------------------------------------------------------
histogram_fig = figure();
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		%   1 : phased SNP ratio data.
		%   2 : unphased SNP ratio data.
		%   3 : phased SNP position data.
		%   4 : unphased SNP position data.
		data_phased   = [];
		data_unphased = [];
		for i = 1:length(chr_SNPdata{chr,1})
			data_phased   = [data_phased   chr_SNPdata{chr,1}{i}];
		end;
		for i = 1:length(chr_SNPdata{chr,2})
			data_unphased = [data_unphased chr_SNPdata{chr,2}{i}];
		end;
		if (useHapmap)
		elseif (useParent)
			data_phased   = [data_phased   1-data_phased  ];
			data_unphased = [data_unphased 1-data_unphased];
		else
			data_phased   = [data_phased   1-data_phased  ];
			data_unphased = [data_unphased 1-data_unphased];
		end;
		histogram_phased         = hist([data_phased   0 1],200);
		histogram_unphased       = hist([data_unphased 0 1],200);

		subplot(2,num_chrs,chr);
		hold on;
		if (useHapmap)
			plot(1:200, log(histogram_unphased+1),  'Color',[1.0 0.0 0.0]);
			plot(1:200, log(histogram_phased+1),    'Color',[1/3 1/3 1/3]);
		else
			plot(1:200, log(histogram_unphased+1),  'Color',[1/3 1/3 1/3]);
			plot(1:200, log(histogram_phased+1),    'Color',[1.0 0.0 0.0]);
		end;
		ylim([0 6]);
		title([chr_label{chr}]);
        set(gca,'FontSize',10*(8/num_chrs));
		set(gca,'XTick',[0 50 100 150 200]);
		set(gca,'XTickLabel',{'0','1/4','1/2','3/4','1'});
		ylabel('log(data count)');
		xlabel('allelic ratio');
		hold off;

		subplot(2,num_chrs,chr+num_chrs);
		hold on;
		if (useHapmap)
			plot(1:200, histogram_unphased,  'Color',[1.0 0.0 0.0]);
			plot(1:200, histogram_phased,    'Color',[1/3 1/3 1/3]);
		else
			plot(1:200, histogram_unphased,  'Color',[1/3 1/3 1/3]);
			plot(1:200, histogram_phased,    'Color',[1.0 0.0 0.0]);
		end;
		ylim([0 200]);
		title([chr_label{chr}]);
        set(gca,'FontSize',10*(8/num_chrs));
		set(gca,'XTick',[0 50 100 150 200]);
		set(gca,'XTickLabel',{'0','1/4','1/2','3/4','1'});
		ylabel('data count');
		xlabel('allelic ratio');
		hold off;
	end;
end;
set(histogram_fig,'PaperPosition',[0 0 8*(num_chrs/8)^2 1.5*(num_chrs/8)]*4);
%saveas(histogram_fig, [projectDir 'fig.allelic_fraction_histogram.eps'], 'epsc');
saveas(histogram_fig, [projectDir 'fig.allelic_fraction_histogram.png'], 'png');
delete(histogram_fig);


%%================================================================================================
% Setup for main figure generation.
%-------------------------------------------------------------------------------------------------
% load size definitions
[linear_fig_height,linear_fig_width,Linear_left_start,Linear_chr_gap,Linear_Chr_max_width,Linear_height...
    ,Linear_base,rotate,linear_chr_font_size,linear_axis_font_size,linear_gca_font_size,stacked_fig_height,...
    stacked_fig_width,stacked_chr_font_size,stacked_title_size,stacked_axis_font_size,...
    gca_stacked_font_size,stacked_copy_font_size,max_chrom_label_size] = Load_size_info(chr_in_use,num_chrs,chr_label,chr_size);

fig = figure(1);
largestChr = find(chr_width == max(chr_width));
largestChr = largestChr(1);


%%================================================================================================
% Setup for linear-view figure generation.
%-------------------------------------------------------------------------------------------------
if (Linear_display == true)
	Linear_fig = figure(2);
	Linear_genome_size   = sum(chr_size);
	Linear_TickSize      = -0.01;  %negative for outside, percentage of longest chr figure.
	maxY                 = 1; % ploidyBase*2;
	Linear_left          = Linear_left_start;
	axisLabelPosition_horiz = 0.01125;
end;
axisLabelPosition_vert = 0.01125;


%%================================================================================================
% Make figures
%-------------------------------------------------------------------------------------------------
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
		hold on;
		fprintf(['\tfigposition = [' num2str(left) ' | ' num2str(bottom) ' | ' num2str(width) ' | ' num2str(height) ']\n']);

		% standard : axes labels etc.
		xlim([0,chr_size(chr)/bases_per_bin]);
    
		%% standard : modify y axis limits to show annotation locations if any are provided.
		if (length(annotations) > 0)
			ylim([-maxY/10*1.5,maxY]);
		else
			ylim([0,maxY]);
		end;
		set(gca,'YTick',[]);
		set(gca,'YTickLabel',[]);
		set(gca,'TickLength',[(TickSize*chr_size(largestChr)/chr_size(chr)) 0]); %ensures same tick size on all subfigs.
		text(-50000/5000/2*3, maxY/2,     chr_label{chr}, 'Rotation',90, 'HorizontalAlignment','center', 'VerticalAlign','bottom', 'Fontsize',stacked_chr_font_size);
		set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
		set(gca,'XTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8','2.0','2.2','2.4','2.6','2.8','3.0','3.2'});

		% standard : This section sets the Y-axis labelling.
		axisLabelPosition = -50000/bases_per_bin;
		set(gca,'FontSize',12);
		if (chr == find(chr_posY == max(chr_posY)))	Linear_base          = 0.1;

			title([ project ' allelic fraction map'],'Interpreter','none','FontSize',stacked_title_size);
		end;
		% standard : end axes labels etc.


		% standard : draw colorbars.
		if (useHapmap)   % an experimental dataset vs. a hapmap.
			for chr_bin = 1:ceil(chr_size(chr)/new_bases_per_bin)
				snpColorR = chr_SNPdata_colorsC{chr,1}(chr_bin);
				snpColorG = chr_SNPdata_colorsC{chr,2}(chr_bin);
				snpColorB = chr_SNPdata_colorsC{chr,3}(chr_bin);
				if (snpColorR <= 1) || (snpColorG <= 1) || (snpColorB <= 1)
					plot([chr_bin chr_bin], [0 maxY],'Color',[snpColorR snpColorG snpColorB]);
				end;
			end;
		% Following variation will draw experimetnal vs. reference dataset like the above vs. hapmap figure.
		%	elseif (useParent)   % an experimental dataset vs. a reference dataset.
		%		for chr_bin = 1:ceil(chr_size(chr)/new_bases_per_bin)
		%			snpColorR   = chr_SNPdata_colorsP{chr,1}(chr_bin);
		%			snpColorG   = chr_SNPdata_colorsP{chr,2}(chr_bin);
		%			snpColorB   = chr_SNPdata_colorsP{chr,3}(chr_bin);
		%			if (snpColorR <= 1) || (snpColorG <= 1) || (snpColorB <= 1)
		%				plot([chr_bin chr_bin], [0 maxY],'Color',[snpColorR snpColorG snpColorB]);
		%			end;
		%		end;
		elseif (useParent)   % an experimental dataset vs. a reference dataset.
			for chr_bin = 1:ceil(chr_size(chr)/new_bases_per_bin)
				%       chr_SNPdata{chr,1}{chr_bin} = child SNP ratio data.
				%       chr_SNPdata{chr,2}{chr_bin} = parent SNP ratio data.
				%       chr_SNPdata{chr,3}{chr_bin} = child SNP position data.
				%       chr_SNPdata{chr,4}{chr_bin} = parent SNP position data.

				bin_data_C              = chr_SNPdata{chr,1}{chr_bin};
				for SNP = 1:length(bin_data_C)
					bin_data_C(SNP) = bin_data_C(SNP);
				end;
				allelic_ratio_C         = min(bin_data_C);
				datumY_C                = allelic_ratio_C*maxY;
				plot([chr_bin/2 chr_bin/2], [maxY datumY_C     ],'Color',hom_color);

				bin_data_P              = chr_SNPdata{chr,2}{chr_bin};
				for SNP = 1:length(bin_data_P)
					bin_data_P(SNP) = bin_data_P(SNP);
				end;
				allelic_ratio_P         = min(bin_data_P);
				datumY_P                = allelic_ratio_P*maxY;
				plot([chr_bin/2 chr_bin/2], [0    maxY-datumY_P],'Color',het_color);
			end;
		else   % only an experimental dataset.
			for chr_bin = 1:ceil(chr_size(chr)/new_bases_per_bin)
				%       chr_SNPdata{chr,1}{chr_bin} = SNP ratio data.
				%       chr_SNPdata{chr,3}{chr_bin} = SNP position data.

				bin_data = chr_SNPdata{chr,1}{chr_bin};
				for SNP = 1:length(bin_data)
					bin_data(SNP) = max(bin_data(SNP), 1-bin_data(SNP));
				end;
				allelic_ratio = min(bin_data);
				datumY_C = allelic_ratio*maxY;
				plot([chr_bin/2 chr_bin/2], [maxY datumY_C     ],'Color',het_color);
				plot([chr_bin/2 chr_bin/2], [0    maxY-datumY_C],'Color',het_color);
			end;
		end;
		% standard : end draw colorbars.


		if (displayBREAKS == true) && (show_annotations == true)
			chr_length = ceil(chr_size(chr)/bases_per_bin);
			for segment = 2:length(chr_breaks{chr})-1
				bP = chr_breaks{chr}(segment)*chr_length;
				plot([bP bP], [(-maxY/10*2.5) 0],  'Color',[1 0 0],'LineWidth',2);
			end;
		end;

		% standard : show centromere outlines and horizontal marks.
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
		% standard : end show centromere.
    
		% standard : show annotation locations
		if (show_annotations) && (length(annotations) > 0)
			plot([leftEnd rightEnd], [-maxY/10*1.5 -maxY/10*1.5],'color',[0 0 0]);
			annotation_location = (annotation_start+annotation_end)./2;
			for annoteID = 1:length(annotation_location)
				if (annotation_chr(annoteID) == chr)
					annotationloc = annotation_location(annoteID)/bases_per_bin-0.5*(5000/bases_per_bin);
					annotationStart = annotation_start(annoteID)/bases_per_bin-0.5*(5000/bases_per_bin);
					annotationEnd   = annotation_end(annoteID)/bases_per_bin-0.5*(5000/bases_per_bin);
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
		end;
		% standard : end show annotation locations.


		%% =========================================================================================
		% Draw angleplots to left of main chromosome cartoons.
		%-------------------------------------------------------------------------------------------
		apply_phasing = true;
		angle_plot_subfigures;


%%%%%%%%%%%%%%%% Linear figure draw section.


		%% Linear figure draw section
		if (Linear_display == true)
			figure(Linear_fig);
			Linear_width = Linear_Chr_max_width*chr_size(chr)/Linear_genome_size;
			subplot('Position',[Linear_left Linear_base Linear_width Linear_height]);
			hold on;
			Linear_left = Linear_left + Linear_width + Linear_chr_gap;

			% linear : draw colorbars
			if (useHapmap)   % an experimental dataset vs. a hapmap.
				for chr_bin = 1:ceil(chr_size(chr)/new_bases_per_bin)
					snpColorR   = chr_SNPdata_colorsC{chr,1}(chr_bin);
					snpColorG   = chr_SNPdata_colorsC{chr,2}(chr_bin);
					snpColorB   = chr_SNPdata_colorsC{chr,3}(chr_bin);
					if (snpColorR < 1) || (snpColorG < 1) || (snpColorB < 1)
						plot([chr_bin chr_bin], [0 maxY],'Color',[snpColorR snpColorG snpColorB]);
					end;
				end;
			% Following variation will draw experimetnal vs. reference dataset like the above vs. hapmap figure.
			%	elseif (useParent)   % an experimental dataset vs. a reference dataset.
			%		for i = 1:ceil(chr_size(chr)/new_bases_per_bin)
			%			snpColorR   = chr_SNPdata_colorsC{chr,1}(chr_bin);
			%			snpColorG   = chr_SNPdata_colorsC{chr,2}(chr_bin);
			%			snpColorB   = chr_SNPdata_colorsC{chr,3}(chr_bin);
			%			if (snpColorR < 1) || (snpColorG < 1) || (snpColorB < 1)
			%				plot([i i], [0 maxY],'Color',[snpColorR snpColorG snpColorB]);
			%			end;
			%		end;
			elseif (useParent)   % an experimental dataset vs. a reference dataset.
				for chr_bin = 1:ceil(chr_size(chr)/new_bases_per_bin)
					%       chr_SNPdata{chr,1}{chr_bin} = child SNP ratio data.
					%       chr_SNPdata{chr,2}{chr_bin} = parent SNP ratio data.
					%       chr_SNPdata{chr,3}{chr_bin} = child SNP position data.
					%       chr_SNPdata{chr,4}{chr_bin} = parent SNP position data.

					bin_data_C              = chr_SNPdata{chr,1}{chr_bin};
					for SNP = 1:length(bin_data_C)
						bin_data_C(SNP) = bin_data_C(SNP);
					end;
					allelic_ratio_C         = min(bin_data_C);
					datumY_C                = allelic_ratio_C*maxY;
					plot([chr_bin/2 chr_bin/2], [maxY datumY_C     ],'Color',hom_color);

					bin_data_P              = chr_SNPdata{chr,2}{chr_bin};
					for SNP = 1:length(bin_data_P)
						bin_data_P(SNP) = bin_data_P(SNP);
					end;
					allelic_ratio_P         = min(bin_data_P);
					datumY_P                = allelic_ratio_P*maxY;
					plot([chr_bin/2 chr_bin/2], [0    maxY-datumY_P],'Color',het_color);
				end;
			else   % only an experimental dataset.
				for chr_bin = 1:ceil(chr_size(chr)/new_bases_per_bin)
					%       chr_SNPdata{chr,2}{chr_bin} = SNP ratio data.
					%       chr_SNPdata{chr,4}{chr_bin} = SNP position data.

					bin_data = chr_SNPdata{chr,2}{chr_bin};
					for SNP = 1:length(bin_data)
						bin_data(SNP) = max(bin_data(SNP), 1-bin_data(SNP));
					end;
					allelic_ratio = min(bin_data);
					datumY_C = allelic_ratio*maxY;
					plot([chr_bin/2 chr_bin/2], [maxY datumY_C     ],'Color',het_color);
					plot([chr_bin/2 chr_bin/2], [0    maxY-datumY_C],'Color',het_color);
				end;
			end;
			% linear : end draw colorbars.


			if (Linear_displayBREAKS == true) && (show_annotations == true)
				chr_length = ceil(chr_size(chr)/bases_per_bin);
                                for segment = 2:length(chr_breaks{chr})-1
                                        bP = chr_breaks{chr}(segment)*chr_length;
                                        plot([bP bP], [(-maxY/10*2.5) 0],  'Color',[1 0 0],'LineWidth',2);
                                end;
                        end;

			% linear : show centromere.
			x1 = cen_start(chr)/bases_per_bin;
			x2 = cen_end(chr)/bases_per_bin;
			leftEnd  = 0.5*5000/bases_per_bin;
			rightEnd = (chr_size(chr) - 0.5*5000)/bases_per_bin;

			if (Centromere_format == 0)
				% standard chromosome cartoons in a way which will not cause segfaults when running via commandline.
				dx = cen_tel_Xindent; %5*5000/bases_per_bin;
				dy = cen_tel_Yindent; %maxY/10;
				% draw white triangles at corners and centromere locations.
				c_ = [1.0 1.0 1.0];
				% top left corner.
				x_ = [leftEnd   leftEnd   leftEnd+dx];        y_ = [maxY-dy   maxY      maxY        ];    f = fill(x_,y_,c_);    set(f,'linestyle','none');
				% bottom left corner.     
				x_ = [leftEnd   leftEnd   leftEnd+dx];        y_ = [dy        0         0           ];    f = fill(x_,y_,c_);    set(f,'linestyle','none');
				% top right corner.
				x_ = [rightEnd   rightEnd   rightEnd-dx];     y_ = [maxY-dy    maxY       maxY      ];    f = fill(x_,y_,c_);    set(f,'linestyle','none');
				% bottom right corner.
				x_ = [rightEnd   rightEnd   rightEnd-dx];     y_ = [dy         0          0         ];    f = fill(x_,y_,c_);    set(f,'linestyle','none');
				% top centromere.
				x_ = [x1-dx   x1        x2        x2+dx];     y_ = [maxY    maxY-dy   maxY-dy   maxY];    f = fill(x_,y_,c_);    set(f,'linestyle','none');
				% bottom centromere.
				x_ = [x1-dx   x1   x2   x2+dx];               y_ = [0       dy   dy   0    ];             f = fill(x_,y_,c_);    set(f,'linestyle','none');
				% draw outlines of chromosome cartoon.   (drawn after horizontal lines to that cartoon edges are not interrupted by horiz lines.
				plot([leftEnd   leftEnd   leftEnd+dx   x1-dx   x1        x2        x2+dx   rightEnd-dx   rightEnd   rightEnd   rightEnd-dx   x2+dx   x2   x1   x1-dx   leftEnd+dx   leftEnd],...
				     [dy        maxY-dy   maxY         maxY    maxY-dy   maxY-dy   maxY    maxY          maxY-dy    dy         0             0       dy   dy   0       0            dy],...
				      'Color',[0 0 0]);
			end;
			% linear : end show centromere.

			% linear : show annotation locations
			if (show_annotations) && (length(annotations) > 0)
				plot([leftEnd rightEnd], [-maxY/10*1.5 -maxY/10*1.5],'color',[0 0 0]);
				annotation_location = (annotation_start+annotation_end)./2;
				for annoteID = 1:length(annotation_location)
					if (annotation_chr(annoteID) == chr)
						annotationloc = annotation_location(annoteID)/bases_per_bin-0.5*(5000/bases_per_bin);
						annotationStart = annotation_start(annoteID)/bases_per_bin-0.5*(5000/bases_per_bin);
						annotationEnd   = annotation_end(annoteID)/bases_per_bin-0.5*(5000/bases_per_bin);
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
			end;
			% linear : end show annotation locations.

			% linear : Final formatting stuff.
			xlim([0,chr_size(chr)/bases_per_bin]);
			% modify y axis limits to show annotation locations if any are provided.
			if (length(annotations) > 0)
				ylim([-maxY/10*1.5,maxY]);
			else
				ylim([0,maxY]);
			end;
			set(gca,'YTick',[]);
			set(gca,'YTickLabel',[]);
			set(gca,'TickLength',[(Linear_TickSize*chr_size(largestChr)/chr_size(chr)) 0]); %ensures same tick size on all subfigs.
			set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
			set(gca,'XTickLabel',[]);
			set(gca,'FontSize',linear_gca_font_size);
			% linear : end final reformatting.
			% note: adding title is done in the end since if placed upper
			% in the code somehow the plot function changes the title position	
			% adding title in the middle of the cartoon
			if (rotate == 0 && chr_size(chr) ~= 0 )
				title(chr_label{chr},'Interpreter','none','FontSize',linear_chr_font_size,'Rotation',rotate);
			else
				text((chr_size(chr)/bases_per_bin)/2,maxY+0.5,chr_label{chr},'Interpreter','none','FontSize',linear_chr_font_size,'Rotation',rotate);
			end;

			hold off;
	        
			% shift back to main figure generation.
			figure(fig);
			first_chr = false;
		end;
	end;
end;

%% Save figures.
set(fig,'PaperPosition',[0 0 stacked_fig_width stacked_fig_height]);
%saveas(fig,        [projectDir 'fig.allelic_ratio-map.c1.eps'], 'epsc');
saveas(fig,        [projectDir 'fig.allelic_ratio-map.c1.png'], 'png');
delete(fig);

set(Linear_fig,'PaperPosition',[0 0 linear_fig_width linear_fig_height]);
%saveas(Linear_fig, [projectDir 'fig.allelic_ratio-map.c2.eps'], 'epsc');
saveas(Linear_fig, [projectDir 'fig.allelic_ratio-map.c2.png'], 'png');
delete(Linear_fig);

%%================================================================================================
% end stuff
%=================================================================================================
end

