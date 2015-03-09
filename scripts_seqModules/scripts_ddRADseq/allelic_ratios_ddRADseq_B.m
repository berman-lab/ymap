function [] = allelic_ratios_ddRADseq_B(main_dir,user,genomeUser,project,parent,hapmap,genome,ploidyEstimateString,ploidyBaseString,SNP_verString,LOH_verString,CNV_verString,displayBREAKS);
addpath('../');


%%================================================================================================
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


%%================================================================================================
% Load the colors defined for the hapmap being used.
%-------------------------------------------------------------------------------------------------
if (useHapmap)
	%% Load color names defined for hapmap;
	colorsFile = [hapmapDir 'colors.txt'];
	if (exist(colorsFile,'file') == 2)
		colors_fid = fopen([main_dir 'users/' hapmapUser '/hapmaps/' hapmap '/colors.txt'], 'r');
		% The swapped colors are to correct for a polarity mistake in the python preprocessing steps.
		%    correcting the error there would require reprocessing all current datasets.
		colorA_string = fgetl(colors_fid);
		colorB_string = fgetl(colors_fid);
		fclose(colors_fid);
	else
		colorB_string = 'red';
		colorA_string = 'red';
	end;
	fprintf(['\nHapmap colors:\n\tcolorA = ' colorA_string '\n\tcolorB = ' colorB_string '\n\n']);

	switch colorA_string
		case 'deep pink'
			homolog_a_color = [1.0 0.0 0.5];
		case 'magenta'
			homolog_a_color = [1.0 0.0 1.0];
		case 'electric indigo'
			homolog_a_color = [0.5 0.0 1.0];
		case 'blue'
			homolog_a_color = [0.0 0.0 1.0];
		case 'dodger blue'
			homolog_a_color = [0.0 0.5 1.0];
		case 'cyan'
			homolog_a_color = [0.0 1.0 1.0];
		case 'spring green'
			homolog_a_color = [0.0 1.0 0.5];
		case 'green'
			homolog_a_color = [0.0 1.0 0.0];
		case 'chartreuse'
			homolog_a_color = [0.5 1.0 0.0];
		case 'yellow'
			homolog_a_color = [1.0 1.0 0.0];
		case 'dark orange'
			homolog_a_color = [1.0 0.5 0.0];
		case 'red'
			homolog_a_color = [1.0 0.0 0.0];
	end;
	switch colorB_string
		case 'deep pink'
			homolog_b_color = [1.0 0.0 0.5];
		case 'magenta'
			homolog_b_color = [1.0 0.0 1.0];
		case 'electric indigo'
			homolog_b_color = [0.5 0.0 1.0];
		case 'blue'
			homolog_b_color = [0.0 0.0 1.0];
		case 'dodger blue'
			homolog_b_color = [0.0 0.5 1.0];
		case 'cyan'
			homolog_b_color = [0.0 1.0 1.0];
		case 'spring green'
			homolog_b_color = [0.0 1.0 0.5];
		case 'green'
			homolog_b_color = [0.0 1.0 0.0];
		case 'chartreuse'
			homolog_b_color = [0.5 1.0 0.0];
		case 'yellow'
			homolog_b_color = [1.0 1.0 0.0];
		case 'dark orange'
			homolog_b_color = [1.0 0.5 0.0];
		case 'red'
			homolog_b_color = [1.0 0.0 0.0];
	end;

	het_color             = [0.66667 0.66667 0.66667]; % heterozygous.
	hom_unphased_color    = [1.0     0.0     0.0    ]; % homozygous, unphased.
	het_unphased_color    = [0.66667 0.66667 0.66667]; % heterozygous.
else
	% Haplotype map is not in use.
	if (strcmp(project,hapmap) == 1)
		% The 'project' is the same as the 'hapmap'/'parent'.
		homolog_a_color       = [0.66667 0.66667 0.66667];
		homolog_b_color       = [0.66667 0.66667 0.66667];
		het_color             = [0.66667 0.66667 0.66667]; % heterozygous.
		hom_unphased_color    = [0.66667 0.66667 0.66667]; % homozygous, unphased.
		het_unphased_color    = [0.66667 0.66667 0.66667]; % heterozygous.
		oddhet_unphased_color = [0.0     1.0     0.0    ]; % non-heterozygous data that isn't 100 hom.
	else
		% The 'project' is different than the 'hapmap'/'parent'.
		homolog_a_color       = [1.0 0.0 0.0];
		homolog_b_color       = [1.0 0.0 0.0];
		het_color             = [0.66667 0.66667 0.66667]; % heterozygous.
		hom_unphased_color    = [1.0     0.0     0.0    ]; % homozygous, unphased.
		het_unphased_color    = [0.66667 0.66667 0.66667]; % heterozygous.
		oddhet_unphased_color = [0.0     1.0     0.0    ]; % non-heterozygous data that isn't 100 hom.
	end;
end;


% phased data colors.
	% haploid colors.
	colorA          = homolog_a_color;
	colorB          = homolog_b_color;
	% diploid colors.
	colorAA         = homolog_a_color;
	colorAB         = het_color;
	colorBB         = homolog_b_color;
	% triploid colors.
	colorAAA        = homolog_a_color;
	colorAAB        = homolog_a_color*2/3 + homolog_b_color*1/3;
	colorABB        = homolog_a_color*1/3 + homolog_b_color*2/3;
	colorBBB        = homolog_b_color;
	% tetraploid colors.
	colorAAAA       = homolog_a_color;
	colorAAAB       = homolog_a_color*3/4 + homolog_b_color*1/4;
	colorAABB       = het_color;
	colorABBB       = homolog_a_color*1/4 + homolog_b_color*3/4;
	colorBBBB       = homolog_b_color;
	% pentaploid colors.
	colorAAAAA      = homolog_a_color;
	colorAAAAB      = homolog_a_color*4/5 + homolog_b_color*1/5;
	colorAAABB      = homolog_a_color*3/5 + homolog_b_color*2/5;
	colorAABBB      = homolog_a_color*2/5 + homolog_b_color*3/5;
	colorABBBB      = homolog_a_color*1/5 + homolog_b_color*4/5;
	colorBBBBB      = homolog_b_color;
	% hexaploid colors.
	colorAAAAAA     = homolog_a_color;
	colorAAAAAB     = homolog_a_color*5/6 + homolog_b_color*1/6;
	colorAAAABB     = homolog_a_color*4/6 + homolog_b_color*2/6;
	colorAAABBB     = het_color;
	colorAABBBB     = homolog_a_color*2/6 + homolog_b_color*4/6;
	colorABBBBB     = homolog_a_color*1/6 + homolog_b_color*5/6;
	colorBBBBBB     = homolog_b_color;
	% heptaploid colors.
	colorAAAAAAA    = homolog_a_color;
	colorAAAAAAB    = homolog_a_color*6/7 + homolog_b_color*1/7;
	colorAAAAABB    = homolog_a_color*5/7 + homolog_b_color*2/7;
	colorAAAABBB    = homolog_a_color*4/7 + homolog_b_color*3/7;
	colorAAABBBB    = homolog_a_color*3/7 + homolog_b_color*4/7;
	colorAABBBBB    = homolog_a_color*2/7 + homolog_b_color*5/7;
	colorABBBBBB    = homolog_a_color*1/7 + homolog_b_color*6/7;
	colorBBBBBBB    = homolog_b_color;
	% octaploid colors.
	colorAAAAAAAA   = homolog_a_color;
	colorAAAAAAAB   = homolog_a_color*7/8 + homolog_b_color*1/8;
	colorAAAAAABB   = homolog_a_color*6/8 + homolog_b_color*2/8;
	colorAAAAABBB   = homolog_a_color*5/8 + homolog_b_color*3/8;
	colorAAAABBBB   = het_color;
	colorAAABBBBB   = homolog_a_color*3/8 + homolog_b_color*5/8;
	colorAABBBBBB   = homolog_a_color*2/8 + homolog_b_color*6/8;
	colorABBBBBBB   = homolog_a_color*1/8 + homolog_b_color*7/8;
	colorBBBBBBBB   = homolog_b_color;
	% nonaploid colors.
	colorAAAAAAAAA  = homolog_a_color;
	colorAAAAAAAAB  = homolog_a_color*8/9 + homolog_b_color*1/9;
	colorAAAAAAABB  = homolog_a_color*7/9 + homolog_b_color*2/9;
	colorAAAAAABBB  = homolog_a_color*6/9 + homolog_b_color*3/9;
	colorAAAAABBBB  = homolog_a_color*5/9 + homolog_b_color*4/9;
	colorAAAABBBBB  = homolog_a_color*4/9 + homolog_b_color*5/9;
	colorAAABBBBBB  = homolog_a_color*3/9 + homolog_b_color*6/9;
	colorAABBBBBBB  = homolog_a_color*2/9 + homolog_b_color*7/9;
	colorABBBBBBBB  = homolog_a_color*1/9 + homolog_b_color*8/9;
	colorBBBBBBBBB  = homolog_b_color;

% unphased colors.
	% haploid colors.
	unphased_color_1of1 = hom_unphased_color;
	% diploid colors.
	unphased_color_2of2 = hom_unphased_color;
	unphased_color_1of2 = het_unphased_color;
	% triploid colors.
	unphased_color_3of3 = hom_unphased_color;
	unphased_color_2of3 = hom_unphased_color*2/3 + het_unphased_color*1/3;
	% tetraploid colors.
	unphased_color_4of4 = hom_unphased_color;
	unphased_color_3of4 = hom_unphased_color*3/4 + het_unphased_color*1/4;
	unphased_color_2of4 = het_unphased_color;
	% pentaploid colors.
	unphased_color_5of5 = hom_unphased_color;
	unphased_color_4of5 = hom_unphased_color*4/5 + het_unphased_color*1/5;
	unphased_color_3of5 = hom_unphased_color*3/5 + het_unphased_color*2/5;
	% hexaploid colors.
	unphased_color_6of6 = hom_unphased_color;
	unphased_color_5of6 = hom_unphased_color*5/6 + het_unphased_color*1/6;
	unphased_color_4of6 = hom_unphased_color*4/6 + het_unphased_color*2/6;
	unphased_color_3of6 = het_unphased_color;
	% heptaploid colors.
	unphased_color_7of7 = hom_unphased_color;
	unphased_color_6of7 = hom_unphased_color*6/7 + het_unphased_color*1/7;
	unphased_color_5of7 = hom_unphased_color*5/7 + het_unphased_color*2/7;
	unphased_color_4of7 = hom_unphased_color*4/7 + het_unphased_color*3/7;
	% octaploid colors.
	unphased_color_8of8 = hom_unphased_color;
	unphased_color_7of8 = hom_unphased_color*7/8 + het_unphased_color*1/8;
	unphased_color_6of8 = hom_unphased_color*6/8 + het_unphased_color*2/8;
	unphased_color_5of8 = hom_unphased_color*5/8 + het_unphased_color*3/8;
	unphased_color_4of8 = het_unphased_color;
	% nonaploid colors.
	unphased_color_9of9 = hom_unphased_color;
	unphased_color_8of9 = hom_unphased_color*8/9 + het_unphased_color*1/9;
	unphased_color_7of9 = hom_unphased_color*7/9 + het_unphased_color*2/9;
	unphased_color_6of9 = hom_unphased_color*6/9 + het_unphased_color*3/9;
	unphased_color_5of9 = hom_unphased_color*5/9 + het_unphased_color*4/9;


genomeDir  = [main_dir 'users/' genomeUser '/genomes/' genome '/'];

[centromeres, chr_sizes, figure_details, annotations, ploidy_default] = Load_genome_information(genomeDir);
[Aneuploidy]                                                          = Load_dataset_information(projectDir);
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

%define colors for colorBars plot
colorNoData = [1.0   1.0   1.0  ]; %used when no data is available for the bin.
colorInit   = [0.5   0.5   0.5  ]; %external; used in blending at ends of chr.
colorHET    = [0.0   0.0   0.0  ]; % near 1:1 ratio SNPs
colorOddHET = [0.0   1.0   0.0  ]; % Het, but not near 1:1 ratio SNPs.
colorHOM    = [1.0   0.0   0.0  ]; % Hom SNPs;

colorAB     = [0.667 0.667 0.667]; % heterozygous.
colorA      = [1.0   0.0   1.0  ]; % homozygous a:magenta.
colorB      = [0.0   1.0   1.0  ]; % homozygous b:cyan.

fprintf(['\nGenerating LOH-map figure from ''' project ''' vs. (hapmap)''' hapmap ''' data.\n']);

% Initializes vectors used to hold allelic ratios for each chromosome segment.
if (useHapmap)
	new_bases_per_bin = bases_per_bin;
else
	new_bases_per_bin = bases_per_bin/2;
end;
for chr = 1:length(chr_sizes)
	%   1 : experimental : phased ratio data.
	%   2 : experimental : unphased ratio data.
	%   3 : reference    : phased ratio data.
	%   4 : reference    : unphased ratio data.
	chr_length = ceil(chr_size(chr)/new_bases_per_bin);
	for j = 1:4
		chr_SNPdata{chr,j} = ones(chr_length,1);
	end;
	% Colors used to illustrate SNP/LOH data.
	%    chr_SNPdata_colorsC           : colors scheme defined by hapmap or red for unspecified LOH.
	%    chr_SNPdata_colorsC_alternate : color scheme defined to accentuate difference between homozygous and skewed heterozygous data.
	for j = 1:3
		% Track the RGB value sum per standard bin, then divide by the count to reach the average color per standard genome bin.
		chr_SNPdata_colorsC{chr,j}           = zeros(chr_length,1);   
		chr_SNPdata_colorsP{chr,j}           = zeros(chr_length,1);
		chr_SNPdata_colorsC_alternate{chr,j} = zeros(chr_length,1);
	end;

	% Track the number of SNP colors per standard bin.
	chr_SNPdata_countC{chr} = zeros(chr_length,1);
	chr_SNPdata_countP{chr} = zeros(chr_length,1); 
end;


%%================================================================================================
% Load SNP/LOH data.
%-------------------------------------------------------------------------------------------------
if (useHapmap)
%
% Only run when compared vs. a hapmap.
%
	if (exist([projectDir 'SNP_' SNP_verString '.all3.mat'],'file') == 0)
		fprintf('\nAllelic fraction MAT file 3 not found, generating.\n');
		process_2dataset_hapmap_allelicRatios(projectDir, projectDir, hapmapDir, chr_size, chr_name, chr_in_use, SNP_verString);
	else
		fprintf('\nAllelic fraction MAT file 3 found, loading.\n');
	end;
	load([projectDir 'SNP_' SNP_verString '.all3.mat']);
	% child data:  'C_chr_SNP_data_positions','C_chr_SNP_data_ratios','C_chr_count','C_chr_baseCall','C_chr_SNP_homologA','C_chr_SNP_homologB','C_chr_SNP_flipHomologs'
	% parent data: 'P_chr_SNP_data_positions','P_chr_SNP_data_ratios','P_chr_count','P_chr_baseCall','P_chr_SNP_homologA','P_chr_SNP_homologB','P_chr_SNP_flipHomologs'
	%
	% C_chr_SNP_data_positions = coordinate of SNP.
	% C_chr_SNP_data_ratios    = allelic ratio of SNP.
	% C_chr_count              = number of reads at SNP coordinate.
	% C_chr_baseCall           = majority basecall of SNP.
	% C_chr_SNP_homologA       = hapmap homolog a basecall.
	% C_chr_SNP_homologB       = hapmap homolog b basecall.
	% C_chr_SNP_flipHomologs   = does hapmap entry need flipped?
else
%
% Run when compared vs. a parent dataset or vs. itself. 
%
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
% Load 'Common_CNV.mat' file containing CNV estimates per standard genome bin.
%-------------------------------------------------------------------------------------------------
fprintf('\nLoading "Common_CNV" data file for ddRADseq project.');
load([projectDir 'Common_CNV.mat']);   % 'CNVplot2', 'genome_CNV'
[chr_breaks, chrCopyNum, ploidyAdjust] = FindChrSizes_4(Aneuploidy,CNVplot2,ploidy,num_chrs,chr_in_use);


%%================================================================================================
% Save workspace variables for use in "allelic_ratios_ddRADseq_D.m"
%-------------------------------------------------------------------------------------------------
save([projectDir 'allelic_ratios_ddRADseq_B.workspace_variables.mat']);


%%================================================================================================
% Process SNP/hapmap data to determine colors for presentation.
%-------------------------------------------------------------------------------------------------
if (useHapmap)
%
% Only run when compared vs. a hapmap.
%
	for chr = 1:num_chrs
		if (chr_in_use(chr) == 1)
			if (length(C_chr_count{chr}) > 1)
				%
				% Determining colors for each SNP coordinate.
				%
				for i = 1:length(C_chr_count{chr})
					pos                             = ceil(C_chr_SNP_data_positions{chr}(i)/new_bases_per_bin);
					localCopyEstimate               = round(CNVplot2{chr}(pos)*ploidy*ploidyAdjust);
					chr_SNPdata{chr,2}(pos)         = C_chr_SNP_data_ratios{ chr}(i);
					baseCall                        = C_chr_baseCall{        chr}{i};
					homologA                        = C_chr_SNP_homologA{    chr}{i};
					homologB                        = C_chr_SNP_homologB{    chr}{i};
					flipper                         = C_chr_SNP_flipHomologs{chr}(i);
					if (flipper)
						temp                    = homologA;
						homologA                = homologB;
						homologB                = temp;
					end;
					allelicFraction                 = C_chr_SNP_data_ratios{chr}(i);
					if (localCopyEstimate <= 0)
						colorList = [1 1 1];
					elseif (localCopyEstimate == 1)
						if (baseCall == homologA)
							colorList = colorA;
						elseif (baseCall == homologB)
							colorList = colorB;
						else
							colorList = hom_unphased_color;
						end;
					elseif (localCopyEstimate == 2)
						if (baseCall == homologA)
							if (allelicFraction > 3/4)
								colorList = colorAA;
							else
								colorList = colorAB;
							end;
						elseif (baseCall == homologB)
							if (allelicFraction > 3/4)
								colorList = colorBB;
							else
								colorList = colorAB;
							end;
						else
							if (allelicFraction > 3/4)
								colorList = hom_unphased_color;
							else
								colorList = het_unphased_color;
							end;
						end;
					elseif (localCopyEstimate == 3)
						if (baseCall == homologA)
							if (allelicFraction > 5/6)
								colorList = colorAAA;
							else
								colorList = colorAAB;
							end;
						elseif (baseCall == homologB)
							if (allelicFraction > 5/6)
								colorList = colorBBB;
							else
								colorList = colorABB;
							end;
						else
							if (allelicFraction > 5/6)
								colorList = unphased_color_3of3;
							else
								colorList = unphased_color_2of3;
							end;
						end;
					elseif (localCopyEstimate == 4)
						if (baseCall == homologA)
							if (allelicFraction > 7/8)
								colorList = colorAAAA;
							elseif (allelicFraction > 5/8)
								colorList = colorAAAB;
							else
								colorList = colorAABB;
							end;
						elseif (baseCall == homologB)
							if (allelicFraction > 7/8)
								colorList = colorBBBB;
							elseif (allelicFraction > 5/8)
								colorList = colorABBB;
							else
								colorList = colorAABB;
							end;
						else
							if (allelicFraction > 7/8)
								colorList = unphased_color_4of4;
							elseif (allelicFraction > 5/8)
								colorList = unphased_color_3of4;
							else
								colorList = unphased_color_2of4;
							end;
						end;
					elseif (localCopyEstimate == 5)
						if (baseCall == homologA)
							if (allelicFraction > 9/10)
								colorList = colorAAAAA;
							elseif (allelicFraction > 7/10)
								colorList = colorAAAAB;
							else
								colorList = colorAAABB;
							end;
						elseif (baseCall == homologB)
							if (allelicFraction > 9/10)
								colorList = colorBBBBB;
							elseif (allelicFraction > 7/10)
								colorList = colorABBBB;
							else
								colorList = colorAABBB;
							end;
						else
							if (allelicFraction > 9/10)
								colorList = unphased_color_5of5;
							elseif (allelicFraction > 7/10)
								colorList = unphased_color_4of5;
							else
								colorList = unphased_color_3of5;
							end;
						end;
					elseif (localCopyEstimate == 6)
						if (baseCall == homologA)
							if (allelicFraction > 11/12)
								colorList = colorAAAAAA;
							elseif (allelicFraction > 9/12)
								colorList = colorAAAAAB;
							elseif (allelicFraction > 7/12)
								colorList = colorAAAABB;
							else
								colorList = colorAAABBB;
							end;
						elseif (baseCall == homologB)
							if (allelicFraction > 11/12)
								colorList = colorBBBBBB;
							elseif (allelicFraction > 9/12)
								colorList = colorABBBBB;
							elseif (allelicFraction > 7/12)
								colorList = colorAABBBB;
							else
								colorList = colorAAABBB;
							end;
						else
							if (allelicFraction > 11/12)
								colorList = unphased_color_6of6;
							elseif (allelicFraction > 9/12)
								colorList = unphased_color_5of6;
							elseif (allelicFraction > 7/12)
								colorList = unphased_color_4of6;
							else
								colorList = unphased_color_3of6;
							end;
						end;
					elseif (localCopyEstimate == 7)
						if (baseCall == homologA)
							if (allelicFraction > 13/14)
								colorList = colorAAAAAAA;
							elseif (allelicFraction > 11/14)
								colorList = colorAAAAAAB;
							elseif (allelicFraction > 9/14)
								colorList = colorAAAAABB;
							else
								colorList = colorAAAABBB;
							end;
						elseif (baseCall == homologB)
							if (allelicFraction > 13/14)
								colorList = colorBBBBBBB;
							elseif (allelicFraction > 11/14)
								colorList = colorABBBBBB;
							elseif (allelicFraction > 9/14)
								colorList = colorAABBBBB;
							else
								colorList = colorAAABBBB;
							end;
						else
							if (allelicFraction > 13/14)
								colorList = unphased_color_7of7;
							elseif (allelicFraction > 11/14)
								colorList = unphased_color_6of7;
							elseif (allelicFraction > 9/14)
								colorList = unphased_color_5of7;
							else
								colorList = unphased_color_4of7;
							end;
						end;
					elseif (localCopyEstimate == 8)
						if (baseCall == homologA)
							if (allelicFraction > 15/16)
								colorList = colorAAAAAAAA;
							elseif (allelicFraction > 13/16)
								colorList = colorAAAAAAAB;
							elseif (allelicFraction > 11/16)
								colorList = colorAAAAAABB;
							elseif (allelicFraction > 9/16)
								colorList = colorAAAAABBB;
							else
								colorList = colorAAAABBBB;
							end;
						elseif (baseCall == homologB)
							if (allelicFraction > 15/16)
								colorList = colorBBBBBBBB;
							elseif (allelicFraction > 13/16)
								colorList = colorABBBBBBB;
							elseif (allelicFraction > 11/16)
								colorList = colorAABBBBBB;
							elseif (allelicFraction > 9/16)
								colorList = colorAAABBBBB;
							else
								colorList = colorAAAABBBB;
							end;
						else
							if (allelicFraction > 15/16)
								colorList = unphased_color_8of8;
							elseif (allelicFraction > 13/16)
								colorList = unphased_color_7of8;
							elseif (allelicFraction > 11/16)
								colorList = unphased_color_6of8;
							elseif (allelicFraction > 9/16)
								colorList = unphased_color_5of8;
							else
								colorList = unphased_color_4of8;
							end;
						end;
					elseif (localCopyEstimate >= 9)
						if (baseCall == homologA)
							if (allelicFraction > 17/18)
								colorList = colorAAAAAAAAA;
							elseif (allelicFraction > 15/18)
								colorList = colorAAAAAAAAB;
							elseif (allelicFraction > 13/18)
								colorList = colorAAAAAAABB;
							elseif (allelicFraction > 11/18)
								colorList = colorAAAAAABBB;
							else
								colorList = colorAAAAABBBB;
							end;
						elseif (baseCall == homologB)
							if (allelicFraction > 17/18)
								colorList = colorBBBBBBBBB;
							elseif (allelicFraction > 15/18)
								colorList = colorABBBBBBBB;
							elseif (allelicFraction > 13/18)
								colorList = colorAABBBBBBB;
							elseif (allelicFraction > 11/18)
								colorList = colorAAABBBBBB;
							else
								colorList = colorAAAABBBBB;
							end;
						else
							if (allelicFraction > 17/18)
								colorList = unphased_color_9of9;
							elseif (allelicFraction > 15/18)
								colorList = unphased_color_8of9;
							elseif (allelicFraction > 13/18)
								colorList = unphased_color_7of9;
							elseif (allelicFraction > 11/18)
								colorList = unphased_color_6of9;
							else
								colorList = unphased_color_5of9;
							end;
						end;
					end;

					chr_SNPdata_colorsC{chr,1}(pos) = chr_SNPdata_colorsC{chr,1}(pos) + colorList(1);
					chr_SNPdata_colorsC{chr,2}(pos) = chr_SNPdata_colorsC{chr,2}(pos) + colorList(2);
					chr_SNPdata_colorsC{chr,3}(pos) = chr_SNPdata_colorsC{chr,3}(pos) + colorList(3);
					chr_SNPdata_countC{ chr  }(pos) = chr_SNPdata_countC{ chr  }(pos) + 1;
				end;

				%
				% Average color per bin.
				%
				for pos = 1:length(chr_SNPdata_countC{chr})
					if (chr_SNPdata_countC{chr}(pos) > 0)
						chr_SNPdata_colorsC{chr,1}(pos) = chr_SNPdata_colorsC{chr,1}(pos)/chr_SNPdata_countC{chr}(pos);
						chr_SNPdata_colorsC{chr,2}(pos) = chr_SNPdata_colorsC{chr,2}(pos)/chr_SNPdata_countC{chr}(pos);
						chr_SNPdata_colorsC{chr,3}(pos) = chr_SNPdata_colorsC{chr,3}(pos)/chr_SNPdata_countC{chr}(pos);
					else
						chr_SNPdata_colorsC{chr,1}(pos) = 1.0;
						chr_SNPdata_colorsC{chr,2}(pos) = 1.0;
						chr_SNPdata_colorsC{chr,3}(pos) = 1.0;
					end;
				end;
			end;
		end;
	end;
else
%
% Run when compared vs. a parent dataset or vs. itself.
%
	for chr = 1:num_chrs
		if (chr_in_use(chr) == 1)
			if (length(C_chr_count{chr}) > 1)
				for i = 1:length(C_chr_count{chr})
					pos = ceil(C_chr_SNP_data_positions{chr}(i)/new_bases_per_bin);
					if (C_chr_SNP_data_ratios{chr}(i) < chr_SNPdata{chr,2}(pos))
						chr_SNPdata{chr,2}(pos)      = C_chr_SNP_data_ratios{chr}(i);
						colorList                    = [1.0 1.0 1.0];
						chr_SNPdata_colorsC{chr,1}(i) = colorList(1);
						chr_SNPdata_colorsC{chr,2}(i) = colorList(2);
						chr_SNPdata_colorsC{chr,3}(i) = colorList(3);
					end;
				end;
			end;
			if (length(P_chr_count{chr}) > 1)
				for i = 1:length(P_chr_count{chr})
					pos = ceil(P_chr_SNP_data_positions{chr}(i)/new_bases_per_bin);
					if (P_chr_SNP_data_ratios{chr}(i) < chr_SNPdata{chr,4}(pos))
						chr_SNPdata{chr,4}(pos)      = P_chr_SNP_data_ratios{chr}(i);
						colorList                    = [1.0 1.0 1.0];
						chr_SNPdata_colorsP{chr,1}(i) = colorList(1);
						chr_SNPdata_colorsP{chr,2}(i) = colorList(2);
						chr_SNPdata_colorsP{chr,3}(i) = colorList(3);
					end;
				end;
			end;
		end;
	end;
end;

save([projectDir 'SNP_' SNP_verString '.reduced.mat'],'chr_SNPdata','new_bases_per_bin','chr_SNPdata_colorsC','chr_SNPdata_colorsP');



%%================================================================================================
% Make histogram of 'homozygous' data.
%-------------------------------------------------------------------------------------------------
histogram_fig = figure();
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		data_C = [];
		data_P = [];
		for i = 1:length(chr_SNPdata{chr,2})
			data_C = [data_C chr_SNPdata{chr,2}(i)];
		end;
		for i = 1:length(chr_SNPdata{chr,4})
			data_P = [data_P chr_SNPdata{chr,4}(i)];
		end;
		histogram_C = hist(data_C,100);	
		histogram_P = hist(data_P,100);
		final_histogram_C = [fliplr(histogram_C) histogram_C];
		final_histogram_P = [fliplr(histogram_P) histogram_P];

		subplot(2,num_chrs,chr);
		hold on;
		plot(1:200, log(final_histogram_C+1),'Color',[1.0 0.0 0.0]);
		plot(1:200, log(final_histogram_P+1),'Color',[1/3 1/3 1/3]);
		ylim([0 6]);
		title([chr_label{chr}]);
		set(gca,'XTick',[0 50 100 150 200]);
		set(gca,'XTickLabel',{'0','1/4','1/2','3/4','1'});
		ylabel('log(data count)');
		xlabel('allelic ratio');
		hold off;

		subplot(2,num_chrs,chr+num_chrs);
		hold on;
		plot(1:200, final_histogram_C,'Color',[1.0 0.0 0.0]);
		plot(1:200, final_histogram_P,'Color',[1/3 1/3 1/3]);
		ylim([0 200]);
		title([chr_label{chr}]);
		set(gca,'XTick',[0 50 100 150 200]);
		set(gca,'XTickLabel',{'0','1/4','1/2','3/4','1'});
		ylabel('data count');
		xlabel('allelic ratio');
		hold off;
	end;
end;
set(histogram_fig,'PaperPosition',[0 0 8 1.5]*4);
saveas(histogram_fig, [projectDir 'fig.allelic_fraction_histogram.eps'], 'epsc');
saveas(histogram_fig, [projectDir 'fig.allelic_fraction_histogram.png'], 'png');
delete(histogram_fig);


%%================================================================================================
% Setup for main figure generation.
%-------------------------------------------------------------------------------------------------
fig = figure(1);
set(gcf, 'Position', [0 70 1024 600]);
largestChr = find(chr_width == max(chr_width));


%%================================================================================================
% Setup for linear-view figure generation.
%-------------------------------------------------------------------------------------------------
if (Linear_display == true)
	Linear_fig = figure(2);
	Linear_genome_size   = sum(chr_size);

	Linear_Chr_max_width = 0.91;               % width for all chromosomes across figure.  1.00 - leftMargin - rightMargin - subfigure gaps.
	Linear_left_start    = 0.02;               % left margin (also right margin).
	Linear_left_chr_gap  = 0.07/(num_chrs-1);  % gaps between chr subfigures.

	Linear_height        = 0.6;
	Linear_base          = 0.1;
	Linear_TickSize      = -0.01;  %negative for outside, percentage of longest chr figure.
	maxY                 = 1; % ploidyBase*2;
	Linear_left          = Linear_left_start;

	axisLabelPosition_horiz = -50000/bases_per_bin;
	axisLabelPosition_horiz = 0.01125;
end;

axisLabelPosition_vert = -50000/bases_per_bin;
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
		text(-50000/5000/2*3, maxY/2,     chr_label{chr}, 'Rotation',90, 'HorizontalAlignment','center', 'VerticalAlign','bottom', 'Fontsize',20);
		set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
		set(gca,'XTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8','2.0','2.2','2.4','2.6','2.8','3.0','3.2'});

		% standard : This section sets the Y-axis labelling.
		axisLabelPosition = -50000/bases_per_bin;
		set(gca,'FontSize',12);
		if (chr == find(chr_posY == max(chr_posY)))
			title([ project ' allelic fraction map'],'Interpreter','none','FontSize',24);
		end;
		% standard : end axes labels etc.

		% standard : draw colorbars.
		if (useHapmap)   % an experimental dataset vs. a hapmap.
			for i = 1:ceil(chr_size(chr)/new_bases_per_bin)
				colorR   = chr_SNPdata_colorsC{chr,1}(i);
				colorG   = chr_SNPdata_colorsC{chr,2}(i);
				colorB   = chr_SNPdata_colorsC{chr,3}(i);
				if (colorR < 1) || (colorG < 1) || (colorB < 1)
					plot([i i], [0 maxY],'Color',[colorR colorG colorB]);
				end;
			end;
		% Following variation will draw experimetnal vs. reference dataset like the above vs. hapmap figure.
		%	elseif (useParent)   % an experimental dataset vs. a reference dataset.
		%		for i = 1:ceil(chr_size(chr)/new_bases_per_bin)
		%			colorR   = chr_SNPdata_colorsP{chr,1}(i);
		%			colorG   = chr_SNPdata_colorsP{chr,2}(i);
		%			colorB   = chr_SNPdata_colorsP{chr,3}(i);
		%			if (colorR < 1) || (colorG < 1) || (colorB < 1)
		%				plot([i i], [0 maxY],'Color',[colorR colorG colorB]);
		%			end;
		%		end;
		elseif (useParent)   % an experimental dataset vs. a reference dataset.
			for i = 1:ceil(chr_size(chr)/new_bases_per_bin)
				datumY_C = chr_SNPdata{chr,2}(i)*maxY;
				datumY_P = chr_SNPdata{chr,4}(i)*maxY;
				plot([i/2 i/2], [maxY datumY_C     ],'Color',[1.0 0.0 0.0]);
				plot([i/2 i/2], [0    maxY-datumY_P],'Color',[1/3 1/3 1/3]);
			end;
		else   % only a reference dataset.
			for i = 1:ceil(chr_size(chr)/new_bases_per_bin)
				datumY_P = chr_SNPdata{chr,4}(i)*maxY;
				plot([i/2 i/2], [maxY datumY_P     ],'Color',[1/3 1/3 1/3]);
				plot([i/2 i/2], [0    maxY-datumY_P],'Color',[1/3 1/3 1/3]);
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
		end;
		% standard : end show annotation locations.
		hold off;

		%% Linear figure draw section
		if (Linear_display == true)
			figure(Linear_fig);
			Linear_width = Linear_Chr_max_width*chr_size(chr)/Linear_genome_size;
			subplot('Position',[Linear_left Linear_base Linear_width Linear_height]);
			hold on;
			Linear_left = Linear_left + Linear_width + Linear_left_chr_gap;
			title(chr_label{chr},'Interpreter','none','FontSize',20);

			% linear : draw colorbars
			if (useHapmap)   % an experimental dataset vs. a hapmap.
				for i = 1:ceil(chr_size(chr)/new_bases_per_bin)
					colorR   = chr_SNPdata_colorsC{chr,1}(i);
					colorG   = chr_SNPdata_colorsC{chr,2}(i);
					colorB   = chr_SNPdata_colorsC{chr,3}(i);
					if (colorR < 1) || (colorG < 1) || (colorB < 1)
						plot([i i], [0 maxY],'Color',[colorR colorG colorB]);
					end;
				end;
			% Following variation will draw experimetnal vs. reference dataset like the above vs. hapmap figure.
			%	elseif (useParent)   % an experimental dataset vs. a reference dataset.
			%		for i = 1:ceil(chr_size(chr)/new_bases_per_bin)
			%			colorR   = chr_SNPdata_colorsC{chr,1}(i);
			%			colorG   = chr_SNPdata_colorsC{chr,2}(i);
			%			colorB   = chr_SNPdata_colorsC{chr,3}(i);
			%			if (colorR < 1) || (colorG < 1) || (colorB < 1)
			%				plot([i i], [0 maxY],'Color',[colorR colorG colorB]);
			%			end;
			%		end;
			elseif (useParent)   % an experimental dataset vs. a reference dataset.
				for i = 1:ceil(chr_size(chr)/new_bases_per_bin)
					datumY_C = chr_SNPdata{chr,2}(i)*maxY;
					datumY_P = chr_SNPdata{chr,4}(i)*maxY;
					plot([i/2 i/2], [maxY datumY_C     ],'Color',[1.0 0.0 0.0]);
					plot([i/2 i/2], [0    maxY-datumY_P],'Color',[1/3 1/3 1/3]);
				end;
			else   % only a reference dataset.
				for i = 1:ceil(chr_size(chr)/new_bases_per_bin)
					datumY_C = chr_SNPdata{chr,2}(i)*maxY;
					datumY_P = chr_SNPdata{chr,4}(i)*maxY;
					plot([i/2 i/2], [maxY datumY_P     ],'Color',[1/3 1/3 1/3]);
					plot([i/2 i/2], [0    maxY-datumY_P],'Color',[1/3 1/3 1/3]);
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
			set(gca,'FontSize',12);
			% linear : end final reformatting.

			hold off;
	        
			% shift back to main figure generation.
			figure(fig);
			first_chr = false;
		end;
	end;
end;

%% Save figures.
set(fig,'PaperPosition',[0 0 8 6]*2);
saveas(fig,        [projectDir 'fig.allelic_ratio-map.c1.eps'], 'epsc');
saveas(fig,        [projectDir 'fig.allelic_ratio-map.c1.png'], 'png');
delete(fig);

set(Linear_fig,'PaperPosition',[0 0 8 0.62222222]*2);
saveas(Linear_fig, [projectDir 'fig.allelic_ratio-map.c2.eps'], 'epsc');
saveas(Linear_fig, [projectDir 'fig.allelic_ratio-map.c2.png'], 'png');
delete(Linear_fig);

%%================================================================================================
% end stuff
%=================================================================================================
end

