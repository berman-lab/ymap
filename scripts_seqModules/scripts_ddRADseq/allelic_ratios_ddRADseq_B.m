function [] = allelic_ratios_ddRADseq_B(main_dir,user,genomeUser,project,parent,hapmap,genome,ploidyEstimateString,ploidyBaseString,SNP_verString,LOH_verString,CNV_verString,displayBREAKS);
addpath('../');

workingDir      = [main_dir 'users/' user '/projects/' project '/'];

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


%% =========================================================================================
% Define colors for figure generation.
%-------------------------------------------------------------------------------------------
phased_and_unphased_color_definitions;


genomeDir  = [main_dir 'users/' genomeUser '/genomes/' genome '/'];
[centromeres, chr_sizes, figure_details, annotations, ploidy_default] = Load_genome_information(genomeDir);
[Aneuploidy]                                                          = Load_dataset_information(projectDir);
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
		chr_label {figure_details(detailID).chr} = figure_details(detailID).label;
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

fprintf(['\nGenerating LOH-map figure from ''' project ''' vs. (hapmap)''' hapmap ''' data.\n']);
% Initializes vectors used to hold allelic ratios for each chromosome segment.
if (useHapmap)
	new_bases_per_bin = bases_per_bin;
else
	new_bases_per_bin = bases_per_bin/2;
end;
for chr = 1:length(chr_sizes)
	%   1 : phased SNP ratio data.
	%   2 : unphased SNP ratio data.
	%   3 : phased SNP position data.
	%   4 : unphased SNP position data.
	chr_length = ceil(chr_size(chr)/new_bases_per_bin);
	for j = 1:4
		for i = 1:chr_length
			chr_SNPdata{chr,j}{i} = [];
		end;
	end;
	% Colors used to illustrate SNP/LOH data.
	%    chr_SNPdata_colorsC           : colors scheme defined by hapmap or red for unspecified LOH.
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
	%
	% C_chr_SNP_data_positions = coordinate of SNP.
	% C_chr_SNP_data_ratios    = allelic ratio of SNP.
	% C_chr_count              = number of reads at SNP coordinate.
	% C_chr_baseCall           = majority basecall of SNP.
	% C_chr_SNP_homologA       = hapmap homolog a basecall.
	% C_chr_SNP_homologB       = hapmap homolog b basecall.
	% C_chr_SNP_flipHomologs   = does hapmap entry need flipped?   0 = 'correct phase, no', 1 = 'incorrect phase, yes', 10 = 'no phasing info'.
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
	%%%% chr_SNPdata{chr,1}(SNP) = phased SNP ratio data.
	%%%% chr_SNPdata{chr,2}(SNP) = unphased SNP ratio data.
	%%%% chr_SNPdata{chr,3}(SNP) = phased SNP position data.
	%%%% chr_SNPdata{chr,4}(SNP) = unphased SNP position data.
	for chr = 1:num_chrs
		if (chr_in_use(chr) == 1)
			if (length(C_chr_count{chr}) > 1)
				%
				% Old method for determining colors for SNP coordinates.  Building SNP data sctructure 'chr_SNPdata'.
				%
				for SNP = 1:length(C_chr_SNP_data_positions{chr})  % 'length(C_chr_SNP_data_positions)' is the number of SNPs per chromosome.
					coordinate                      = C_chr_SNP_data_positions{chr}(SNP);
					pos                             = ceil(coordinate/new_bases_per_bin);
					localCopyEstimate               = round(CNVplot2{chr}(pos)*ploidy*ploidyAdjust);
					baseCall                        = C_chr_baseCall{        chr}{SNP};
					homologA                        = C_chr_SNP_homologA{    chr}{SNP};
					homologB                        = C_chr_SNP_homologB{    chr}{SNP};
					flipper                         = C_chr_SNP_flipHomologs{chr}(SNP);
					allelic_ratio                   = C_chr_SNP_data_ratios{ chr}(SNP);
					if (flipper == 10)                         % Variable 'flipper' value of '10' indicates no phasing information is available in the hapmap.
						baseCall                = 'Z';     % Variable 'baseCall' value of 'Z' will prevent either hapmap allele from matching and so unphased ratio colors will be used in the following section.
						chr_SNPdata{chr,2}{pos} = [chr_SNPdata{chr,2}{pos} C_chr_SNP_data_ratios{chr}(SNP) 1-C_chr_SNP_data_ratios{chr}(SNP)];
						chr_SNPdata{chr,4}{pos} = [chr_SNPdata{chr,4}{pos} coordinate coordinate];
					elseif (flipper == 1)
						temp                    = homologA;
						homologA                = homologB;
						homologB                = temp;
						chr_SNPdata{chr,1}{pos} = [chr_SNPdata{chr,1}{pos} 1-allelic_ratio];
						chr_SNPdata{chr,3}{pos} = [chr_SNPdata{chr,3}{pos} coordinate];
					else % (flipper == 0)
						chr_SNPdata{chr,1}{pos} = [chr_SNPdata{chr,1}{pos} allelic_ratio];
						chr_SNPdata{chr,3}{pos} = [chr_SNPdata{chr,3}{pos} coordinate];
					end;
                    
					allelicFraction                 = C_chr_SNP_data_ratios{chr}(SNP);
					if (localCopyEstimate <= 0);                        colorList = colorNoData;
					elseif (localCopyEstimate == 1)
						if (baseCall == homologB);                  colorList = colorA;
						elseif (baseCall == homologA);              colorList = colorB;
						else                                        colorList = unphased_color_1of1;
						end;
					elseif (localCopyEstimate == 2)
						if (baseCall == homologB)
							if (allelicFraction > 3/4);         colorList = colorAA;
							else                                colorList = colorAB;
							end;
						elseif (baseCall == homologA)
							if (allelicFraction > 3/4);         colorList = colorBB;
							else                                colorList = colorAB;
							end;
						else
							if (allelicFraction > 3/4);         colorList = unphased_color_2of2;
							else                                colorList = unphased_color_1of2;
							end;
						end;
					elseif (localCopyEstimate == 3)
						if (baseCall == homologB)
							if (allelicFraction > 5/6);         colorList = colorAAA;
							else                                colorList = colorAAB;
							end;
						elseif (baseCall == homologA)
							if (allelicFraction > 5/6);         colorList = colorBBB;
							else                                colorList = colorABB;
							end;
						else
							if (allelicFraction > 5/6);         colorList = unphased_color_3of3;
							else                                colorList = unphased_color_2of3;
							end;
						end;
					elseif (localCopyEstimate == 4)
						if (baseCall == homologB)
							if (allelicFraction > 7/8);         colorList = colorAAAA;
							elseif (allelicFraction > 5/8);     colorList = colorAAAB;
							else                                colorList = colorAABB;
							end;
						elseif (baseCall == homologA)
							if (allelicFraction > 7/8);         colorList = colorBBBB;
							elseif (allelicFraction > 5/8);     colorList = colorABBB;
							else                                colorList = colorAABB;
							end;
						else
							if (allelicFraction > 7/8);         colorList = unphased_color_4of4;
							elseif (allelicFraction > 5/8);     colorList = unphased_color_3of4;
							else                                colorList = unphased_color_2of4;
							end;
						end;
					elseif (localCopyEstimate == 5)
						if (baseCall == homologB)
							if (allelicFraction > 9/10);        colorList = colorAAAAA;
							elseif (allelicFraction > 7/10);    colorList = colorAAAAB;
							else                                colorList = colorAAABB;
							end;
						elseif (baseCall == homologA)
							if (allelicFraction > 9/10);        colorList = colorBBBBB;
							elseif (allelicFraction > 7/10);    colorList = colorABBBB;
							else                                colorList = colorAABBB;
							end;
						else
							if (allelicFraction > 9/10);        colorList = unphased_color_5of5;
							elseif (allelicFraction > 7/10);    colorList = unphased_color_4of5;
							else                                colorList = unphased_color_3of5;
							end;
						end;
					elseif (localCopyEstimate == 6)
						if (baseCall == homologB)
							if (allelicFraction > 11/12);       colorList = colorAAAAAA;
							elseif (allelicFraction > 9/12);    colorList = colorAAAAAB;
							elseif (allelicFraction > 7/12);    colorList = colorAAAABB;
							else                                colorList = colorAAABBB;
							end;
						elseif (baseCall == homologA)
							if (allelicFraction > 11/12);       colorList = colorBBBBBB;
							elseif (allelicFraction > 9/12);    colorList = colorABBBBB;
							elseif (allelicFraction > 7/12);    colorList = colorAABBBB;
							else                                colorList = colorAAABBB;
							end;
						else
							if (allelicFraction > 11/12);       colorList = unphased_color_6of6;
							elseif (allelicFraction > 9/12);    colorList = unphased_color_5of6;
							elseif (allelicFraction > 7/12);    colorList = unphased_color_4of6;
							else                                colorList = unphased_color_3of6;
							end;
						end;
					elseif (localCopyEstimate == 7)
						if (baseCall == homologB)
							if (allelicFraction > 13/14);       colorList = colorAAAAAAA;
							elseif (allelicFraction > 11/14);   colorList = colorAAAAAAB;
							elseif (allelicFraction > 9/14);    colorList = colorAAAAABB;
							else                                colorList = colorAAAABBB;
							end;
						elseif (baseCall == homologA)
							if (allelicFraction > 13/14);       colorList = colorBBBBBBB;
							elseif (allelicFraction > 11/14);   colorList = colorABBBBBB;
							elseif (allelicFraction > 9/14);    colorList = colorAABBBBB;
							else                                colorList = colorAAABBBB;
							end;
						else
							if (allelicFraction > 13/14);       colorList = unphased_color_7of7;
							elseif (allelicFraction > 11/14);   colorList = unphased_color_6of7;
							elseif (allelicFraction > 9/14);    colorList = unphased_color_5of7;
							else                                colorList = unphased_color_4of7;
							end;
						end;
					elseif (localCopyEstimate == 8)
						if (baseCall == homologB)
							if (allelicFraction > 15/16);       colorList = colorAAAAAAAA;
							elseif (allelicFraction > 13/16);   colorList = colorAAAAAAAB;
							elseif (allelicFraction > 11/16);   colorList = colorAAAAAABB;
							elseif (allelicFraction > 9/16);    colorList = colorAAAAABBB;
							else                                colorList = colorAAAABBBB;
							end;
						elseif (baseCall == homologA)
							if (allelicFraction > 15/16);       colorList = colorBBBBBBBB;
							elseif (allelicFraction > 13/16);   colorList = colorABBBBBBB;
							elseif (allelicFraction > 11/16);   colorList = colorAABBBBBB;
							elseif (allelicFraction > 9/16);    colorList = colorAAABBBBB;
							else                                colorList = colorAAAABBBB;
							end;
						else
							if (allelicFraction > 15/16);       colorList = unphased_color_8of8;
							elseif (allelicFraction > 13/16);   colorList = unphased_color_7of8;
							elseif (allelicFraction > 11/16);   colorList = unphased_color_6of8;
							elseif (allelicFraction > 9/16);    colorList = unphased_color_5of8;
							else                                colorList = unphased_color_4of8;
							end;
						end;
					elseif (localCopyEstimate >= 9)
						if (baseCall == homologB)
							if (allelicFraction > 17/18);       colorList = colorAAAAAAAAA;
							elseif (allelicFraction > 15/18);   colorList = colorAAAAAAAAB;
							elseif (allelicFraction > 13/18);   colorList = colorAAAAAAABB;
							elseif (allelicFraction > 11/18);   colorList = colorAAAAAABBB;
							else                                colorList = colorAAAAABBBB;
							end;
						elseif (baseCall == homologA)
							if (allelicFraction > 17/18);       colorList = colorBBBBBBBBB;
							elseif (allelicFraction > 15/18);   colorList = colorABBBBBBBB;
							elseif (allelicFraction > 13/18);   colorList = colorAABBBBBBB;
							elseif (allelicFraction > 11/18);   colorList = colorAAABBBBBB;
							else                                colorList = colorAAAABBBBB;
							end;
						else
							if (allelicFraction > 17/18);       colorList = unphased_color_9of9;
							elseif (allelicFraction > 15/18);   colorList = unphased_color_8of9;
							elseif (allelicFraction > 13/18);   colorList = unphased_color_7of9;
							elseif (allelicFraction > 11/18);   colorList = unphased_color_6of9;
							else                                colorList = unphased_color_5of9;
							end;
						end;
					end;
					chr_SNPdata_colorsC{chr,1}(pos) = chr_SNPdata_colorsC{chr,1}(pos) + colorList(1);
					chr_SNPdata_colorsC{chr,2}(pos) = chr_SNPdata_colorsC{chr,2}(pos) + colorList(2);
					chr_SNPdata_colorsC{chr,3}(pos) = chr_SNPdata_colorsC{chr,3}(pos) + colorList(3);
					chr_SNPdata_countC{ chr  }(pos) = chr_SNPdata_countC{ chr  }(pos) + 1;
				end;
			end;
		end;
	end;

	for chr = 1:num_chrs
		if (chr_in_use(chr) == 1)
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
    
	for chr = 1:num_chrs
		if (chr_in_use(chr) == 1)
			chr_length = chr_size(chr);
			%
			% Determining SNP ratio cutoffs for each chromosome segment.
			%
			for segment = 1:(length(chrCopyNum{chr}))
				histAll_a = [];
				histAll_b = [];
				histAll2  = [];
				%% Construct and smooth a histogram of alleleic fraction data in the segment of interest.
				% phased data is stored into arrays 'histAll_a' and 'histAll_b', since proper phasing is known.
				% unphased data is stored inverted into the second array, since proper phasing is not known.
				for chr_bin = 1:length(CNVplot2{chr})
					%   1 : phased SNP ratio data.
					%   2 : unphased SNP ratio data.
					%   3 : phased SNP position data.
					%   4 : unphased SNP position data.
					ratioData_phased        = chr_SNPdata{chr,1}{chr_bin};
					ratioData_unphased      = chr_SNPdata{chr,2}{chr_bin};
					coordinateData_phased   = chr_SNPdata{chr,3}{chr_bin};
					coordinateData_unphased = chr_SNPdata{chr,4}{chr_bin};
					if (length(ratioData_phased) > 0)
						for SNP_in_bin = 1:length(ratioData_phased)
							if ((coordinateData_phased(SNP_in_bin) > chr_breaks{chr}(segment)*chr_length) && (coordinateData_phased(SNP_in_bin) <= chr_breaks{chr}(segment+1)*chr_length))
								% Ratio data is phased, so it is added twice in its proper orientation (to match density of unphased data below).
								histAll_a = [histAll_a ratioData_phased(SNP_in_bin)];
								histAll_b = [histAll_b ratioData_phased(SNP_in_bin)];
							end;
						end;
					end;
					if (length(ratioData_unphased) > 0)
						for SNP_in_bin = 1:length(ratioData_unphased)
							if ((coordinateData_unphased(SNP_in_bin) > chr_breaks{chr}(segment)*chr_length) && (coordinateData_unphased(SNP_in_bin) <= chr_breaks{chr}(segment+1)*chr_length))
								% Ratio data is unphased, so it is added evenly in both orientations.
								histAll_a = [histAll_a ratioData_unphased(SNP_in_bin)  ];
								histAll_b = [histAll_b 1-ratioData_unphased(SNP_in_bin)];
							end;
						end;
					end;
				end;
				% make a histogram of SNP allelic fractions in segment, then smooth for display.
				histAll                    = [histAll_a histAll_b];
				histAll(histAll == -1)     = [];
				histAll(length(histAll)+1) = 0;
				histAll(length(histAll)+1) = 1;
				% Invert histogram values;
				histAll                    = 1-histAll;
				% generate the histogram.
				data_hist                  = hist(histAll,200);
				% log-scale the histogram.
				data_hist                  = log(data_hist+1);
				data_hist                  = log(data_hist+1);
				% smooth the histogram.
				smoothed                   = smooth_gaussian(data_hist,10,30);
				% make a smoothed version of just the endpoints used to ensure histogram bounds.
				histAll2(1)                = 0;
				histAll2(2)                = 1;
				smoothed2                  = smooth_gaussian(hist(histAll2,200),3,20);
				% subtract the smoothed endpoints from the histogram to remove the influence of the added endpoints.
				smoothed                   = (smoothed-smoothed2);
				if (max(smoothed) > 0)
					smoothed           = smoothed/max(smoothed);
				end;
				smoothed__{chr}{segment}   = smoothed;
                
				% Define cutoffs between Gaussian fits.
				segment_copyNum = round(chrCopyNum{chr}(segment));
				saveName = ['allelic_ratios.chr_' num2str(chr) '.segment_' num2str(segment) ];
				[peaks,actual_cutoffs__{chr}{segment},mostLikelyGaussians__{chr}{segment}] = FindGaussianCutoffs_3(workingDir,saveName, chr,segment, segment_copyNum,smoothed__{chr}{segment}, false);

				ratio_cutoffs = (actual_cutoffs__{chr}{segment}-1)/199;
				fprintf(['chr                   = ' num2str(chr) '; segment = ' num2str(segment) '.\n']);
				fprintf(['segment copy num      = ' num2str(segment_copyNum) '.\n']);
				fprintf(['allelic ratio cutoffs = [' num2str(ratio_cutoffs)                       ']\n']);
				fprintf(['most likely Gaussians = [' num2str(mostLikelyGaussians__{chr}{segment}) ']\n\n']);
			end;
		end;
	end;
    
	if (true)
		%%
		%% Reset tracking vectors.
		%%
		chr_SNPdata_colorsC = [];
		chr_SNPdata_colorsP = [];
		chr_SNPdata_countC  = [];
		chr_SNPdata_countP  = [];
		for chr = 1:num_chrs
			if (chr_in_use(chr) == 1)
				chr_length = ceil(chr_size(chr)/bases_per_bin);
				for j = 1:3
					chr_SNPdata_colorsC{chr,j} = zeros(chr_length,1);
					chr_SNPdata_colorsP{chr,j} = zeros(chr_length,1);
				end;
				chr_SNPdata_countC{chr}        = zeros(chr_length,1);
				chr_SNPdata_countP{chr}        = zeros(chr_length,1);
			end;
		end;
		%%
		%% Define new colors for SNPs.
		%%
		for chr = 1:num_chrs
			if (chr_in_use(chr) == 1)
				if (length(C_chr_count{chr}) > 1)
					%
					% Determining colors for each SNP coordinate from calculated cutoffs.
					%
					saveName                                = ['allelic_ratios.chr_' num2str(chr) '.segment_' num2str(segment) ];
					for SNP = 1:length(C_chr_SNP_data_positions{chr})  % 'length(C_chr_SNP_data_positions)' is the number of SNPs per chromosome.
						coordinate                      = C_chr_SNP_data_positions{chr}(SNP);
						pos                             = ceil(coordinate/new_bases_per_bin);
						localCopyEstimate               = round(CNVplot2{chr}(pos)*ploidy*ploidyAdjust);
						baseCall                        = C_chr_baseCall{        chr}{SNP};
						homologA                        = C_chr_SNP_homologA{    chr}{SNP};
						homologB                        = C_chr_SNP_homologB{    chr}{SNP};
						flipper                         = C_chr_SNP_flipHomologs{chr}(SNP);
						allelic_ratio                   = C_chr_SNP_data_ratios{ chr}(SNP);
						if (flipper == 10)                         % Variable 'flipper' value of '10' indicates no phasing information is available in the hapmap.
							baseCall                = 'Z';     % Variable 'baseCall' value of 'Z' will prevent either hapmap allele from matching and so unphased ratio colors will be used in the following section.
							chr_SNPdata{chr,2}{pos} = [chr_SNPdata{chr,2}{pos} allelic_ratio 1-allelic_ratio];
							chr_SNPdata{chr,4}{pos} = [chr_SNPdata{chr,4}{pos} coordinate    coordinate     ];
						elseif (flipper == 1)
							temp                    = homologA;
							homologA                = homologB;
							homologB                = temp;
							chr_SNPdata{chr,1}{pos} = [chr_SNPdata{chr,1}{pos} 1-allelic_ratio];
							chr_SNPdata{chr,3}{pos} = [chr_SNPdata{chr,3}{pos} coordinate];
							allelic_ratio           = 1-allelic_ratio;
						else % (flipper == 0)
							chr_SNPdata{chr,1}{pos} = [chr_SNPdata{chr,1}{pos} allelic_ratio];
							chr_SNPdata{chr,3}{pos} = [chr_SNPdata{chr,3}{pos} coordinate];
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
						segment_copyNum                 = round(chrCopyNum{chr}(segmentID));
						actual_cutoffs                  = actual_cutoffs__{chr}{segmentID};
						mostLikelyGaussians             = mostLikelyGaussians__{chr}{segmentID};

						% Convert the allelic ratio to the range[1..200].
						if (flipper == 1)
							SNPratio_int            = (allelic_ratio)*199+1;
						else
							SNPratio_int            = (1-allelic_ratio)*199+1;
						end;

						% Identify the allelic ratio region containing the SNP.
						cutoffs                         = [1 actual_cutoffs 200];
						ratioRegionID                   = 0;
						for ratioRegion = 1:(length(cutoffs)-1)
							cutoff_start            = cutoffs(ratioRegion  );
							cutoff_end              = cutoffs(ratioRegion+1);
							if (SNPratio_int > cutoff_start) && (SNPratio_int <= cutoff_end)
								ratioRegionID   = ratioRegion;
							end;
						end;

						if (segment_copyNum <= 0);                          colorList = colorNoData;
						elseif (segment_copyNum == 1)
							% allelic fraction cutoffs: [0.50000] => [A B]
							if (baseCall == homologA);
								if (ratioRegionID == 2);            colorList = colorA;
								else                                colorList = colorB;
								end;
							elseif (baseCall == homologB);
								if (ratioRegionID == 2);            colorList = colorB;
								else                                colorList = colorA;
								end;
							else
								if (ratioRegionID == 2);            colorList = unphased_color_1of1;
								else                                colorList = unphased_color_1of1;
								end;
							end;
						elseif (segment_copyNum == 2)
							%   allelic fraction cutoffs: [0.25000 0.75000] => [AA AB BB]
							%   cutoffs_1 from Gaussian fittings = [0.18342 0.19347 0.18844 0.16332 0.16332 0.15327 0.12814] => 0.1676;
							%   cutoffs_2 from Gaussian fittings = [0.82663 0.82161 0.80151 0.84171 0.83668 0.84673 0.84673] => 0.8317;
							% fprintf(['allelic_ratio = ' num2str(allelic_ratio) '\tchr = ' num2str(chr) '\tbp = ' num2str(coordinate) '\n']);
							if (baseCall == homologA)
								if (ratioRegionID == 3);            colorList = colorAA;
								elseif (ratioRegionID == 2);        colorList = colorAB;
								else                                colorList = colorBB;
								end;
							elseif (baseCall == homologB)
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
							if (baseCall == homologA)
								if (ratioRegionID == 4);            colorList = colorAAA;
								elseif (ratioRegionID == 3);        colorList = colorAAB;
								elseif (ratioRegionID == 2);        colorList = colorABB;
								else                                colorList = colorBBB;
								end;
							elseif (baseCall == homologB)
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
							if (baseCall == homologA)
								if (ratioRegionID == 5);            colorList = colorAAAA;
								elseif (ratioRegionID == 4);        colorList = colorAAAB;
								elseif (ratioRegionID == 3);        colorList = colorAABB;
								elseif (ratioRegionID == 2);        colorList = colorABBB;
								else                                colorList = colorBBBB;
								end;
							elseif (baseCall == homologB)
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
							if (baseCall == homologA)
								if (ratioRegionID == 6);            colorList = colorAAAAA;
								elseif (ratioRegionID == 5);        colorList = colorAAAAB;
								elseif (ratioRegionID == 4);        colorList = colorAAABB;
								elseif (ratioRegionID == 3);        colorList = colorAABBB;
								elseif (ratioRegionID == 2);        colorList = colorABBBB;
								else                                colorList = colorBBBBB;
								end;
							elseif (baseCall == homologB)
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
							if (baseCall == homologA)
								if (ratioRegionID == 7);            colorList = colorAAAAAA;
								elseif (ratioRegionID == 6);        colorList = colorAAAAAB;
								elseif (ratioRegionID == 5);        colorList = colorAAAABB;
								elseif (ratioRegionID == 4);        colorList = colorAAABBB;
								elseif (ratioRegionID == 3);        colorList = colorAABBBB;
								elseif (ratioRegionID == 2);        colorList = colorABBBBB;
								else                                colorList = colorBBBBBB;
								end;
							elseif (baseCall == homologB)
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
							if (baseCall == homologA)
								if (ratioRegionID == 8);            colorList = colorAAAAAAA;
								elseif (ratioRegionID == 7);        colorList = colorAAAAAAB;
								elseif (ratioRegionID == 6);        colorList = colorAAAAABB;
								elseif (ratioRegionID == 5);        colorList = colorAAAABBB;
								elseif (ratioRegionID == 4);        colorList = colorAAABBBB;
								elseif (ratioRegionID == 3);        colorList = colorAABBBBB;
								elseif (ratioRegionID == 2);        colorList = colorABBBBBB;
								else                                colorList = colorBBBBBBB;
								end;
							elseif (baseCall == homologB)
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
							if (baseCall == homologA)
								if (ratioRegionID == 9);            colorList = colorAAAAAAAA;
								elseif (ratioRegionID == 8);        colorList = colorAAAAAAAB;
								elseif (ratioRegionID == 7);        colorList = colorAAAAAABB;
								elseif (ratioRegionID == 6);        colorList = colorAAAAABBB;
								elseif (ratioRegionID == 5);        colorList = colorAAAABBBB;
								elseif (ratioRegionID == 4);        colorList = colorAAABBBBB;
								elseif (ratioRegionID == 3);        colorList = colorAABBBBBB;
								elseif (ratioRegionID == 2);        colorList = colorABBBBBBB;
								else                                colorList = colorBBBBBBBB;
								end;
							elseif (baseCall == homologB)
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
							% allelic fraction cutoffs: [0.05556 0.16667 0.27778 0.38889 0.50000 0.61111 0.72222 0.83333 0.94444] => [AAAAAAAAA AAAAAAAAB AAAAAAABB AAAAAABBB AAAAABBBB AAAABBBBB AAABBBBBB AABBBBBBB ABBBBBBBB BBBBBBBBB]
							if (baseCall == homologA)
								if (ratioRegionID == 10);           colorList = colorAAAAAAAAA;
								elseif (ratioRegionID == 9);        colorList = colorAAAAAAAAB;
								elseif (ratioRegionID == 8);        colorList = colorAAAAAAABB;
								elseif (ratioRegionID == 7);        colorList = colorAAAAAABBB;
								elseif (ratioRegionID == 6);        colorList = colorAAAAABBBB;
								elseif (ratioRegionID == 5);        colorList = colorAAAABBBBB;
								elseif (ratioRegionID == 4);        colorList = colorAAABBBBBB;
								elseif (ratioRegionID == 3);        colorList = colorAABBBBBBB;
								elseif (ratioRegionID == 2);        colorList = colorABBBBBBBB;
								else                                colorList = colorBBBBBBBBB;
								end;
							elseif (baseCall == homologB)
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
						chr_SNPdata_colorsC{chr,1}(pos) = chr_SNPdata_colorsC{chr,1}(pos) + colorList(1);
						chr_SNPdata_colorsC{chr,2}(pos) = chr_SNPdata_colorsC{chr,2}(pos) + colorList(2);
						chr_SNPdata_colorsC{chr,3}(pos) = chr_SNPdata_colorsC{chr,3}(pos) + colorList(3);
						chr_SNPdata_countC{ chr  }(pos) = chr_SNPdata_countC{ chr  }(pos) + 1;
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
	 end;
else
%
% Run when compared vs. a parent dataset or vs. itself.
%
%%%% chr_SNPdata{chr,2}(SNP) = child data.
%%%% chr_SNPdata{chr,4}(SNP) = parent data.
	for chr = 1:num_chrs
		if (chr_in_use(chr) == 1)
			if (length(C_chr_count{chr}) > 1)
				for SNP = 1:length(C_chr_count{chr})
					pos = ceil(C_chr_SNP_data_positions{chr}(SNP)/new_bases_per_bin);
					if (C_chr_SNP_data_ratios{chr}(SNP) < chr_SNPdata{chr,2}(pos))
						chr_SNPdata{chr,2}(pos)      = C_chr_SNP_data_ratios{chr}(SNP);
						colorList                    = [1.0 1.0 1.0];
						chr_SNPdata_colorsC{chr,1}(SNP) = colorList(1);
						chr_SNPdata_colorsC{chr,2}(SNP) = colorList(2);
						chr_SNPdata_colorsC{chr,3}(SNP) = colorList(3);
					end;
				end;
			end;
			if (length(P_chr_count{chr}) > 1)
				for i = 1:length(P_chr_count{chr})
					pos = ceil(P_chr_SNP_data_positions{chr}(SNP)/new_bases_per_bin);
					if (P_chr_SNP_data_ratios{chr}(SNP) < chr_SNPdata{chr,4}(pos))
						chr_SNPdata{chr,4}(pos)      = P_chr_SNP_data_ratios{chr}(SNP);
						colorList                    = [1.0 1.0 1.0];
						chr_SNPdata_colorsP{chr,1}(SNP) = colorList(1);
						chr_SNPdata_colorsP{chr,2}(SNP) = colorList(2);
						chr_SNPdata_colorsP{chr,3}(SNP) = colorList(3);
					end;
				end;
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
			data_phased = [data_phased chr_SNPdata{chr,1}{i}];
		end;
		for i = 1:length(chr_SNPdata{chr,2})
			data_unphased = [data_unphased chr_SNPdata{chr,2}{i}];
		end;
		histogram_phased         = hist([data_phased   0 1],200);
		histogram_unphased       = hist([data_unphased 0 1],200);

		subplot(2,num_chrs,chr);
		hold on;
		plot(1:200, log(histogram_unphased+1),  'Color',[1.0 0.0 0.0]);
		plot(1:200, log(histogram_phased+1),'Color',[1/3 1/3 1/3]);
		ylim([0 6]);
		title([chr_label{chr}]);
		set(gca,'XTick',[0 50 100 150 200]);
		set(gca,'XTickLabel',{'0','1/4','1/2','3/4','1'});
		ylabel('log(data count)');
		xlabel('allelic ratio');
		hold off;

		subplot(2,num_chrs,chr+num_chrs);
		hold on;
		plot(1:200, histogram_unphased,  'Color',[1.0 0.0 0.0]);
		plot(1:200, histogram_phased,'Color',[1/3 1/3 1/3]);
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
%saveas(histogram_fig, [projectDir 'fig.allelic_fraction_histogram.eps'], 'epsc');
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
			for chr_bin = 1:ceil(chr_size(chr)/new_bases_per_bin)
				colorR   = chr_SNPdata_colorsC{chr,1}(chr_bin);
				colorG   = chr_SNPdata_colorsC{chr,2}(chr_bin);
				colorB   = chr_SNPdata_colorsC{chr,3}(chr_bin);
				if (colorR <= 1) || (colorG <= 1) || (colorB <= 1)
					plot([chr_bin chr_bin], [0 maxY],'Color',[colorR colorG colorB]);
				end;
			end;
		% Following variation will draw experimetnal vs. reference dataset like the above vs. hapmap figure.
		%	elseif (useParent)   % an experimental dataset vs. a reference dataset.
		%		for chr_bin = 1:ceil(chr_size(chr)/new_bases_per_bin)
		%			colorR   = chr_SNPdata_colorsP{chr,1}(chr_bin);
		%			colorG   = chr_SNPdata_colorsP{chr,2}(chr_bin);
		%			colorB   = chr_SNPdata_colorsP{chr,3}(chr_bin);
		%			if (colorR <= 1) || (colorG <= 1) || (colorB <= 1)
		%				plot([chr_bin chr_bin], [0 maxY],'Color',[colorR colorG colorB]);
		%			end;
		%		end;
		elseif (useParent)   % an experimental dataset vs. a reference dataset.
			for chr_bin = 1:ceil(chr_size(chr)/new_bases_per_bin)
				datumY_C = chr_SNPdata{chr,2}(chr_bin)*maxY;
				datumY_P = chr_SNPdata{chr,4}(chr_bin)*maxY;
				plot([chr_bin/2 chr_bin/2], [maxY datumY_C     ],'Color',[1.0 0.0 0.0]);
				plot([chr_bin/2 chr_bin/2], [0    maxY-datumY_P],'Color',[1/3 1/3 1/3]);
			end;
		else   % only a reference dataset.
			for chr_bin = 1:ceil(chr_size(chr)/new_bases_per_bin)
				datumY_P = chr_SNPdata{chr,4}(chr_bin)*maxY;
				plot([chr_bin/2 chr_bin/2], [maxY datumY_P     ],'Color',[1/3 1/3 1/3]);
				plot([chr_bin/2 chr_bin/2], [0    maxY-datumY_P],'Color',[1/3 1/3 1/3]);
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
				for chr_bin = 1:ceil(chr_size(chr)/new_bases_per_bin)
					colorR   = chr_SNPdata_colorsC{chr,1}(chr_bin);
					colorG   = chr_SNPdata_colorsC{chr,2}(chr_bin);
					colorB   = chr_SNPdata_colorsC{chr,3}(chr_bin);
					if (colorR < 1) || (colorG < 1) || (colorB < 1)
						plot([chr_bin chr_bin], [0 maxY],'Color',[colorR colorG colorB]);
					end;
				end;
			% Following variation will draw experimetnal vs. reference dataset like the above vs. hapmap figure.
			%	elseif (useParent)   % an experimental dataset vs. a reference dataset.
			%		for i = 1:ceil(chr_size(chr)/new_bases_per_bin)
			%			colorR   = chr_SNPdata_colorsC{chr,1}(chr_bin);
			%			colorG   = chr_SNPdata_colorsC{chr,2}(chr_bin);
			%			colorB   = chr_SNPdata_colorsC{chr,3}(chr_bin);
			%			if (colorR < 1) || (colorG < 1) || (colorB < 1)
			%				plot([i i], [0 maxY],'Color',[colorR colorG colorB]);
			%			end;
			%		end;
			elseif (useParent)   % an experimental dataset vs. a reference dataset.
				for chr_bin = 1:ceil(chr_size(chr)/new_bases_per_bin)
					datumY_C = chr_SNPdata{chr,2}(chr_bin)*maxY;
					datumY_P = chr_SNPdata{chr,4}(chr_bin)*maxY;
					plot([chr_bin/2 chr_bin/2], [maxY datumY_C     ],'Color',[1.0 0.0 0.0]);
					plot([chr_bin/2 chr_bin/2], [0    maxY-datumY_P],'Color',[1/3 1/3 1/3]);
				end;
			else   % only a reference dataset.
				for chr_bin = 1:ceil(chr_size(chr)/new_bases_per_bin)
					datumY_C = chr_SNPdata{chr,2}(chr_bin)*maxY;
					datumY_P = chr_SNPdata{chr,4}(chr_bin)*maxY;
					plot([chr_bin/2 chr_bin/2], [maxY datumY_P     ],'Color',[1/3 1/3 1/3]);
					plot([chr_bin/2 chr_bin/2], [0    maxY-datumY_P],'Color',[1/3 1/3 1/3]);
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
%saveas(fig,        [projectDir 'fig.allelic_ratio-map.c1.eps'], 'epsc');
saveas(fig,        [projectDir 'fig.allelic_ratio-map.c1.png'], 'png');
delete(fig);

set(Linear_fig,'PaperPosition',[0 0 8 0.62222222]*2);
%saveas(Linear_fig, [projectDir 'fig.allelic_ratio-map.c2.eps'], 'epsc');
saveas(Linear_fig, [projectDir 'fig.allelic_ratio-map.c2.png'], 'png');
delete(Linear_fig);

%%================================================================================================
% end stuff
%=================================================================================================
end

