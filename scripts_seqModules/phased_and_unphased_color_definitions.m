colorNoData = [1.0     1.0     1.0    ]; % used when no data is available for the bin.
colorInit   = [0.5     0.5     0.5    ]; % external; used in blending at ends of chr.
colorHET    = [0.66667 0.66667 0.66667]; % near 1:1 ratio SNPs
colorHOM    = [1.0     0.0     0.0    ]; % Hom SNPs;

if (useHapmap == true)
	%% Load color names defined for hapmap;
	colorsFile = [hapmapDir 'colors.txt'];
	if (exist(colorsFile,'file') == 2)
		colors_fid = fopen([main_dir 'users/' hapmapUser '/hapmaps/' hapmap '/colors.txt'], 'r');
		% The swapped colors are to correct for a polarity mistake in the python preprocessing steps.
		%    correcting the error there would require reprocessing all current datasets.
		colorB_string = fgetl(colors_fid);
		colorA_string = fgetl(colors_fid);
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
else
	% Haplotype map is not in use.
	if (strcmp(project,hapmap) == 1)
		% The 'project' is the same as the 'hapmap'/'parent'.
		homolog_a_color = [0.66667 0.66667 0.66667];
		homolog_b_color = [0.66667 0.66667 0.66667];
	else
		% The 'project' is different than the 'hapmap'/'parent'.
		homolog_a_color = [1.0 0.0 0.0];
		homolog_b_color = [1.0 0.0 0.0];
	end;
end;
hom_color    = [1.0     0.0     0.0    ]; % completely homozygous.
het_color    = [0.66667 0.66667 0.66667]; % heterozygous.
oddHet_color = [0.0     1.0     0.0    ]; % non-heterozygous data that isn't 100 hom.

%%%%%%%% phased data colors.

% haploid colors.
colorA  = homolog_a_color;
colorB  = homolog_b_color;
% diploid colors.
colorAA = homolog_a_color;
colorAB = het_color;
colorBB = homolog_b_color;
% triploid colors.
colorAAA= homolog_a_color;
colorAAB= homolog_a_color*2/3 + homolog_b_color*1/3;
colorABB= homolog_a_color*1/3 + homolog_b_color*2/3;
colorBBB= homolog_b_color;
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

%%%%%%%% unphased colors.

% haploid colors.
unphased_color_1of1 = hom_color;
% diploid colors.
unphased_color_2of2 = hom_color;
unphased_color_1of2 = het_color;
% triploid colors.
unphased_color_3of3 = hom_color;
unphased_color_2of3 = hom_color*2/3 + het_color*1/3;
% tetraploid colors.
unphased_color_4of4 = hom_color;
unphased_color_3of4 = hom_color*3/4 + het_color*1/4;
unphased_color_2of4 = het_color;
% pentaploid colors.
unphased_color_5of5 = hom_color;
unphased_color_4of5 = hom_color*4/5 + het_color*1/5;
unphased_color_3of5 = hom_color*3/5 + het_color*2/5;
% hexaploid colors.
unphased_color_6of6 = hom_color;
unphased_color_5of6 = hom_color*5/6 + het_color*1/6;
unphased_color_4of6 = hom_color*4/6 + het_color*2/6;
unphased_color_3of6 = het_color;
% heptaploid colors.
unphased_color_7of7 = hom_color;
unphased_color_6of7 = hom_color*6/7 + het_color*1/7;
unphased_color_5of7 = hom_color*5/7 + het_color*2/7;
unphased_color_4of7 = hom_color*4/7 + het_color*3/7;
% octaploid colors.
unphased_color_8of8 = hom_color;
unphased_color_7of8 = hom_color*7/8 + het_color*1/8;
unphased_color_6of8 = hom_color*6/8 + het_color*2/8;
unphased_color_5of8 = hom_color*5/8 + het_color*3/8;
unphased_color_4of8 = het_color;
% nonaploid colors.
unphased_color_9of9 = hom_color;
unphased_color_8of9 = hom_color*8/9 + het_color*1/9;
unphased_color_7of9 = hom_color*7/9 + het_color*2/9;
unphased_color_6of9 = hom_color*6/9 + het_color*3/9;
unphased_color_5of9 = hom_color*5/9 + het_color*4/9;
