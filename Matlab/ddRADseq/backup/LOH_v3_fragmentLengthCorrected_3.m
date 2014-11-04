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

Centromere_format           = 0;
Chr_max_width               = 0.8;
colorBars                   = true;
blendColorBars              = false;
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

% bases_per_bin			= 5000;
bases_per_bin			= max(chr_size)/700;
chr_length_scale_multiplier	= 1/bases_per_bin;
CGD_bases_per_bin		= 1000;
CGD_chr_length_scale_multiplier	= 1/CGD_bases_per_bin;

%%=========================================================================
%%= No further control variables below. ===================================
%%=========================================================================
if (strcmp(parent_name,child_name) == 1)
    fprintf(['\nGenerating SNP-map figure from ' parent_name ' genome data.\n']);
else
    fprintf(['\nGenerating LOH-map figure from ' parent_name ' and ' child_name ' genome data.\n']);
end;

%% Sanitize user input of euploid state.
% ploidyBase = round(str2num(ploidyBaseString));
% if (ploidyBase > 4);   ploidyBase = 4;   end;
% if (ploidyBase < 1);   ploidyBase = 1;   end;
ploidyBase = 4;
fprintf(['\nEuploid base = "' num2str(ploidyBase) '"\n']);
    
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

%% Load pre-processed ddRADseq fragment repetitiveness data for genome.
if (exist([workingDir 'main_script_dir/repetitiveness_files/' genome  '.fragment_repet_data.mat'],'file') == 0)
    fprintf('## Process genome for repetitiveness first.');
    error('## Process genome for repetitiveness first.');
else
    load([workingDir 'main_script_dir/repetitiveness_files/' genome  '.fragment_repet_data.mat']);
end;


if (exist([workingDir 'matlab_dir/' child_name '.Lowess1.mat'],'file') == 0)
    fprintf('## Process project for CNV first.');
    error('## Process project for CNV first.');
else
    fprintf('\tLOWESS fitting (project) : ave read count vs. fragment length.\n');
    load([workingDir 'matlab_dir/' child_name '.Lowess1.mat']);
end;
newX1_proj = newX1;
newY1_proj = newY1;

if (strcmp(parent_name,child_name) == 0)
    if (exist([workingDir 'matlab_dir/' parent_name '.Lowess1.mat'],'file') == 0)
	fprintf('## Process reference for CNV first.');
	error('## Process reference for CNV first.');
    else
	fprintf('\tLOWESS fitting (reference) : ave read count vs. fragment length.\n');
	load([workingDir 'matlab_dir/' parent_name '.Lowess1.mat']);
    end;
    newX1_ref = newX1;
    newY1_ref = newY1;
else
    newX1_ref = newX1_proj;
    newY1_ref = newY1_proj;
end;


%% Standardize 'fragment_data' data structure to hold results.
fragment_data = project_fragments_CNV;
numFragments  = length(fragment_data);
for frag = 1:numFragments
    fragment_data(frag).ref_read_count = reference_fragments_CNV(frag).read_count;  
    fragment_data(frag).ref_read_max   = reference_fragments_CNV(frag).read_max;
    fragment_data(frag).ref_read_ave   = reference_fragments_CNV(frag).read_ave;
    fragment_data(frag).repet_count    = fragments_repet(frag).repet_count;
    fragment_data(frag).repet_max      = fragments_repet(frag).repet_max;
    fragment_data(frag).repet_ave      = fragments_repet(frag).repet_ave;
end;


%%================================================================================================
% Load SNP/LOH data.
%-------------------------------------------------------------------------------------------------

%% Load pre-processed ddRADseq fragment SNP data for project.
if (exist([workingDir 'matlab_dir/' child_name '.fragment_SNP_data.mat'],'file') == 0)
    % ### chr_num     bp_start        be_end  copyNum binData
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
        test = sscanf(tline, '%s',1);
        if (strcmp(test,'###') == 0) && (strcmp(test,'##') == 0)	% If the first space-delimited string is not '###' and is not '##', then process it as a valid data line.
            % The number of valid lines found so far...  the number of restriction fragments with data so far.
            count = count + 1;

            % Chr ID of bp.
            chr_num = sscanf(tline, '%s',1);  
    
            % bp coordinate of fragment start.
            new_string = sscanf(tline, '%s',  2 );
            for i = 1:size(sscanf(tline,'%s', 1 ),2);   new_string(1) = [];   end;
            bp_start = new_string;
        
            % bp coordinate of fragment end.
            new_string = sscanf(tline, '%s',  3 );
            for i = 1:size(sscanf(tline,'%s', 2 ),2);   new_string(1) = [];   end;
            bp_end = new_string;

            % corrected copy number estimate for fragment.
	    new_string = sscanf(tline, '%s',  4 );
            for i = 1:size(sscanf(tline,'%s', 3 ),2);   new_string(1) = [];   end;
            copyNum = new_string;
    
            % ratio data bins, separated by ':'.   if first character is '*', then no ratio bins were calculated.
            new_string = sscanf(tline, '%s',  5 );
            for i = 1:size(sscanf(tline,'%s', 4 ),2);   new_string(1) = [];   end;
            data_bins = new_string;  
        
            % Add fragment data to data structure.
            fragments_SNP(count).chr     = str2num(chr_num);
            fragments_SNP(count).startbp = str2num(bp_start);
            fragments_SNP(count).endbp   = str2num(bp_end);
            fragments_SNP(count).length  = str2num(bp_end) - str2num(bp_start) + 1;
            fragments_SNP(count).copyNum = str2num(copyNum);
            fragments_SNP(count).data    = data_bins;
        end;
    end;
    fclose(data_RADseq);
    save([workingDir 'matlab_dir/' child_name '.fragment_SNP_data.mat'], 'fragments_SNP');
else
    load([workingDir 'matlab_dir/' child_name '.fragment_SNP_data.mat']);
end;
child_fragments_SNP = fragments_SNP;
clear fragments_SNP;

if (strcmp(parent_name,child_name) == 0)
    %% Load pre-processed ddRADseq fragment SNP data for reference.
    if (exist([workingDir 'matlab_dir/' parent_name '.fragment_SNP_data.mat'],'file') == 0)
	% ### chr_num     bp_start        be_end  copyNum binData
	% 1       14001   14291   2       20:0:0
	% 1       15900   16123   3       18:0:0
	% 1       21632   21834   2       14:0:0
	fprintf('Loading reference results from Python script, which pre-processed the putative_SNPs relative to genome restriction fragments and corrected fragment copy number estimates.\n');
	datafile_RADseq  = [workingDir 'pileup_dir/' parent_name '_RADseq_digest_analysis_SNP.txt'];
	data_RADseq      = fopen(datafile_RADseq);
	count            = 0;
	fragments_SNP    = [];
	while ~feof(data_RADseq)
	    % Load fragment data from pre-processed text file, single line.
	    tline = fgetl(data_RADseq);
            
	    % check if line is a comment.
	    test = sscanf(tline, '%s',1);
	    if (strcmp(test,'###') == 0) && (strcmp(test,'##') == 0)        % If the first space-delimited string is not '###' and is not '##', then process it as a valid data line.
		% The number of valid lines found so far...  the number of restriction fragments with data so far.
		count = count + 1;

		% Chr ID of bp.
		chr_num = sscanf(tline, '%s',1);

		% bp coordinate of fragment start.
		new_string = sscanf(tline, '%s',  2 );
		for i = 1:size(sscanf(tline,'%s', 1 ),2);   new_string(1) = [];   end;
		bp_start = new_string;

		% bp coordinate of fragment end.
		new_string = sscanf(tline, '%s',  3 );
		for i = 1:size(sscanf(tline,'%s', 2 ),2);   new_string(1) = [];   end;
		bp_end = new_string;

		% corrected copy number estimate for fragment.
		new_string = sscanf(tline, '%s',  4 );
		for i = 1:size(sscanf(tline,'%s', 3 ),2);   new_string(1) = [];   end;
		copyNum = new_string;

		% ratio data bins, separated by ':'.   if first character is '*', then no ratio bins were calculated.
		new_string = sscanf(tline, '%s',  5 );
		for i = 1:size(sscanf(tline,'%s', 4 ),2);   new_string(1) = [];   end;
		data_bins = new_string;

		% Add fragment data to data structure.
		fragments_SNP(count).chr     = str2num(chr_num);
		fragments_SNP(count).startbp = str2num(bp_start);
		fragments_SNP(count).endbp   = str2num(bp_end);
		fragments_SNP(count).length  = str2num(bp_end) - str2num(bp_start) + 1;
		fragments_SNP(count).copyNum = str2num(copyNum);
		fragments_SNP(count).data    = data_bins;		% data bins separated by ":".   first character is "*" for unusable data.
	    end;
	end;
	fclose(data_RADseq);
	save([workingDir 'matlab_dir/' parent_name '.fragment_SNP_data.mat'], 'fragments_SNP');
    else
	load([workingDir 'matlab_dir/' parent_name '.fragment_SNP_data.mat']);
    end;
    parent_fragments_SNP = fragments_SNP;
    clear fragments_SNP;
else
    parent_fragments_SNP = child_fragments_SNP;
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

fprintf(['all fragments    = ' num2str(numFragments) ':' num2str(length(fragment_data)) '\n']);
numSNPfragments = 0;
for frag = 1:numFragments;
    if (fragment_data(frag).usable == 1)
	numSNPfragments = numSNPfragments + 1;
    end;
end;
fprintf(['usable fragments = ' num2str(numSNPfragments) ':' num2str(length(child_fragments_SNP)) '\n']);


% Initializes vectors used to hold SNP number data.
for chr = 1:num_chr   % number of chrs.
    for j = 1:2       % 2 categories tracked : total read counts in bin and number of data entries in region (size of bin in base-pairs).
	chr_SNPdata_RADseq{chr,j} = zeros(1,ceil(chr_size(chr)/bases_per_bin));
    end;
end;
SNPfrag = 0;
%% Gather SNP data into bins for map display.
if (exist([workingDir 'matlab_dir/' child_name '.corrected_SNP_' SNP_verString '.mat'],'file') == 0)
    fprintf('\nMAT file containing project SNP information not found, regenerating from prior data files.\n');

    % fragments_SNP(count).data    = data_bins;           % data bins separated by ":".   first character is "*" for unusable data.

    % Convert bias corrected SNP ratio data per restriction digest fragment into bin-aligned data for plotting
    fprintf('\n# Adding fragment corrected SNP ratio values to map bins.');

    fprintf(['\nnumFragments = ' num2str(numFragments) '\n']);
    for frag = 1:numFragments
        if (fragment_data(frag).usable == 1)
	    SNPfrag    = SNPfrag+1;
            % Load important data from fragments data structure.
            chrID      = fragment_data(frag).chr;
            posStart   = fragment_data(frag).startbp;
            posEnd     = fragment_data(frag).endbp;
	    copyNum    = child_fragments_SNP(SNPfrag).copyNum;
	    data       = child_fragments_SNP(SNPfrag).data;
            fragLength = posEnd - posStart + 1;

            % Identify locations of fragment end in relation to plotting bins.
            val1 = ceil(posStart/bases_per_bin);
            val2 = ceil(posEnd  /bases_per_bin);

            if (chrID > 0)
		dataBins = strsplit(data,':');
		numBins  = length(dataBins);

		%% Calculate the SNPvalue for the fragment.
		%	Each putative_SNP contributes to the SNPvalue depending on how close its ratio call is to heterozygous 1:1.
		%	Homozygous calls contribute 0 per base pair.   Heterozygous 1:1 calls contribute 1 per base pair.
		%	Intermediate ratios contribute intermediate values.
		%  The SNPvalue for each fragment is a measure of how heterozygous the fragment is.
		if (strcmp(dataBins{1},'*') == 1)
		    % This bin has too high a copy number to have sensical ratios.
		    SNPvalue = 0;
		else
		    % This bin has a copy number between 1 & 8.
		    switch copyNum
			case 0	% No useful data.
			    SNPvalue = 0;
			case 1	% No useful data.
			    SNPvalue = 0;
			case 2	% aa, ab.
			    % SNPvalue = 0*str2num(dataBins{1}) + 1*str2num(dataBins{2});
			    if (str2num(dataBins{2}) > 0)
				SNPvalue = 1;
			    else
				SNPvalue = 0;
			    end;
			case 3	% aaa, aab, abb.
			    % SNPvalue = 0*str2num(dataBins{1}) + 2/3*str2num(dataBins{2});
			    if (str2num(dataBins{2}) > 0)
                                SNPvalue = 1;
			    else
				SNPvalue = 0;
                            end;
			case 4	% aaaa, aaab, aabb.
			    % SNPvalue = 0*str2num(dataBins{1}) + 1/2*str2num(dataBins{2}) + 1*str2num(dataBins{3});
			    if (str2num(dataBins{3}) > 0)
                                SNPvalue = 1;
			    elseif (str2num(dataBins{2}) > 0)
                                SNPvalue = 1/2;
			    else
				SNPvalue = 0;
                            end;
			case 5	% aaaaa, aaaab, aaabb.
			    % SNPvalue = 0*str2num(dataBins{1}) + 2/5*str2num(dataBins{2}) + 4/5*str2num(dataBins{3});
			    if (str2num(dataBins{3}) > 0)
                                SNPvalue = 4/5;
                            elseif (str2num(dataBins{2}) > 0)
                                SNPvalue = 2/5;
                            else
                                SNPvalue = 0;
                            end;
			case 6	% aaaaaa, aaaaab, aaaabb, aaabbb.
			    % SNPvalue = 0*str2num(dataBins{1}) + 1/3*str2num(dataBins{2}) + 2/3*str2num(dataBins{3}) + 1*str2num(dataBins{4});
			    if (str2num(dataBins{4}) > 0)
                                SNPvalue = 1;
                            elseif (str2num(dataBins{3}) > 0)
                                SNPvalue = 2/3;
			    elseif (str2num(dataBins{2}) > 0)
                                SNPvalue = 1/3; 
                            else
                                SNPvalue = 0;
                            end;
			case 7	% aaaaaaa, aaaaaab, aaaaabb, aaaabbb.
			    % SNPvalue = 0*str2num(dataBins{1}) + 2/7*str2num(dataBins{2}) + 4/7*str2num(dataBins{3}) + 6/7*str2num(dataBins{4});
			    if (str2num(dataBins{4}) > 0)
                                SNPvalue = 6/7;
                            elseif (str2num(dataBins{3}) > 0)
                                SNPvalue = 4/7;
			    elseif (str2num(dataBins{2}) > 0)
                                SNPvalue = 2/7; 
                            else
                                SNPvalue = 0;
                            end;
			case 8	% aaaaaaaa, aaaaaaab, aaaaaabb, aaaaabbb, aaaabbbb.
			    % SNPvalue = 0*str2num(dataBins{1}) + 1/4*str2num(dataBins{2}) + 1/2*str2num(dataBins{3}) + 3/4*str2num(dataBins{4}) + 1*str2num(dataBins{5});
			    if (str2num(dataBins{5}) > 0)
                                SNPvalue = 1;
                            elseif (str2num(dataBins{4}) > 0)
                                SNPvalue = 3/4;
			    elseif (str2num(dataBins{3}) > 0)
                                SNPvalue = 1/2;
			    elseif (str2num(dataBins{2}) > 0)
                                SNPvalue = 1/4; 
                            else
                                SNPvalue = 0;
                            end;
		    end;
		end;

                % 'count' is average read count across fragment, so it will need multiplied by fragment length (or fraction) before
                %     adding to each bin.
                if (val1 == val2)
                    % All of the restriction fragment belongs to one bin.
                    if (val1 <= length(chr_SNPdata_RADseq{chrID,1}))
			chr_SNPdata_RADseq{chrID,1}(val1) = chr_SNPdata_RADseq{chrID,1}(val1)+SNPvalue*fragLength;
			chr_SNPdata_RADseq{chrID,2}(val1) = chr_SNPdata_RADseq{chrID,2}(val1)+fragLength;
                    end;
                else % (val1 < val2)
                    % The restriction fragment belongs partially to two bins, so we must determine fraction assigned to each bin.
                    posEdge     = val1*bases_per_bin;
                    fragLength1 = posEdge-posStart+1;
                    fragLength2 = posEnd-posEdge;
   
                    % Add data to first bin.
                    if (val1 <= length(chr_SNPdata_RADseq{chr,1}))
			chr_SNPdata_RADseq{chrID,1}(val1) = chr_SNPdata_RADseq{chrID,1}(val1)+SNPvalue*fragLength1;
			chr_SNPdata_RADseq{chrID,2}(val1) = chr_SNPdata_RADseq{chrID,2}(val1)+fragLength1;
                    end;
            
                    % Add data to second bin.
                    if (val2 <= length(chr_SNPdata_RADseq{chr,1}))
			chr_SNPdata_RADseq{chrID,1}(val1) = chr_SNPdata_RADseq{chrID,1}(val1)+SNPvalue*fragLength2;
			chr_SNPdata_RADseq{chrID,2}(val1) = chr_SNPdata_RADseq{chrID,2}(val1)+fragLength2;
                    end;
                end;
            end;
        end;
    end;
    fprintf('\n# Fragment corrected SNP ratio values have been added to map bins.');
                
    save([workingDir 'matlab_dir/' child_name '.corrected_SNP_' SNP_verString '.mat'], 'chr_SNPdata_RADseq');
else
    fprintf('\nMAT file containing project SNP information found, loading.\n');
    load([workingDir 'matlab_dir/' child_name '.corrected_SNP_' SNP_verString '.mat']);
end;
chr_SNPdata_RADseq_child = chr_SNPdata_RADseq;
clear chr_SNPdata_RADseq;

if (strcmp(parent_name,child_name) == 0)
    % Initializes vectors used to hold SNP number data.
    for chr = 1:num_chr   % number of chrs.
	for j = 1:2       % 2 categories tracked : total read counts in bin and number of data entries in region (size of bin in base-pairs).
	    chr_SNPdata_RADseq{chr,j} = zeros(1,ceil(chr_size(chr)/bases_per_bin));
	end;
    end;
    SNPfrag = 0;
    %% Gather SNP data into bins for map display.
    if (exist([workingDir 'matlab_dir/' parent_name '.corrected_SNP_' SNP_verString '.mat'],'file') == 0)
	fprintf('\nMAT file containing reference SNP information not found, regenerating from prior data files.\n');

	% fragments_SNP(count).data    = data_bins;           % data bins separated by ":".   first character is "*" for unusable data.

	% Convert bias corrected SNP ratio data per restriction digest fragment into bin-aligned data for plotting
	fprintf('\n# Adding fragment corrected SNP ratio values to map bins.');

	for frag = 1:numFragments
	    if (fragment_data(frag).usable == 1)
		SNPfrag    = SNPfrag+1;
		% Load important data from fragments data structure.
		chrID      = fragment_data(frag).chr;
		posStart   = fragment_data(frag).startbp;
		posEnd     = fragment_data(frag).endbp;
		copyNum    = parent_fragments_SNP(SNPfrag).copyNum;
		data       = parent_fragments_SNP(SNPfrag).data;
		fragLength = posEnd - posStart + 1;
                            
		% Identify locations of fragment end in relation to plotting bins.
		val1 = ceil(posStart/bases_per_bin);
		val2 = ceil(posEnd  /bases_per_bin);

		if (chrID > 0)
		    dataBins = strsplit(data,':');
		    numBins  = length(dataBins);

		    %% Calculate the SNPvalue for the fragment.
		    %       Each putative_SNP contributes to the SNPvalue depending on how close its ratio call is to heterozygous 1:1.
		    %       Homozygous calls contribute 0 per base pair.   Heterozygous 1:1 calls contribute 1 per base pair.
		    %       Intermediate ratios contribute intermediate values.
		    %  The SNPvalue for each fragment is a measure of how heterozygous the fragment is.
		    if (strcmp(dataBins{1},'*') == 1)
                        % This bin has too high a copy number to have sensical ratios.
                        SNPvalue = 0;
                    else
                        % This bin has a copy number between 1 & 8.
                        switch copyNum
                            case 0  % No useful data.
                                SNPvalue = 0;
                            case 1  % No useful data.
                                SNPvalue = 0;
                            case 2  % aa, ab.   
                                % SNPvalue = 0*str2num(dataBins{1}) + 1*str2num(dataBins{2});
                                if (str2num(dataBins{2}) > 0)
                                    SNPvalue = 1;
                                else
                                    SNPvalue = 0;
                                end;
                            case 3  % aaa, aab, abb.
                                % SNPvalue = 0*str2num(dataBins{1}) + 2/3*str2num(dataBins{2});
                                if (str2num(dataBins{2}) > 0)
                                    SNPvalue = 1;
                                else
                                    SNPvalue = 0;   
                                end;
                            case 4  % aaaa, aaab, aabb.
                                % SNPvalue = 0*str2num(dataBins{1}) + 1/2*str2num(dataBins{2}) + 1*str2num(dataBins{3});
                                if (str2num(dataBins{3}) > 0)
                                    SNPvalue = 1;
                                elseif (str2num(dataBins{2}) > 0)
                                    SNPvalue = 1/2;
                                else
                                    SNPvalue = 0;
                                end;
                            case 5  % aaaaa, aaaab, aaabb.
                                % SNPvalue = 0*str2num(dataBins{1}) + 2/5*str2num(dataBins{2}) + 4/5*str2num(dataBins{3});
                                if (str2num(dataBins{3}) > 0)
                                    SNPvalue = 4/5;
                                elseif (str2num(dataBins{2}) > 0) 
                                    SNPvalue = 2/5;
                                else
                                    SNPvalue = 0;
                                end;
                            case 6  % aaaaaa, aaaaab, aaaabb, aaabbb.
                                % SNPvalue = 0*str2num(dataBins{1}) + 1/3*str2num(dataBins{2}) + 2/3*str2num(dataBins{3}) + 1*str2num(dataBins{4});
                                if (str2num(dataBins{4}) > 0)
                                    SNPvalue = 1;
                                elseif (str2num(dataBins{3}) > 0)
                                    SNPvalue = 2/3;
                                elseif (str2num(dataBins{2}) > 0)
                                    SNPvalue = 1/3;
                                else
                                    SNPvalue = 0;
                                end;
                            case 7  % aaaaaaa, aaaaaab, aaaaabb, aaaabbb.
                                % SNPvalue = 0*str2num(dataBins{1}) + 2/7*str2num(dataBins{2}) + 4/7*str2num(dataBins{3}) + 6/7*str2num(dataBins{4});
                                if (str2num(dataBins{4}) > 0)
                                    SNPvalue = 6/7;
                                elseif (str2num(dataBins{3}) > 0)
                                    SNPvalue = 4/7;
                                elseif (str2num(dataBins{2}) > 0)
                                    SNPvalue = 2/7;
                                else
                                    SNPvalue = 0;
                                end;
                            case 8  % aaaaaaaa, aaaaaaab, aaaaaabb, aaaaabbb, aaaabbbb.
                                % SNPvalue = 0*str2num(dataBins{1}) + 1/4*str2num(dataBins{2}) + 1/2*str2num(dataBins{3}) + 3/4*str2num(dataBins{4}) + 1*str2num(dataBins{5});
                                if (str2num(dataBins{5}) > 0)
                                    SNPvalue = 1;
                                elseif (str2num(dataBins{4}) > 0)
                                    SNPvalue = 3/4;
                                elseif (str2num(dataBins{3}) > 0)
                                    SNPvalue = 1/2;
				elseif (str2num(dataBins{2}) > 0)
                                    SNPvalue = 1/4;
                                else
                                    SNPvalue = 0;
                                end;
                        end;
                    end;

		    % 'count' is average read count across fragment, so it will need multiplied by fragment length (or fraction) before
		    %     adding to each bin.
		    if (val1 == val2)
			% All of the restriction fragment belongs to one bin.
			if (val1 <= length(chr_SNPdata_RADseq{chrID,1}))
			    chr_SNPdata_RADseq{chrID,1}(val1) = chr_SNPdata_RADseq{chrID,1}(val1)+SNPvalue*fragLength;
			    chr_SNPdata_RADseq{chrID,2}(val1) = chr_SNPdata_RADseq{chrID,2}(val1)+fragLength;
			end;
		    else % (val1 < val2)
			% The restriction fragment belongs partially to two bins, so we must determine fraction assigned to each bin.
			posEdge     = val1*bases_per_bin;
			fragLength1 = posEdge-posStart+1;
			fragLength2 = posEnd-posEdge;
                            
			% Add data to first bin.
			if (val1 <= length(chr_SNPdata_RADseq{chr,1}))
			    chr_SNPdata_RADseq{chrID,1}(val1) = chr_SNPdata_RADseq{chrID,1}(val1)+SNPvalue*fragLength1;
			    chr_SNPdata_RADseq{chrID,2}(val1) = chr_SNPdata_RADseq{chrID,2}(val1)+fragLength1;
			end;
			% Add data to second bin.
			if (val2 <= length(chr_SNPdata_RADseq{chr,1}))
			    chr_SNPdata_RADseq{chrID,1}(val1) = chr_SNPdata_RADseq{chrID,1}(val1)+SNPvalue*fragLength2;
			    chr_SNPdata_RADseq{chrID,2}(val1) = chr_SNPdata_RADseq{chrID,2}(val1)+fragLength2;
			end;
		    end;
		end;
	    end;
	end;
	fprintf('\n# Fragment corrected SNP ratio values have been added to map bins.');

	save([workingDir 'matlab_dir/' parent_name '.corrected_SNP_' SNP_verString '.mat'], 'chr_SNPdata_RADseq');                
    else
	fprintf('\nMAT file containing reference SNP information found, loading.\n');
	load([workingDir 'matlab_dir/' parent_name '.corrected_SNP_' SNP_verString '.mat']);
    end;
    chr_SNPdata_RADseq_parent = chr_SNPdata_RADseq;
    clear chr_SNPdata_RADseq;
else
    chr_SNPdata_RADseq_parent = chr_SNPdata_RADseq_child;
end;

% basic plot parameters not defined per genome.
TickSize        = -0.005;  %negative for outside, percentage of longest chr figure.
maxY            = 2; % ploidyBase*2;


%% -----------------------------------------------------------------------------------------
% Make figures
%-------------------------------------------------------------------------------------------
fig = figure(1);
set(gcf, 'Position', [0 70 1024 600]);
totalSNPdata = [];
for chr = 1:num_chr
    for pos = 1:length(chr_SNPdata_RADseq_child{chr,1})
	if (chr_SNPdata_RADseq_child{chr,2}(pos) == 0)
	    % No data elements => null value is plotted.
	    SNPplot{chr}(pos) = 0;
	else
	    % Sum of data elements is divided by the number of data elements.
	    SNPplot{chr}(pos) = chr_SNPdata_RADseq_child{chr,1}(pos)/chr_SNPdata_RADseq_child{chr,2}(pos);
	end;
    end;
    totalSNPdata = [totalSNPdata SNPplot{chr}];
end;
max_count    = max(totalSNPdata);
median_count = median(totalSNPdata);
fprintf(['\n max data    = ' num2str(max_count)]);
fprintf(['\n median data = ' num2str(median_count) '\n\n']);

for chr = 1:num_chr
%     SNPplot2{chr} = SNPplot{chr}/median_count;  
    SNPplot2{chr} = SNPplot{chr};
end;
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
    maxY                   = 1; % = ploidyBase*2;
    
    Linear_left = Linear_left_start;
end;

%% -----------------------------------------------------------------------------------------
% Make figures
%-------------------------------------------------------------------------------------------
for chr = 1:num_chr
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
    fprintf(['chr' num2str(chr) ':' num2str(length(SNPplot2{chr})) '\n']);
    for i = 1:length(SNPplot2{chr}); 
        x_ = [i i i-1 i-1];
        SNPhistValue = SNPplot2{chr}(i);
    
        % The SNP-histogram values were normalized to a median value of 1.
        % The ratio of 'ploidy' to 'ploidyBase' determines where the data is displayed relative to the median line.
        startY = 0;
        endY   = SNPhistValue;
        y_ = [startY endY endY startY];
    
        % makes a blackbar for each bin.
        f = fill(x_,y_,c_);
        set(f,'linestyle','none');
    end;
    x2 = chr_size(chr)*chr_length_scale_multiplier;
    
    %% draw lines across plots for easier interpretation of SNP regions.
    switch ploidyBase
	case 1
	    line([0 x2], [maxY/2*1 maxY/2*1],'Color',[0.85 0.85 0.85]);
        case 2
            line([0 x2], [maxY/4*1 maxY/4*1],'Color',[0.85 0.85 0.85]);
	    line([0 x2], [maxY/4*2 maxY/4*2],'Color',[0.85 0.85 0.85]);
            line([0 x2], [maxY/4*3 maxY/4*3],'Color',[0.85 0.85 0.85]);
        case 3
            line([0 x2], [maxY/6*1 maxY/6*1],'Color',[0.85 0.85 0.85]);
            line([0 x2], [maxY/6*2 maxY/6*2],'Color',[0.85 0.85 0.85]);
	    line([0 x2], [maxY/6*3 maxY/6*3],'Color',[0.85 0.85 0.85]);
            line([0 x2], [maxY/6*4 maxY/6*4],'Color',[0.85 0.85 0.85]);
            line([0 x2], [maxY/6*5 maxY/6*5],'Color',[0.85 0.85 0.85]);
        case 4
            line([0 x2], [maxY/8*1 maxY/8*1],'Color',[0.85 0.85 0.85]);
            line([0 x2], [maxY/8*2 maxY/8*2],'Color',[0.85 0.85 0.85]);   
            line([0 x2], [maxY/8*3 maxY/8*3],'Color',[0.85 0.85 0.85]);
	    line([0 x2], [maxY/8*4 maxY/8*4],'Color',[0.85 0.85 0.85]);
            line([0 x2], [maxY/8*5 maxY/8*5],'Color',[0.85 0.85 0.85]);
            line([0 x2], [maxY/8*6 maxY/8*6],'Color',[0.85 0.85 0.85]);
            line([0 x2], [maxY/8*7 maxY/8*7],'Color',[0.85 0.85 0.85]);
    end;
    %% end cgh plot section.
        
    %axes labels etc.
    hold off;
    xlim([0,chr_size(chr)*chr_length_scale_multiplier]);
        
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
            
    % This section sets the Y-axis labelling.
    switch ploidyBase
        case 1
            set(gca,'YTick',[0 maxY/2 maxY]);
            set(gca,'YTickLabel',{'','1','2'});
        case 2
            set(gca,'YTick',[0 maxY/4 maxY/2 maxY/4*3 maxY]);
            set(gca,'YTickLabel',{'','1','2','3','4'});
        case 3
            set(gca,'YTick',[0 maxY/6 maxY/3 maxY/2 maxY/3*2 maxY/6*5 maxY]);
            set(gca,'YTickLabel',{'','','','3','','','6'});
        case 4
            set(gca,'YTick',[0 maxY/8 maxY/4 maxY/8*3 maxY/2 maxY/8*5 maxY/4*3 maxY/8*7 maxY]);
            set(gca,'YTickLabel',{'','','2','','4','','6','','8'});
    end;
    
    set(gca,'FontSize',6);
    if (chr == find(chr_posY == max(chr_posY)))
        title([child_name ' SNP map'],'Interpreter','none','FontSize',12);
    end;
    
    hold on;
    %end axes labels etc.
    
    %show segmental anueploidy breakpoints.
    if (displayBREAKS == true)
	for segment = 2:length(chr_breaks{chr})-1
            bP = chr_breaks{chr}(segment)*length(SNPplot2{chr});
            c_ = [0 0 1];
            x_ = [bP bP bP-1 bP-1];
            y_ = [0 maxY maxY 0];
            f = fill(x_,y_,c_);
            set(f,'linestyle','none');
        end;
    end;
        
    %show centromere.
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
        plot([leftEnd   leftEnd   leftEnd+dx   x1-dx   x1        x2        x2+dx   rightEnd-dx   rightEnd   rightEnd   rightEnd-dx x2+dx   x2   x1   x1-dx   leftEnd+dx   leftEnd],...
             [dy        maxY-dy   maxY         maxY    maxY-dy   maxY-dy   maxY    maxY          maxY-dy    dy         0           0       dy   dy   0       0            dy     ],...
            'Color',[0 0 0]);
    end;
    %end show centromere.
    %show annotation locations
    if (show_annotations) && (length(annotations) > 0)
        plot([leftEnd rightEnd], [-1.5 -1.5],'color',[0 0 0]);
        hold on;
        annotation_location = (annotation_start+annotation_end)./2;
        for i = 1:length(annotation_location)
            if (annotation_chr(i) == chr)
                annotationloc = annotation_location(i)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
		annotationStart = annotation_start(i)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
                annotationEnd   = annotation_end(i)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
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
	    %subplot('Position',[(0.15 + chr_width(chr) + 0.005)+width*(segment-1) bottom width height]);
            subplot('Position',[(left+chr_width(chr)+0.005)+width*(segment-1) bottom width height]);
        
            % The SNP-histogram values were normalized to a median value of 1.
            for i = round(1+length(SNPplot2{chr})*chr_breaks{chr}(segment)):round(length(SNPplot2{chr})*chr_breaks{chr}(segment+1))
                if (Low_quality_ploidy_estimate == true)
                    histAll{segment}(i) = SNPplot2{chr}(i)*ploidy*ploidyAdjust;
                else
                    histAll{segment}(i) = SNPplot2{chr}(i)*ploidy;
                end;
            end;
                     
            % make a histogram of CGH data, then smooth it for display.
            histAll{segment}(histAll{segment}<=0) = [];
            histAll{segment}(histAll{segment}>ploidyBase*2+2)= ploidyBase*2+2;
            histAll{segment}(length(histAll{segment})+1) = 0;   % endpoints added to ensure histogram bounds.
            histAll{segment}(length(histAll{segment})+1) = ploidyBase*2+2;
            smoothed{segment} = smooth_gaussian(hist(histAll{segment},(ploidyBase*2+2)*50),5,20);
            % make a smoothed version of just the endpoints used to ensure histogram bounds.
            histAll2{segment}(1) = 0;
            histAll2{segment}(2) = ploidyBase*2+2;
            smoothed2{segment} = smooth_gaussian(hist(histAll2{segment},(ploidyBase*2+2)*50),5,20)*4;
            % subtract the smoothed endpoints from the histogram to remove the influence of the added endpoints.
            smoothed{segment} = (smoothed{segment}-smoothed2{segment});
            smoothed{segment} = smoothed{segment}/max(smoothed{segment});
        
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
        
        %% cgh plot section. 
        c_ = [0 0 0];
        fprintf(['chr' num2str(chr) ':' num2str(length(SNPplot2{chr})) '\n']);
        for i = 1:length(SNPplot2{chr});
            x_ = [i i i-1 i-1];
            SNPhistValue = SNPplot2{chr}(i);
        
            % The SNP-histogram values were normalized to a median value of 1.
            % The ratio of 'ploidy' to 'ploidyBase' determines where the data is displayed relative to the median line.
            startY = 0;
            endY = SNPhistValue;
            y_ = [startY endY endY startY];
    
            % makes a blackbar for each bin.
            f = fill(x_,y_,c_);
            set(f,'linestyle','none');
        end;
        x2 = chr_size(chr)*chr_length_scale_multiplier;
        
        %% draw lines across plots for easier interpretation of SNP regions.
        switch ploidyBase    
            case 1
		line([0 x2], [maxY/2*1 maxY/2*1],'Color',[0.85 0.85 0.85]);
            case 2
                line([0 x2], [maxY/4*1 maxY/4*1],'Color',[0.85 0.85 0.85]);
		line([0 x2], [maxY/4*2 maxY/4*2],'Color',[0.85 0.85 0.85]);
                line([0 x2], [maxY/4*3 maxY/4*3],'Color',[0.85 0.85 0.85]);
            case 3
                line([0 x2], [maxY/6*1 maxY/6*1],'Color',[0.85 0.85 0.85]);
                line([0 x2], [maxY/6*2 maxY/6*2],'Color',[0.85 0.85 0.85]);
		line([0 x2], [maxY/6*3 maxY/6*3],'Color',[0.85 0.85 0.85]);
                line([0 x2], [maxY/6*4 maxY/6*4],'Color',[0.85 0.85 0.85]);
                line([0 x2], [maxY/6*5 maxY/6*5],'Color',[0.85 0.85 0.85]);
            case 4
                line([0 x2], [maxY/8*1 maxY/8*1],'Color',[0.85 0.85 0.85]);
                line([0 x2], [maxY/8*2 maxY/8*2],'Color',[0.85 0.85 0.85]);   
                line([0 x2], [maxY/8*3 maxY/8*3],'Color',[0.85 0.85 0.85]);
		line([0 x2], [maxY/8*4 maxY/8*4],'Color',[0.85 0.85 0.85]);
                line([0 x2], [maxY/8*5 maxY/8*5],'Color',[0.85 0.85 0.85]);
                line([0 x2], [maxY/8*6 maxY/8*6],'Color',[0.85 0.85 0.85]);
                line([0 x2], [maxY/8*7 maxY/8*7],'Color',[0.85 0.85 0.85]);
        end;
        %% end cgh plot section.

	%show segmental anueploidy breakpoints.
        if (displayBREAKS == true)
            for segment = 2:length(chr_breaks{chr})-1
                bP = chr_breaks{chr}(segment)*length(SNPplot2{chr});
                c_ = [0 0 1];
                x_ = [bP bP bP-1 bP-1];
                y_ = [0 maxY maxY 0];
                f = fill(x_,y_,c_);
                set(f,'linestyle','none');
            end;
        end;
            
        %show centromere.
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
            
            % draw outlines of chromosome cartoon.   (drawn after horizontal lines to that cartoon edges are not interrupted by
            % horiz lines.
	    plot([leftEnd   leftEnd   leftEnd+dx   x1-dx   x1        x2        x2+dx   rightEnd-dx   rightEnd   rightEnd ...
                  rightEnd-dx   x2+dx   x2   x1   x1-dx   leftEnd+dx  leftEnd],...
                 [dy        maxY-dy   maxY         maxY    maxY-dy   maxY-dy   maxY    maxY          maxY-dy    dy       ...
                  0             0       dy   dy   0       0           dy  ],...
                 'Color',[0 0 0]);
        end;
        %end show centromere.
            
        %show annotation locations
        if (show_annotations) && (length(annotations) > 0)
            plot([leftEnd rightEnd], [-1.5 -1.5],'color',[0 0 0]);
            hold on;
            annotation_location = (annotation_start+annotation_end)./2;
            for i = 1:length(annotation_location)
                if (annotation_chr(i) == chr)
                    annotationloc = annotation_location(i)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
		    annotationStart = annotation_start(i)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
                    annotationEnd   = annotation_end(i)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
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
        set(gca,'XTick',0:(40*(5000/bases_per_bin)):(650*(5000/bases_per_bin)));
        set(gca,'XTickLabel',[]);
        if (chr == 1)
            ylabel(child_name, 'Rotation', 0, 'HorizontalAlign', 'right', 'VerticalAlign', 'bottom','Interpreter','none','FontSize',5);
        
            % This section sets the Y-axis labelling.
            switch ploidyBase
                case 1
                    set(gca,'YTick',[0 maxY/2 maxY]);
                    set(gca,'YTickLabel',{'','1','2'});
                case 2
                    set(gca,'YTick',[0 maxY/4 maxY/2 maxY/4*3 maxY]);
                    set(gca,'YTickLabel',{'','1','2','3','4'});
                case 3
                    set(gca,'YTick',[0 maxY/6 maxY/3 maxY/2 maxY/3*2 maxY/6*5 maxY]);
                    set(gca,'YTickLabel',{'','','','3','','','6'});
                case 4
                    set(gca,'YTick',[0 maxY/8 maxY/4 maxY/8*3 maxY/2 maxY/8*5 maxY/4*3 maxY/8*7 maxY]);
                    set(gca,'YTickLabel',{'','','2','','4','','6','','8'});
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

% Save original arrangement of chromosomes.
saveas(fig, [figureDir child_name '.corrected-SNP-map.1.eps'], 'epsc');
delete(fig);
        
% Save horizontally arranged chromosomes.
set(Linear_fig,'PaperPosition',[0 0 8 0.62222222]);
saveas(Linear_fig, [figureDir child_name '.corrected-SNP-map.2.eps'], 'epsc');
delete(Linear_fig);
        










%{
% basic plot parameters not defined per genome.
TickSize        = -0.005;  %negative for outside, percentage of longest chr figure.
maxY            = 1; % = ploidyBase*2;

%define colors for colorBars plot
colorNoData = [1.0   1.0   1.0  ]; %used when no data is available for the bin.
colorInit   = [0.5   0.5   0.5  ]; %external; used in blending at ends of chr.
colorHET    = [0.0   0.0   0.0  ]; % near 1:1 ratio SNPs
colorOddHET = [0.0   1.0   0.0  ]; % Het, but not near 1:1 ratio SNPs.
colorHOM    = [1.0   0.0   0.0  ]; % Hom SNPs;

%% -----------------------------------------------------------------------------------------
% Setup for main figure generation.
%------------------------------------------------------------------------------------------
fig = figure(1);
set(gcf, 'Position', [0 70 1024 600]);
data_mode = 3;
for chr = 1:num_chr
    if (data_mode == 1)
        % Regenerate chr plot data if the save file does not exist.
        TOTplot{chr}                           = chr_SNPdata{chr,1}+chr_SNPdata{chr,2}+chr_SNPdata{chr,3};  % TOT data
        TOTave{chr}                            = sum(TOTplot{chr})/length(TOTplot{chr});
        TOTplot2{chr}                          = TOTplot{chr}/ave_copy_num;
        TOTplot2{chr}(TOTplot2{chr} > 1)       = 1;
        TOTave2{chr}                           = sum(TOTplot2{chr})/length(TOTplot2{chr});
    
        HETplot{chr}                           = chr_SNPdata{chr,1};  % HET data
        HETave{chr}                            = sum(HETplot{chr})/length(HETplot{chr});
        HETplot2{chr}                          = HETplot{chr}/ave_copy_num;
        HETplot2{chr}(HETplot2{chr} > 1)       = 1;
        HETave2{chr}                           = sum(HETplot2{chr})/length(HETplot2{chr});
    
        oddHETplot{chr}                        = chr_SNPdata{chr,2};  % oddHET data
        oddHETave{chr}                         = sum(oddHETplot{chr})/length(oddHETplot{chr});
        oddHETplot2{chr}                       = oddHETplot{chr}/ave_copy_num;
        oddHETplot2{chr}(oddHETplot2{chr} > 1) = 1;
        oddHETave2{chr}                        = sum(oddHETplot2{chr})/length(oddHETplot2{chr});
    
        HOMplot{chr}                           = chr_SNPdata{chr,3};  % HOM data
        HOMave{chr}                            = sum(HOMplot{chr})/length(HOMplot{chr});
        HOMplot2{chr}                          = HOMplot{chr}/ave_copy_num;
        HOMplot2{chr}(HOMplot2{chr} > 1)       = 1;
        HOMave2{chr}                           = sum(HOMplot2{chr})/length(HOMplot2{chr});
    elseif (data_mode == 2)
        %% Details from LOH_v2a.m :
        % Regenerate chr plot data if the save file does not exist.
        TOTplot{chr}                                  = chr_SNPdata{chr,1}+chr_SNPdata{chr,2}+chr_SNPdata{chr,3};  % TOT data
        HETplot{chr}                                  = chr_SNPdata{chr,1};  % HET data
        oddHETplot{chr}                               = chr_SNPdata{chr,2};  % oddHET data
        HOMplot{chr}                                  = chr_SNPdata{chr,3};  % HOM data

        TOTave{chr}                                   = sum(TOTplot{chr})/length(TOTplot{chr});
        TOTplot2{chr}                                 = TOTplot{chr}/ave_copy_num;
        TOTplot2{chr}(TOTplot2{chr} > 1)              = 1;
        TOTave2{chr}                                  = sum(TOTplot2{chr})/length(TOTplot2{chr});
    
        HETave{chr}                                   = sum(HETplot{chr})/length(HETplot{chr});
        HETplot2{chr}                                 = HETplot{chr}/ave_copy_num;
        HETplot2{chr}(HETplot2{chr} > 1)              = 1;
        HETave2{chr}                                  = sum(HETplot2{chr})/length(HETplot2{chr});
        oddHETave{chr}                                = sum(oddHETplot{chr})/length(oddHETplot{chr});
        oddHETplot2{chr}                              = oddHETplot{chr}/ave_copy_num;
        oddHETplot2{chr}(oddHETplot2{chr} > 1)        = 1;
        oddHETave2{chr}                               = sum(oddHETplot2{chr})/length(oddHETplot2{chr});
        HOMave{chr}                                   = sum(HOMplot{chr})/length(HOMplot{chr});
        HOMplot2{chr}                                 = HOMplot{chr}/ave_copy_num;
        HOMplot2{chr}(HOMplot2{chr} > 1)              = 1;
        HOMave2{chr}                                  = sum(HOMplot2{chr})/length(HOMplot2{chr});
    elseif (data_mode == 3)
        %% Details from LOH_v3a.m :
        % Regenerate chr plot data if the save file does not exist.
        TOTplot{chr}                                  = chr_SNPdata{chr,1}+chr_SNPdata{chr,2}+chr_SNPdata{chr,3};  % TOT data
        TOTave{chr}                                   = sum(TOTplot{chr})/length(TOTplot{chr});
        TOTplot2{chr}                                 = TOTplot{chr}/ave_copy_num;
        TOTplot2{chr}(TOTplot2{chr} > 1)              = 1;
        TOTave2{chr}                                  = sum(TOTplot2{chr})/length(TOTplot2{chr});
    
        HETplot{chr}                                  = chr_SNPdata{chr,1};  % HET data
        HETave{chr}                                   = sum(HETplot{chr})/length(HETplot{chr});
        HETplot2{chr}                                 = HETplot{chr}/ave_copy_num;
        HETplot2{chr}(HETplot2{chr} > 1)              = 1;
        HETave2{chr}                                  = sum(HETplot2{chr})/length(HETplot2{chr});
    
        oddHETplot{chr}                               = chr_SNPdata{chr,2};  % oddHET data
        oddHETave{chr}                                = sum(oddHETplot{chr})/length(oddHETplot{chr});
        oddHETplot2{chr}                              = oddHETplot{chr}/ave_copy_num;
        oddHETplot2{chr}(oddHETplot2{chr} > 1)        = 1;
    
        HOMplot{chr}                                  = chr_SNPdata{chr,3};  % HOM data
        HOMave{chr}                                   = sum(HOMplot{chr})/length(HOMplot{chr});
        HOMplot2{chr}                                 = HOMplot{chr}/ave_copy_num;
        HOMplot2{chr}(HOMplot2{chr} > 1)              = 1;
    end;
end;
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
    maxY                   = ploidyBase*2;

    Linear_left = Linear_left_start;
end;

%% -----------------------------------------------------------------------------------------
% Make figures
%-------------------------------------------------------------------------------------------
for chr = 1:num_chr
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
    for i = 1:length(TOTplot2{chr})+1;
        if (i-1 < length(TOTplot2{chr}))
            c_tot_post = TOTplot2{chr}(i)+TOTplot2{chr}(i);
            if (c_tot_post == 0)
                c_post = colorNoData;
            else
                %c_post =   colorHET*HETplot2{chr}(i) + ...
                %           colorHOM*HOMplot2{chr}(i) + ...
                %           colorNoData*(1-min([HETplot2{chr}(i)+HOMplot2{chr}(i) 1]));
                %colorMix = colorHET   *HETplot2   {chr}(i)/TOTplot2{chr}(i) + ...
                %           colorOddHET*oddHETplot2{chr}(i)/TOTplot2{chr}(i) + ...
                %           colorHOM   *HOMplot2   {chr}(i)/TOTplot2{chr}(i);
                colorMix = colorHET   *   HETplot2{chr}(i)/TOTplot2{chr}(i) + ...
                           colorOddHET*oddHETplot2{chr}(i)/TOTplot2{chr}(i) + ...
                           colorHOM   *   HOMplot2{chr}(i)/TOTplot2{chr}(i);
                c_post =   colorMix   *   min(1,TOTplot2{chr}(i)) + ...
                           colorNoData*(1-min(1,TOTplot2{chr}(i)));
                %colorNoData*(1-min([HETplot2{chr}(i)+oddHETplot2{chr}(i)+HOMplot2{chr}(i) 1]));
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
    xlim([0,chr_size(chr)*chr_length_scale_multiplier]);
    
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
    set(gca,'YTick',[maxY/4 maxY/2 maxY/4*3 maxY]);
    set(gca,'YTickLabel',{'1','2','3','4'});
    set(gca,'FontSize',6);
    if (chr == find(chr_posY == max(chr_posY)))
        if (strcmp(parent_name,child_name) == 1)
            title([ parent_name ' SNP map'],'Interpreter','none','FontSize',12);
        else
            title([ parent_name ' vs. ' child_name ' SNP/LOH map'],'Interpreter','none','FontSize',12);
        end;
    end;
    hold on;
    %end axes labels etc.
    
    %show centromere outlines and horizontal marks.
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
    end;
    %end show centromere.
    
    %show annotation locations
    if (show_annotations) && (length(annotations) > 0)
        plot([leftEnd rightEnd], [-1.5 -1.5],'color',[0 0 0]);
        hold on;
        annotation_location = (annotation_start+annotation_end)./2;
        for i = 1:length(annotation_location)
            if (annotation_chr(i) == chr)
                annotationloc = annotation_location(i)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
                plot(annotationloc,-1.5,'k:o','MarkerEdgeColor',annotation_edgecolor{i},'MarkerFaceColor',annotation_fillcolor{i},'MarkerSize',annotation_size(i));
                % plot(annotationloc,-1.5,'k:o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);
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
        title(chr_label{chr},'Interpreter','none','FontSize',10);

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
                bP = chr_breaks{chr}(segment)*length(SNPplot2{chr});
                c_ = [0 0 1];
                x_ = [bP bP bP-1 bP-1];
                y_ = [0 maxY maxY 0];
                f = fill(x_,y_,c_);   
                set(f,'linestyle','none');
            end;
        end;

        %show centromere.
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
                  [dy        maxY-dy   maxY         maxY    maxY-dy   maxY-dy   maxY    maxY          maxY-dy    dy         0             0       dy   dy   0       0            dy],...
                  'Color',[0 0 0]);
        end;
        %end show centromere.

        %show annotation locations
        if (show_annotations) && (length(annotations) > 0)
            plot([leftEnd rightEnd], [-1.5 -1.5],'color',[0 0 0]);
            hold on;
            annotation_location = (annotation_start+annotation_end)./2;
            for i = 1:length(annotation_location)
                if (annotation_chr(i) == chr)
                    annotationloc = annotation_location(i)*chr_length_scale_multiplier-0.5*(5000/bases_per_bin);
                    plot(annotationloc,-1.5,'k:o','MarkerEdgeColor',annotation_edgecolor{i},'MarkerFaceColor',annotation_fillcolor{i},'MarkerSize',annotation_size(i));
                    % plot(annotationloc,-1.5,'k:o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);
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
            if (strcmp(parent_name,child_name) == 1)
                ylabel(parent_name, 'Rotation', 0, 'HorizontalAlign', 'right', 'VerticalAlign', 'bottom','Interpreter','none','FontSize',5);
            else
                ylabel({parent_name;'vs.';child_name}, 'Rotation', 0, 'HorizontalAlign', 'right', 'VerticalAlign', 'bottom','Interpreter','none','FontSize',5);
            end;
            set(gca,'YTick',[maxY/4 maxY/2 maxY/4*3 maxY]);
            set(gca,'YTickLabel',{'1','2','3','4'});
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
if (strcmp(parent_name,child_name) == 1)
    saveas(fig,        [figureDir parent_name '.SNP-map.1.eps'], 'epsc');
    set(Linear_fig,'PaperPosition',[0 0 8 0.62222222]);
    saveas(Linear_fig, [figureDir parent_name '.SNP-map.2.eps'], 'epsc');
else
    saveas(fig,        [figureDir parent_name '->' child_name '.SNP-map.1.eps'], 'epsc');
    set(Linear_fig,'PaperPosition',[0 0 8 0.62222222]);
    saveas(Linear_fig, [figureDir parent_name '->' child_name '.SNP-map.2.eps'], 'epsc');
end;
delete(fig);
delete(Linear_fig);
%}

%% ========================================================================
% end stuff
%==========================================================================
end
