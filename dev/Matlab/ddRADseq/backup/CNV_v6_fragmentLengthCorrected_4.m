function [] = CNV_v6_fragmentLengthCorrected_4(projectName,genome_data,referenceName,genome_ref,ploidyString,ploidyBaseString, ...
                                               CNV_verString,rDNA_verString,workingDir,figureDir,displayBREAKS, referenceCHR)
%% ========================================================================
% Generate CGH-type figures from RADseq data, using a reference dataset to correct for genome position-dependant biases.
%==========================================================================
Centromere_format_default   = 0;
Yscale_nearest_even_ploidy  = true;
HistPlot                    = true;
ChrNum                      = true;
show_annotations            = true;
Linear_display              = true;
Low_quality_ploidy_estimate = true;

%%=========================================================================
% Control variables.
%--------------------------------------------------------------------------
% Defines chr sizes in bp. (diploid total=28,567,7888)
% Defines centromere locations in bp.
% Defines annotation locations in bp.
[centromeres, chr_sizes, figure_details, annotations, ploidy_default] = Load_genome_information_1(workingDir,figureDir,genome_data);
[Aneuploidy]                                                          = Load_dataset_information_1(projectName,workingDir);

for i = 1:length(chr_sizes)
    chr_size(chr_sizes(i).chr)    = chr_sizes(i).size;
end;
for i = 1:length(centromeres)
    cen_start(centromeres(i).chr) = centromeres(i).start;
    cen_end(centromeres(i).chr)   = centromeres(i).end;
end;
if (length(annotations) > 0)
    fprintf(['\nAnnotations for ' genome_data '.\n']);
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
	chr_in_use(figure_details(i).chr) = str2num(figure_details(i).useChr);
    end;
end;
num_chr = length(chr_size);
% bases_per_bin			= 5000;
bases_per_bin			= max(chr_size)/700;
chr_length_scale_multiplier	= 1/bases_per_bin;

%%=========================================================================
%%= No further control variables below. ===================================
%%=========================================================================


%-------------------------------------------------------------------------------------------------
% Sanitize user input of euploid state.
%-------------------------------------------------------------------------------------------------
ploidyBase = round(str2num(ploidyBaseString));
if (ploidyBase > 4);   ploidyBase = 4;   end;
if (ploidyBase < 1);   ploidyBase = 1;   end;
fprintf(['\nEuploid base = "' num2str(ploidyBase) '"\n']);


%-------------------------------------------------------------------------------------------------
% Determine reference genome FASTA file in use.
%-------------------------------------------------------------------------------------------------
%  Read in and parse : "links_dir/main_script_dir/genome_specific/[genome]/reference.txt"
reference_file   = [workingDir 'main_script_dir/genomeSpecific/' genome_data '/reference.txt'];
refernce_fid     = fopen(reference_file, 'r');
refFASTA         = fgetl(refernce_fid);
fclose(refernce_fid);


%-------------------------------------------------------------------------------------------------
% Load pre-processed ddRADseq fragment copy number data for project.
%-------------------------------------------------------------------------------------------------
if (exist([workingDir 'matlab_dir/' projectName '.fragment_CNV_data.mat'],'file') == 0)
    %  Including : [chrNum, bpStart, bpEnd, maxReads, AveReads, fragmentLength]
    %	1	9638	10115	2	0	478
    %	1	10116	10123	2	1	8
    %	1	13170	13841	0	0	672
    fprintf('Loading results from Python script, which pre-processed the dataset relative to genome restriction fragments : project\n');
    datafile_RADseq  = [workingDir 'pileup_dir/' projectName '_RADseq_digest_analysis_CNV.txt'];
    data_RADseq      = fopen(datafile_RADseq);
    count            = 0;
    fragments_CNV    = [];
    while ~feof(data_RADseq)
	% Load fragment data from pre-processed text file, single line.
	tline = fgetl(data_RADseq);

	% check if line is a comment.
	test = sscanf(tline, '%s',1);
	if (strcmp(test,'###') == 0)   % If the first space-delimited string is not '###', then process it as a valid data line.
	    % The number of valid lines found so far...  the number of usable restriction fragments with data so far.
	    count = count + 1;

	    % Each valid data line consists of six tab-delimited collumns. 
	    % [chrNum, bpStart, bpEnd, Max, Ave, Length]

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

	    % highest read count along fragment.
	    new_string = sscanf(tline, '%s',  4 );
	    for i = 1:size(sscanf(tline,'%s', 3 ),2);   new_string(1) = [];   end;
	    data_max = new_string;

	    % average read count along fragment.
	    new_string = sscanf(tline, '%s',  5 );
	    for i = 1:size(sscanf(tline,'%s', 4 ),2);   new_string(1) = [];   end;
	    data_ave  = new_string;

	    % length of fragment in bp.
	    new_string = sscanf(tline, '%s',  6 );
	    for i = 1:size(sscanf(tline,'%s', 5 ),2);   new_string(1) = [];   end;
	    fragment_length = new_string;

	    % Add fragment data to data structure.
	    fragments_CNV(count).chr        = str2num(chr_num);
	    fragments_CNV(count).startbp    = str2num(bp_start);
	    fragments_CNV(count).endbp      = str2num(bp_end);
	    fragments_CNV(count).length     = str2num(fragment_length);
	    fragments_CNV(count).read_count = 0;
	    fragments_CNV(count).read_max   = str2num(data_max);
	    fragments_CNV(count).read_ave   = str2num(data_ave);
	    fragments_CNV(count).usable     = 1;
	end;
    end;
    fclose(data_RADseq);
    save([workingDir 'matlab_dir/' projectName '.fragment_CNV_data.mat'], 'fragments_CNV');
else
    load([workingDir 'matlab_dir/' projectName '.fragment_CNV_data.mat']);
end;
project_fragments_CNV = fragments_CNV;
clear fragments_CNV;

%-------------------------------------------------------------------------------------------------
% Load pre-processed ddRADseq fragment copy number data for reference.
%-------------------------------------------------------------------------------------------------
if (exist([workingDir 'matlab_dir/' referenceName '.fragment_CNV_data.mat'],'file') == 0)
    %  Including : [chrNum, bpStart, bpEnd, maxReads, AveReads, fragmentLength]
    %   1       9638    10115   2       0       478
    %   1       10116   10123   2       1       8
    %   1       13170   13841   0       0       672
    fprintf('Loading results from Python script, which pre-processed the dataset relative to genome restriction fragments : reference\n');
    datafile_RADseq   = [workingDir 'pileup_dir/' referenceName '_RADseq_digest_analysis_CNV.txt'];
    data_RADseq       = fopen(datafile_RADseq);
    count             = 0;
    fragments_CNV     = [];
    while ~feof(data_RADseq)
        % Load fragment data from pre-processed text file, single line.
        tline = fgetl(data_RADseq);

        % check if line is a comment.
        test = sscanf(tline, '%s',1);
        if (strcmp(test,'###') == 0)   % If the first space-delimited string is not '###', then process it as a valid data line.
            % The number of valid lines found so far...  the number of usable restriction fragments with data so far.
            count = count + 1;

            % Each valid data line consists of six tab-delimited collumns.
            % [chrNum, bpStart, bpEnd, Max, Ave, Length]

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

            % highest read count along fragment.
            new_string = sscanf(tline, '%s',  4 );
            for i = 1:size(sscanf(tline,'%s', 3 ),2);   new_string(1) = [];   end;
            data_max = new_string;

            % average read count along fragment.
            new_string = sscanf(tline, '%s',  5 );
            for i = 1:size(sscanf(tline,'%s', 4 ),2);   new_string(1) = [];   end;
            data_ave  = new_string;

            % length of fragment in bp.
            new_string = sscanf(tline, '%s',  6 );
            for i = 1:size(sscanf(tline,'%s', 5 ),2);   new_string(1) = [];   end;
            fragment_length = new_string;

            % Add fragment data to data structure.
            fragments_CNV(count).chr        = str2num(chr_num);
            fragments_CNV(count).startbp    = str2num(bp_start);
            fragments_CNV(count).endbp      = str2num(bp_end);
            fragments_CNV(count).length     = str2num(fragment_length);
            fragments_CNV(count).read_count = 0;
            fragments_CNV(count).read_max   = str2num(data_max);
            fragments_CNV(count).read_ave   = str2num(data_ave);
            fragments_CNV(count).usable     = 1;
        end;
    end;
    fclose(data_RADseq);
    save([workingDir 'matlab_dir/' referenceName '.fragment_CNV_data.mat'], 'fragments_CNV');
else
    load([workingDir 'matlab_dir/' referenceName '.fragment_CNV_data.mat']);
end;
reference_fragments_CNV = fragments_CNV;
clear fragments_CNV;


%-------------------------------------------------------------------------------------------------
% Load pre-processed ddRADseq fragment repetitiveness data for genome.
%-------------------------------------------------------------------------------------------------
if (exist([workingDir 'main_script_dir/repetitiveness_files/' genome_data  '.fragment_repet_data.mat'],'file') == 0)
    %  Including : [chrNum, bpStart, bpEnd, maxRepet, AveRepet, fragmentLength]
    %	1	9638	10115	4796	731	478
    %	1	10116	10123	1283	1082	8
    %	1	13170	13841	55175	1626	672
    fprintf('Loading repetitiveness per restriction fragment.\n');
    datafile_RADseq  = [workingDir 'main_script_dir/repetitiveness_files/' genome_data '_repetitiveness.MfeI_MboI.txt'];
    data_RADseq      = fopen(datafile_RADseq);
    count            = 0;
    fragments_repet  = [];
    while ~feof(data_RADseq)
	% Load fragment data from pre-processed text file, single line.
	tline = fgetl(data_RADseq);

	% check if line is a comment.
	test = sscanf(tline, '%s',1);
	if (strcmp(test,'###') == 0)   % If the first space-delimited string is not '###', then process it as a valid data line.
	    % The number of valid lines found so far...  the number of usable restriction fragments with data so far.
	    count = count + 1;
        
	    % Each valid data line consists of six tab-delimited collumns.
	    % [chrNum, bpStart, bpEnd, Max, Ave, Length]
        
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
    
	    % highest read count along fragment.
	    new_string = sscanf(tline, '%s',  4 );
	    for i = 1:size(sscanf(tline,'%s', 3 ),2);   new_string(1) = [];   end;
	    data_max = new_string;
        
	    % average read count along fragment.
	    new_string = sscanf(tline, '%s',  5 );
	    for i = 1:size(sscanf(tline,'%s', 4 ),2);   new_string(1) = [];   end;
	    data_ave  = new_string;
        
	    % length of fragment in bp.
	    new_string = sscanf(tline, '%s',  6 );
	    for i = 1:size(sscanf(tline,'%s', 5 ),2);   new_string(1) = [];   end;
	    fragment_length = new_string;

	    % Add fragment data to data structure.
	    fragments_repet(count).chr         = str2num(chr_num);
	    fragments_repet(count).startbp     = str2num(bp_start);
	    fragments_repet(count).endbp       = str2num(bp_end);
	    fragments_repet(count).length      = str2num(fragment_length);
	    fragments_repet(count).repet_count = 0;
	    fragments_repet(count).repet_max   = str2num(data_max);
	    fragments_repet(count).repet_ave   = str2num(data_ave);
	    fragments_repet(count).usable      = 1;
	end;
    end;
    fclose(data_RADseq);
    save([workingDir 'main_script_dir/repetitiveness_files/' genome_data  '.fragment_repet_data.mat'], 'fragments_repet');
else
    load([workingDir 'main_script_dir/repetitiveness_files/' genome_data  '.fragment_repet_data.mat']);
end;


%-------------------------------------------------------------------------------------------------
% Standardize 'fragment_data' data structure to hold results.
%-------------------------------------------------------------------------------------------------
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

%% Figuring out which fragments have usable data.
usable_list = zeros(1,numFragments);
fprintf('Gathering [fragment length], [CNV], & [repetitiveness] data for plotting bias and correction.\n');
for frag = 1:numFragments
    if (fragment_data(frag).length == 0) || (fragment_data(frag).repet_ave == 0) || (fragment_data(frag).read_ave == 0) || (fragment_data(frag).ref_read_ave == 0)
	fragment_data(frag).usable = 0;
    else
	usable_list(frag)          = 1;
    end;
end;

%% clean up empty entries for fitting process.
X_length    = zeros(1,numFragments);
X_repet     = zeros(1,numFragments);
Y_reads     = zeros(1,numFragments);
Y_reads_ref = zeros(1,numFragments);
for frag = 1:numFragments
    X_length(frag)    = fragment_data(frag).length; 
    X_repet(frag)     = fragment_data(frag).repet_ave;
    Y_reads(frag)     = fragment_data(frag).read_ave;
    Y_reads_ref(frag) = fragment_data(frag).ref_read_ave;
end;
fprintf(['\tnumFragments              = ' num2str(numFragments)              '\n']);
fprintf(['\tlength(fragment_data)     = ' num2str(length(fragment_data))     '\n']);
fprintf(['\tlength(X_length)          = ' num2str(length(X_length))          '\n']);
fprintf(['\tlength(Y_reads)           = ' num2str(length(Y_reads))           '\n']);
fprintf(['\tlength(Y_reads_ref)       = ' num2str(length(Y_reads_ref))       '\n']);
fprintf(['\tlength(X_repet)           = ' num2str(length(X_repet))           '\n']);
% trim datasets.
Y_reads    (X_length > 1000) = [];
Y_reads_ref(X_length > 1000) = [];
X_repet    (X_length > 1000) = [];
X_length   (X_length > 1000) = [];
fprintf('\n');
fprintf(['\tlength(X_length)          = ' num2str(length(X_length))          '\n']);
fprintf(['\tlength(Y_reads)           = ' num2str(length(Y_reads))           '\n']);
fprintf(['\tlength(Y_reads_ref)       = ' num2str(length(Y_reads_ref))       '\n']);
fprintf(['\tlength(X_repet)           = ' num2str(length(X_repet))           '\n']);
fprintf('\n');

% subplot handle list.
sh=zeros(6,1);

%-------------------------------------------------------------------------------------------------
% Make figure illustrating CNV correction steps.
%-------------------------------------------------------------------------------------------------
fprintf('Generating figure with plot of length vs. coverage data.\n');
fig = figure(1);

fprintf('Subplot 1/6 : [EXPERIMENT] (Ave read count) vs. (Fragment length).\n');
sh(1) = subplot(5,4,[1 4]);
if (exist([workingDir 'matlab_dir/' projectName '.Lowess1.mat'],'file') == 0)
    fprintf('\tLOWESS fitting to project data.\n');
    [newX1, newY1] = optimize_mylowess(X_length,Y_reads,10);
    fprintf('\tLOWESS fitting to project data complete.\n');
    save([workingDir 'matlab_dir/' projectName '.Lowess1.mat'], 'newX1','newY1');
else
    fprintf('\tLOWESS fitting to project data loading.\n');
    load([workingDir 'matlab_dir/' projectName '.Lowess1.mat']);
    fprintf('\tLOWESS fitting to project data loaded.\n');
end;
plot(X_length,Y_reads,'.', 'color', [0.0, 0.66667, 0.33333], 'MarkerSize', 1);
hold on;
plot(newX1, newY1,'color','k','linestyle','-', 'linewidth',1);
hold off;
axis normal;
h = title('Examine bias between ddRADseq fragment length and read count. [EXPERIMENT]');   set(h, 'FontSize', 5);
h = ylabel('Average Reads');     set(h, 'FontSize', 5);
h = xlabel('Fragment Length');   set(h, 'FontSize', 5);
set(gca,'FontSize',5);
ylim([0 800]);
xlim([1 1000]);
newX1_proj = newX1;
newY1_proj = newY1;

fprintf('Subplot 2/6 : [REFERENCE] (Ave read count) vs. (Fragment length).\n');
sh(2) = subplot(5,4,[5 8]);
if (exist([workingDir 'matlab_dir/' referenceName '.Lowess1.mat'],'file') == 0)   
    fprintf('\tLOWESS fitting to reference data.\n');
    [newX1, newY1] = optimize_mylowess(X_length,Y_reads_ref,10);
    fprintf('\tLOWESS fitting to referemce data complete.\n');
    save([workingDir 'matlab_dir/' referenceName '.Lowess1.mat'], 'newX1','newY1');
else
    fprintf('\tLOWESS fitting to reference data loading.\n');
    load([workingDir 'matlab_dir/' referenceName '.Lowess1.mat']);
    fprintf('\tLOWESS fitting to reference data loaded.\n');
end;
plot(X_length,Y_reads_ref,'.', 'color', [0.0, 0.66667, 0.33333], 'MarkerSize', 1);
hold on;
plot(newX1, newY1,'color','k','linestyle','-', 'linewidth',1);
hold off;
axis normal;
h = title('Examine bias between ddRADseq fragment length and read count. [REFERENCE]');   set(h, 'FontSize', 5);
h = ylabel('Average Reads');     set(h, 'FontSize', 5);
h = xlabel('Fragment Length');   set(h, 'FontSize', 5);
set(gca,'FontSize',5);
ylim([0 800]);
xlim([1 1000]);
newX1_ref = newX1;
newY1_ref = newY1;

%% Redefine data vectors.
X_length    = zeros(1,numFragments);
X_repet     = zeros(1,numFragments);
Y_reads     = zeros(1,numFragments);
Y_reads_ref = zeros(1,numFragments);
for frag = 1:numFragments
    X_length(frag)    = fragment_data(frag).length;
    X_repet(frag)     = fragment_data(frag).repet_ave;
    Y_reads(frag)     = fragment_data(frag).read_ave;
    Y_reads_ref(frag) = fragment_data(frag).ref_read_ave;
end;

fprintf('Subplot 3/6 : [EXPERIMENT] (Ave read count, with bias corrected) vs. (Fragment length).\n');
sh(3) = subplot(5,4,[9 12]);
hold on;
Y_target = 1;
% Calculate bias_corrected ave_read_copy data for plotting and later analysis.
Y_fitCurve1       = interp1(newX1_proj,newY1_proj,X_length,'spline');  
Y_reads_corrected = Y_reads./Y_fitCurve1*Y_target;
fprintf(['\tlength(Y_fitCurve1)       = ' num2str(length(Y_fitCurve1)) '\n']);
fprintf(['\tlength(Y_reads_corrected) = ' num2str(length(Y_reads_corrected)) '\n']);
fprintf('\t');
for frag = 1:numFragments
    % Add corrected data_ave to fragment data structure.
    fragment_data(frag).reads_ave_corrected = Y_reads_corrected(frag);

    % Define data as no useful if correction fit term falls below 5 reads.
    if (Y_fitCurve1(frag)       <  1   );   fragment_data(frag).usable = 0;   end;
    if (Y_reads_corrected(frag) <= 0   );   fragment_data(frag).usable = 0;   end;
    if (X_length(frag)          <  50  );   fragment_data(frag).usable = 0;   end;
    if (X_length(frag)          >  1000);   fragment_data(frag).usable = 0;   end;

    X_datum = fragment_data(frag).length;			% X_Length
    Y_datum = fragment_data(frag).reads_ave_corrected;		% Y_reads_corrected

    % Plot each corrected data point, colored depending on if data is useful or not.
    if (fragment_data(frag).usable == 1)
	plot(X_datum,Y_datum,'.', 'color', [0.0, 0.66667, 0.33333], 'MarkerSize',0.2);
    else
        plot(X_datum,Y_datum,'.', 'color', [1.0, 0.0, 0.0], 'MarkerSize',0.2);
    end;
    if (mod(frag,1000) == 0)
	fprintf('.');
    end;
end;
fprintf('\n');
plot(newX1,Y_target,'color','k','linestyle','-', 'linewidth',1);
hold off;
axis normal;
h = title('Corrected Read count vs. ddRADseq fragment length. [EXPERIMENT]');   set(h, 'FontSize', 5);
h = ylabel('Corrected Reads');   set(h, 'FontSize', 5);
h = xlabel('Fragment Length');   set(h, 'FontSize', 5);
set(gca,'FontSize',5);
ylim([0 4]);
xlim([1 1000]);


fprintf('Subplot 4/6 : [REFERENCE] (Ave read count, with bias corrected) vs. (Fragment length).\n');
sh(4) = subplot(5,4,[13 16]);
hold on;
Y_target = 1;
% Calculate bias_corrected ave_read_copy data for plotting and later analysis.
Y_fitCurve1_ref       = interp1(newX1_ref,newY1_ref,X_length,'spline');   
Y_reads_corrected_ref = Y_reads_ref./Y_fitCurve1_ref*Y_target;
fprintf(['\tlength(Y_fitCurve1_ref)       = ' num2str(length(Y_fitCurve1_ref))       '\n']);
fprintf(['\tlength(Y_reads_corrected_ref) = ' num2str(length(Y_reads_corrected_ref)) '\n']);
fprintf('\t');
for frag = 1:numFragments
    % Add corrected data_ave to fragment data structure.
    fragment_data(frag).reads_ave_corrected_ref = Y_reads_corrected_ref(frag);

    % Define data as no useful if correction fit term falls below 5 reads.
    fragment_data(frag).usable_ref = 1;
    if (Y_fitCurve1_ref(frag)       <  1   );   fragment_data(frag).usable_ref = 0;   end;
    if (Y_reads_corrected_ref(frag) == 0   );   fragment_data(frag).usable_ref = 0;   end;
    if (X_length(frag)              <  50  );   fragment_data(frag).usable_ref = 0;   end;
    if (X_length(frag)              >  1000);   fragment_data(frag).usable_ref = 0;   end;

    X_datum = fragment_data(frag).length;                    % X_Length
    Y_datum = fragment_data(frag).reads_ave_corrected_ref;   % Y_reads_corrected_ref

    % Plot each corrected data point, colored depending on if data is useful or not.
    if (fragment_data(frag).usable_ref == 1)
	plot(X_datum,Y_datum,'.', 'color', [0.0, 0.66667, 0.33333], 'MarkerSize',0.2);
    else
        plot(X_datum,Y_datum,'.', 'color', [1.0, 0.0, 0.0], 'MarkerSize',0.2);
    end;
    if (mod(frag,1000) == 0)
        fprintf('.');
    end;
end;
fprintf('\n');
plot(newX1,Y_target,'color','k','linestyle','-', 'linewidth',1);
hold off;
axis normal;
h = title('Corrected Read count vs. ddRADseq fragment length. [REFERENCE]');   set(h, 'FontSize', 5);
h = ylabel('Corrected Reads');   set(h, 'FontSize', 5);
h = xlabel('Fragment Length');   set(h, 'FontSize', 5);
set(gca,'FontSize',5);
ylim([0 4]);
xlim([1 1000]);


fprintf(['\tlength(X_length)              = ' num2str(length(X_length))              '\n']);
fprintf(['\tlength(Y_reads)               = ' num2str(length(Y_reads))               '\n']);
fprintf(['\tlength(Y_reads_ref)           = ' num2str(length(Y_reads_ref))           '\n']);
fprintf(['\tlength(Y_reads_corrected)     = ' num2str(length(Y_reads_corrected))     '\n']);
fprintf(['\tlength(Y_reads_corrected_ref) = ' num2str(length(Y_reads_corrected_ref)) '\n']);
fprintf(['\tlength(X_repet)               = ' num2str(length(X_repet))               '\n']);
fprintf('\n');


fprintf('Subplot 5/6 : [EXPERIMENT] (Corrected copyNum) vs. (Fragment repetitiveness).\n');
if (exist([workingDir 'matlab_dir/' projectName '.Lowess2.mat'],'file') == 0)
    fprintf('\tPreparing for LOWESS fitting.\n');
    X_data   = zeros(numFragments,1);
    Y_data   = zeros(numFragments,1);
    tracking = zeros(numFragments,1);
    for frag = 1:numFragments;
        if (fragment_data(frag).usable == 1)
	    tracking(frag) = 1;
	    X_data(frag)   = X_repet(frag);
	    Y_data(frag)   = Y_reads_corrected(frag);
        end;
    end;
    X_data(tracking == 1) = [];
    Y_data(tracking == 1) = [];
    fprintf('\tLOWESS fitting to project data.\n');
%    [newX2, newY2] = optimize_mylowess(X_data,Y_data,10);
    fprintf('\tLOWESS fitting to project data complete.\n');
%    save([workingDir 'matlab_dir/' projectName '.Lowess2.mat'], 'newX2','newY2');
else
    fprintf('\tLOWESS fitting to project data loading.\n');
%    load([workingDir 'matlab_dir/' projectName '.Lowess2.mat']);
    fprintf('\tLOWESS fitting to project data loaded.\n');
end;
sh(5) = subplot(5,4,17);
hold on;
for frag = 1:numFragments
    X_datum = X_repet(frag);
    Y_datum = Y_reads_corrected(frag);
    if (fragment_data(frag).usable == 1)
        plot(X_datum,Y_datum,'.', 'color', [0.0, 0.66667, 0.33333], 'MarkerSize',0.1);
%    else
%        plot(X_datum,Y_datum,'.', 'color', [1.0, 0.0, 0.0], 'MarkerSize',0.1);
    end;
end;
% plot(newX2,newY2,'color','k','linestyle','-', 'linewidth',1);
hold off;
axis normal;
h = title('Relationship of repetitiveness and corrected reads. [EXPERIMENT]');   set(h, 'FontSize', 5);
h = ylabel('Corrected Reads');           set(h, 'FontSize', 5);
h = xlabel('Fragment Repetitiveness');   set(h, 'FontSize', 5);
set(gca,'FontSize',5);
xlim([0 5000]);
ylim([0 4]);

sh(6) = subplot(5,4,18);
hold on;
for frag = 1:numFragments
    X_datum = X_repet(frag);
    Y_datum = Y_reads_corrected(frag);
    if (fragment_data(frag).usable == 1)
%        plot(X_datum,Y_datum,'.', 'color', [0.0, 0.66667, 0.33333], 'MarkerSize',0.1);
%    else
        plot(X_datum,Y_datum,'.', 'color', [1.0, 0.0, 0.0], 'MarkerSize',0.1);
    end;
end;
hold off;
axis normal;
h = title('Green = usable; Red = unusable.');   set(h, 'FontSize', 5);
set(gca,'FontSize',5);
xlim([0 5000]);
ylim([0 4]);


fprintf('Subplot 6/6 : [EXPERIMENT vs. REFERENCE] (Corrected ave read count) vs. (Corrected ave read count).\n');
% if (exist([workingDir 'matlab_dir/' projectName '.Lowess2.mat'],'file') == 0)
%     fprintf('\tLOWESS fitting to project data.\n');
%     [newX2, newY2] = optimize_mylowess(X_repet,Y_reads_corrected,10);
%     fprintf('\tLOWESS fitting to project data complete.\n');
%     save([workingDir 'matlab_dir/' projectName '.Lowess2.mat'], 'newX2','newY2');
% else
%     fprintf('\tLOWESS fitting to project data loading.\n');
%     load([workingDir 'matlab_dir/' projectName '.Lowess2.mat']);
%     fprintf('\tLOWESS fitting to project data loaded.\n');
% end;
sh(7) = subplot(5,4,19);
plot(Y_reads_corrected,Y_reads_corrected_ref,'.', 'color', [0.0, 0.66667, 0.33333], 'MarkerSize', 1);
% hold on;
% plot(newX2,newY2,'color','k','linestyle','-', 'linewidth',1);
% hold off;
axis normal;
h = title('Examine corrected reads. [EXPERIMENT vs. REFERENCE]');   set(h, 'FontSize', 5);
h = ylabel({'Corrected Reads';'(Reference)'});   set(h, 'FontSize', 5);
h = xlabel({'Corrected Reads';'(Project)'});     set(h, 'FontSize', 5);
set(gca,'FontSize',5);
axis(sh(7),'square');
xlim([0 10]);
ylim([0 10]);

fprintf('Saving figure.\n');
saveas(fig, [figureDir projectName '.examine_bias.eps'], 'epsc');
delete(fig);

%
%%
%%%
%%%%
%%%%
%%%%
%%%
%%
%

%% Save fragment_data for project again, now that we've assigned usable vs. not usable.
fragments_CNV = fragment_data;
save([workingDir 'matlab_dir/' projectName '.fragment_CNV_data.mat'], 'fragments_CNV');

%%=========================================================================
%%= Analyze corrected CNV data for RADseq datasets. =======================
%%=========================================================================

fprintf(['\nGenerating CNV figure from ''' projectName ''' sequence data corrected by removing restriction fragment length bias.\n']);
    
% Initializes vectors used to hold copy number data.
for chr = 1:num_chr   % number of chrs.
    for j = 1:2       % 2 categories tracked : total read counts in bin and number of data entries in region (size of bin in base-pairs).
        chr_CNVdata_RADseq_proj{chr,j} = zeros(1,ceil(chr_size(chr)/bases_per_bin));
	chr_CNVdata_RADseq_ref{chr,j}  = zeros(1,ceil(chr_size(chr)/bases_per_bin));
    end;
end;

% Output chromosome lengths to log file.
for i = 1:length(chr_name)
    fprintf(['\nchr' num2str(i) ' = ''' chr_name{i} '''.\tCGHlength = ' num2str(length(chr_CNVdata_RADseq_proj{i,1}))]);
end;


%-------------------------------------------------------------------------------------------------
% Distribute project data into standard bins.
%-------------------------------------------------------------------------------------------------
chr_CNVdata_RADseq = chr_CNVdata_RADseq_proj;
if (exist([workingDir 'matlab_dir/' projectName '.corrected_CNV_' CNV_verString '.ploidy_' ploidyString '.mat'],'file') == 0)
    fprintf('\nMAT file containing project CNV information not found, regenerating from prior data files.\n');

    % Convert bias corrected ave_copy_number per restriction digest fragment into CNV data for plotting
    fprintf('\n# Adding fragment corrected_ave_read_copy values to map bins.');
    for frag = 1:numFragments
	if (fragment_data(frag).usable == 1)
            % Load important data from fragments data structure.
            chrID      = fragment_data(frag).chr;
            posStart   = fragment_data(frag).startbp;
            posEnd     = fragment_data(frag).endbp;
	    count      = fragment_data(frag).reads_ave_corrected;
	    fragLength = fragment_data(frag).length;

            % Identify locations of fragment end in relation to plotting bins.
            val1 = ceil(posStart/bases_per_bin);
            val2 = ceil(posEnd  /bases_per_bin);

	    if (chrID > 0) && (count > 0)
		% 'count' is average read count across fragment, so it will need multiplied by fragment length (or fraction) before
		%     adding to each bin.
		if (val1 == val2)
		    % All of the restriction fragment belongs to one bin.
		    if (val1 <= length(chr_CNVdata_RADseq{chrID,1}))
			chr_CNVdata_RADseq{chrID,1}(val1) = chr_CNVdata_RADseq{chrID,1}(val1)+count*fragLength;
			chr_CNVdata_RADseq{chrID,2}(val1) = chr_CNVdata_RADseq{chrID,2}(val1)+fragLength;
		    end;
		else % (val1 < val2)
		    % The restriction fragment belongs partially to two bins, so we must determine fraction assigned to each bin.
		    posEdge     = val1*bases_per_bin;
		    fragLength1 = posEdge-posStart+1;
		    fragLength2 = posEnd-posEdge;

		    % Add data to first bin.
		    if (val1 <= length(chr_CNVdata_RADseq{chr,1}))
			chr_CNVdata_RADseq{chrID,1}(val1) = chr_CNVdata_RADseq{chrID,1}(val1)+count*fragLength1;
			chr_CNVdata_RADseq{chrID,2}(val1) = chr_CNVdata_RADseq{chrID,2}(val1)+fragLength1;
		    end;

		    % Add data to second bin.
		    if (val2 <= length(chr_CNVdata_RADseq{chr,1}))
			chr_CNVdata_RADseq{chrID,1}(val2) = chr_CNVdata_RADseq{chrID,1}(val2)+count*fragLength2;
			chr_CNVdata_RADseq{chrID,2}(val2) = chr_CNVdata_RADseq{chrID,2}(val2)+fragLength2;
		    end;
                end;
            end;
	end;
    end;
    fprintf('\n# Fragment corrected_ave_read_copy values have been added to map bins.');

    save([workingDir 'matlab_dir/' projectName '.corrected_CNV_' CNV_verString '.ploidy_' ploidyString '.mat'],...
	 'chr_CNVdata_RADseq');
else
    fprintf('\nProject CNV MAT file found, loading.\n');
    load([workingDir 'matlab_dir/' projectName '.corrected_CNV_' CNV_verString '.ploidy_' ploidyString '.mat']);
end;
chr_CNVdata_RADseq_proj = chr_CNVdata_RADseq;


%-------------------------------------------------------------------------------------------------
% Distribute parent data into standard bins.
%-------------------------------------------------------------------------------------------------
chr_CNVdata_RADseq = chr_CNVdata_RADseq_ref;
if (exist([workingDir 'matlab_dir/' referenceName '.corrected_CNV_' CNV_verString '.ploidy_' num2str(ploidy_default) '.mat'],'file') == 0)
    fprintf('\nMAT file containing reference CNV information not found, regenerating from prior data files.\n');
 
    % Convert bias corrected ave_copy_number per restriction digest fragment into CNV data for plotting
    fprintf('\n# Adding fragment corrected_ave_read_copy values to map bins.');
    for frag = 1:numFragments
        if (fragment_data(frag).usable == 1)
	    % Load important data from fragments data structure.
	    chrID      = fragment_data(frag).chr;    
	    posStart   = fragment_data(frag).startbp;
	    posEnd     = fragment_data(frag).endbp; 
	    count      = fragment_data(frag).reads_ave_corrected_ref;
	    fragLength = fragment_data(frag).length;

            % Identify locations of fragment end in relation to plotting bins.
            val1 = ceil(posStart/bases_per_bin);
            val2 = ceil(posEnd  /bases_per_bin);

            if (chrID > 0) && (count > 0)   
                % 'count' is average read count across fragment, so it will need multiplied by fragment length (or fraction) before
                %     adding to each bin.
                if (val1 == val2)
                    % All of the restriction fragment belongs to one bin.
                    if (val1 <= length(chr_CNVdata_RADseq{chrID,1}))
                        chr_CNVdata_RADseq{chrID,1}(val1) = chr_CNVdata_RADseq{chrID,1}(val1)+count*fragLength;
                        chr_CNVdata_RADseq{chrID,2}(val1) = chr_CNVdata_RADseq{chrID,2}(val1)+fragLength;
                    end;
                else % (val1 < val2)
                    % The restriction fragment belongs partially to two bins, so we must determine fraction assigned to each bin.
                    posEdge     = val1*bases_per_bin;
                    fragLength1 = posEdge-posStart+1;
                    fragLength2 = posEnd-posEdge;
            
                    % Add data to first bin.
                    if (val1 <= length(chr_CNVdata_RADseq{chr,1}))
                        chr_CNVdata_RADseq{chrID,1}(val1) = chr_CNVdata_RADseq{chrID,1}(val1)+count*fragLength1;
                        chr_CNVdata_RADseq{chrID,2}(val1) = chr_CNVdata_RADseq{chrID,2}(val1)+fragLength1;
                    end;
            
                    % Add data to second bin.
                    if (val2 <= length(chr_CNVdata_RADseq{chr,1}))
                        chr_CNVdata_RADseq{chrID,1}(val2) = chr_CNVdata_RADseq{chrID,1}(val2)+count*fragLength2;
                        chr_CNVdata_RADseq{chrID,2}(val2) = chr_CNVdata_RADseq{chrID,2}(val2)+fragLength2;
                    end;
                end;
            end;
	end;
    end;
    fprintf('\n# Fragment corrected_ave_read_copy values have been added to map bins.');

    save([workingDir 'matlab_dir/' referenceName '.corrected_CNV_' CNV_verString '.ploidy_' num2str(ploidy_default) '.mat'],...
         'chr_CNVdata_RADseq');
else
    fprintf('\nRefrence CNV MAT file found, loading.\n');
    load([workingDir 'matlab_dir/' referenceName '.corrected_CNV_' CNV_verString '.ploidy_' num2str(ploidy_default) '.mat']);
end;
chr_CNVdata_RADseq_ref = chr_CNVdata_RADseq;

        
% basic plot parameters not defined per genome.
TickSize        = -0.005;  %negative for outside, percentage of longest chr figure.
maxY            = ploidyBase*2;

%% -----------------------------------------------------------------------------------------
% Make figures
%-------------------------------------------------------------------------------------------
fig = figure(1);
set(gcf, 'Position', [0 70 1024 600]);
for chr = 1:num_chr
    for pos = 1:length(chr_CNVdata_RADseq{chr,1})
            
%	% Plot the sum of the data in each region, divided by the number of data points in each region.
%	if (chr_CNVdata_RADseq{chr,2}(pos) == 0)
%	    % No data elements => null value is plotted.
%	    CNVplot{chr}(pos) = 0;
%	else
%	    % Sum of data elements is divided by the number of data elements.
%	    CNVplot{chr}(pos) = chr_CNVdata_RADseq_proj{chr,1}(pos)/chr_CNVdata_RADseq_proj{chr,2}(pos);
%	end;

        % Plot the sum of the data in each region, divided by the number of data points in each region.
        % Then divided by this value calculated for SC5314 data.
        if (chr_CNVdata_RADseq_proj{chr,2}(pos) == 0) || (chr_CNVdata_RADseq_ref{chr,2}(pos) == 0)
            % No data elements => null value is plotted.
            CNVplot{chr}(pos) = 0;
        else
            % Sum of data elements is divided by the number of data elements.
            normalizationFactor = chr_CNVdata_RADseq_ref{chr,1}(pos)/chr_CNVdata_RADseq_ref{chr,2}(pos);
            if (normalizationFactor == 0)
                CNVplot{chr}(pos) = 0;
            else
                CNVplot{chr}(pos) = chr_CNVdata_RADseq_proj{chr,1}(pos)/chr_CNVdata_RADseq_proj{chr,2}(pos)/normalizationFactor;
            end;
        end;
    end;
    chr_max(chr) = max(CNVplot{chr});
    chr_med(chr) = median(CNVplot{chr});
end;
max_count     = max(chr_max);
median_count  = sum(chr_med)/length(chr_med);
for chr = 1:num_chr
    CNVplot2{chr} = CNVplot{chr}/median_count;
end;

% Save presented CNV data in a file format common across data types being processed.
fprintf('\nSaving "common_CNV" data file.');
genome_CNV = genome_data;
save([workingDir 'matlab_dir/' projectName '.ploidy_' ploidyString '.common_CNV.mat'], 'CNVplot2','genome_CNV');

% totalCNVdata = [];
% for chr = 1:num_chr
%     totalCNVdata = [totalCNVdata CNVplot{chr}];
% end;
% max_count    = max(totalCNVdata);
% min_count    = min(totalCNVdata);
% mean_count   = mean(totalCNVdata);
% median_count = median(totalCNVdata);
% for chr = 1:num_chr
%     CNVplot2{chr} = CNVplot{chr}/median_count;
% end;

ploidy = str2num(ploidyString);
fprintf(['Ploidy string = "' ploidyString '"\n']);
[chr_breaks, chrCopyNum, ploidyAdjust] = FindChrSizes_2(Aneuploidy,CNVplot2,ploidy,num_chr);
largestChr = find(chr_width == max(chr_width));

%% -----------------------------------------------------------------------------------------
% Output corrected CNV data to pileup-type file.
%-------------------------------------------------------------------------------------------
%% Output corrected reads, multiplied by ploidyEstimate, to : [workingDir 'pileup_dir/' projectName '_RADseq_digest_analysis_CNV_corrected.txt']
%  Including : [chrNum, bpStart, bpEnd, correctedReads*ploidy, fragmentLength]
%   1       9638    10115   4796    478
%   1       10116   10123   1283    8
%   1       13170   13841   55175   672
fprintf('\nOutputting fragment read copy after correction.\n');
outfile     = [workingDir 'pileup_dir/' projectName '_RADseq_digest_analysis_CNV_corrected.txt'];
outfile_fid = fopen(outfile,'w');
for frag = 1:numFragments
    string1 = num2str(fragment_data(frag).chr);
    string2 = num2str(fragment_data(frag).startbp);
    string3 = num2str(fragment_data(frag).endbp);
    if (Low_quality_ploidy_estimate == true)
	string4 = num2str(fragment_data(frag).reads_ave_corrected/fragment_data(frag).reads_ave_corrected_ref*ploidy*ploidyAdjust);
    else
	string4 = num2str(fragment_data(frag).reads_ave_corrected/fragment_data(frag).reads_ave_corrected_ref*ploidy);
    end;
    string5 = num2str(fragment_data(frag).length);
    if (fragment_data(frag).usable == 1)
	fprintf(outfile_fid,[string1 '\t' string2 '\t' string3 '\t' string4 '\t' string5 '\t[*]\n']);
    else
	fprintf(outfile_fid,[string1 '\t' string2 '\t' string3 '\t' string4 '\t' string5 '\n']);
    end;
end;
fclose(outfile_fid);

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
% initialize string to contain chromosome copy number estimates for later output.
stringChrCNVs = '';

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
        title([ projectName ' CNV map'],'Interpreter','none','FontSize',12);
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
        
        % draw outlines of chromosome cartoon.   (drawn after horizontal lines to that cartoon edges are not interrupted by horiz lines.
        plot([leftEnd   leftEnd   leftEnd+dx   x1-dx   x1        x2        x2+dx   rightEnd-dx   rightEnd   rightEnd   rightEnd-dx x2+dx   x2   x1   x1-dx   leftEnd+dx   leftEnd],...
             [dy        maxY-dy   maxY         maxY    maxY-dy   maxY-dy   maxY    maxY          maxY-dy    dy         0           0       dy   dy   0       0            dy     ],...
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

	    % DDD
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

	%% draw lines across plots for easier interpretation of CNV regions.   DDD
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

            % draw outlines of chromosome cartoon.   (drawn after horizontal lines to that cartoon edges are not interrupted by
            % horiz lines.
            plot([leftEnd   leftEnd   leftEnd+dx   x1-dx   x1        x2        x2+dx   rightEnd-dx   rightEnd   rightEnd ...
                  rightEnd-dx   x2+dx   x2   x1   x1-dx   leftEnd+dx  leftEnd],...
                 [dy        maxY-dy   maxY         maxY    maxY-dy   maxY-dy   maxY    maxY          maxY-dy    dy       ...
                  0             0       dy   dy   0       0           dy  ],...
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
            ylabel(projectName, 'Rotation', 0, 'HorizontalAlign', 'right', 'VerticalAlign', 'bottom','Interpreter','none','FontSize',5);

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
            
%       % Main figure colors key.
%       left   = key_posX;
%       bottom = key_posY;
%       width  = str2num(key_width);
%       height = key_height;
%       colorNoData = [1.0   1.0   1.0  ]; %used when no data is available for the bin.
%       colorCNV    = [0.0   0.0   0.0  ]; %external; used in blending at ends of chr.
            
%       subplot('Position',[left bottom width height]); 
%       axis off square;
%       xlim([-0.1,1]);
%       ylim([-0.1,1.6]);
%       set(gca,'XTick',[]);
%       set(gca,'YTick',[]);
%       patch([0 0.2 0.2 0], [1.4 1.4 1.5 1.5], colorCNV);      text(0.3,1.05,'Copy number variation (CNV).');
%       patch([0 0.2 0.2 0], [1.0 1.0 1.1 1.1], colorNoData);   text(0.3,1.45,'No CNV.');

% Make figure output have transparent background.
set(gcf, 'color', 'none',...
         'inverthardcopy', 'off');
            
% Save original arrangement of chromosomes.
saveas(fig, [figureDir projectName '.corrected-CNV-map.1.eps'], 'epsc');
delete(fig);
    
% Save horizontally arranged chromosomes.
set(Linear_fig,'PaperPosition',[0 0 8 0.62222222]);
saveas(Linear_fig, [figureDir projectName '.corrected-CNV-map.2.eps'], 'epsc');
delete(Linear_fig);

% Output chromosome copy number estimates.
textFileName = [figureDir projectName '.CNV-map.3.txt'];
fprintf(['Text output of CNVs : "' textFileName '"\n']);
textFileID = fopen(textFileName,'w');
fprintf(textFileID,stringChrCNVs);
fclose(textFileID);

end
