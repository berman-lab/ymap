function [] = CNV_v6_fragmentLengthCorrected_5(main_dir,user,genome,genomeUser,project,parentProject,parentUser,ploidyEstimate,ploidyBase, ...
                                 CNV_verString,displayBREAKS);
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
% Load FASTA file name from 'reference.txt' file for project.
%--------------------------------------------------------------------------
Reference    = [main_dir 'users/' genomeUser '/genomes/' genome '/reference.txt'];
FASTA_string = strtrim(fileread(Reference));
[FastaPath,FastaName,FastaExt] = fileparts(FASTA_string);


%%=========================================================================
% Control variables.
%--------------------------------------------------------------------------
projectDir = [main_dir 'users/' user       '/projects/' project '/'];
parentDir  = [main_dir 'users/' parentUser '/projects/' parentProject '/'];
genomeDir  = [main_dir 'users/' genomeUser '/genomes/'  genome '/'];

fprintf(['\n### projectDir : ' projectDir '\n']);
fprintf([  '### parentDir  : ' parentDir  '\n']);
fprintf([  '### genomeDir  : ' genomeDir  '\n']);
fprintf([  '### genome     : ' genome     '\n']);
fprintf([  '### project    : ' project    '\n']);
[centromeres, chr_sizes, figure_details, annotations, ploidy_default] = Load_genome_information_1(main_dir,genomeDir,genome);
[Aneuploidy] = [];

if (strcmp(project,parentProject) == 1)
	ParentValid = 0;
else
	ParentValid = 1;
end;

for i = 1:length(chr_sizes)
	chr_size(chr_sizes(i).chr)    = chr_sizes(i).size;
end;
for i = 1:length(centromeres)
	cen_start(centromeres(i).chr) = centromeres(i).start;
	cen_end(centromeres(i).chr)   = centromeres(i).end;
end;
if (length(annotations) > 0)
	for i = 1:length(annotations)
		annotation_chr(i)       = annotations(i).chr;
		annotation_type{i}      = annotations(i).type;
		annotation_start(i)     = annotations(i).start;
		annotation_end(i)       = annotations(i).end;
		annotation_fillcolor{i} = annotations(i).fillcolor;
		annotation_edgecolor{i} = annotations(i).edgecolor;
		annotation_size(i)      = annotations(i).size;
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
num_chrs      = length(chr_size);


%%=========================================================================
%%= No further control variables below. ===================================
%%=========================================================================

% Sanitize user input of euploid state.
ploidyBase = round(str2num(ploidyBase));
if (ploidyBase > 4);   ploidyBase = 4;   end;
if (ploidyBase < 1);   ploidyBase = 1;   end;
fprintf(['\nEuploid base = "' num2str(ploidyBase) '"\n']);

% basic plot parameters not defined per genome.
TickSize         = -0.005;  %negative for outside, percentage of longest chr figure.
bases_per_bin    = max(chr_size)/700;
maxY             = ploidyBase*2;
cen_tel_Xindent  = 5;
cen_tel_Yindent  = maxY/5;

%% Determine reference genome FASTA file in use.
%  Read in and parse : "links_dir/main_script_dir/genome_specific/[genome]/reference.txt"
reference_file   = [main_dir 'users/' genomeUser '/genomes/' genome '/reference.txt'];
refernce_fid     = fopen(reference_file, 'r');
refFASTA         = fgetl(refernce_fid);
fclose(refernce_fid);


%%================================================================================================
% Load pre-processed ddRADseq fragment CNV data for project.
%-------------------------------------------------------------------------------------------------
if (exist([main_dir 'users/' user '/projects/' project '/fragment_CNV_data.mat'],'file') == 0)
    %  Including : [chrNum, bpStart, bpEnd, maxReads, AveReads, fragmentLength]
    %	1	9638	10115	2	0	478
    %	1	10116	10123	2	1	8
    %	1	13170	13841	0	0	672
    fprintf('Loading results from Python script, which pre-processed the dataset relative to genome restriction fragments : project\n');
    datafile_RADseq  = [main_dir 'users/' user '/projects/' project '/preprocessed_CNVs.ddRADseq.txt'];
    data_RADseq      = fopen(datafile_RADseq);
    count            = 0;
    fragments_CNV    = [];
    while ~feof(data_RADseq)
		% Load fragment data from pre-processed text file, single line.
		dataLine = fgetl(data_RADseq);
		if (length(dataLine) > 0)
			if (dataLine(1) ~= '#')   % If the first space-delimited string is not '#', then process it as a valid data line.
			    % The number of valid lines found so far...  the number of usable restriction fragments with data so far.
			    count = count + 1;

			    % Each valid data line consists of four tab-delimited collumns : [chrNum, bpStart, bpEnd, Ave_reads]
				chr_num         = str2num(sscanf(dataLine, '%s',1));
				bp_start        = sscanf(dataLine, '%s',  2 );    for i = 1:size(sscanf(dataLine,'%s', 1 ),2);   bp_start(1)  = [];   end;   bp_start  = str2num(bp_start);
				bp_end          = sscanf(dataLine, '%s',  3 );    for i = 1:size(sscanf(dataLine,'%s', 2 ),2);   bp_end(1)    = [];   end;   bp_end    = str2num(bp_end);
				ave_reads       = sscanf(dataLine, '%s',  4 );    for i = 1:size(sscanf(dataLine,'%s', 3 ),2);   ave_reads(1) = [];   end;   ave_reads = str2num(ave_reads);
				fragment_length = bp_end - bp_start + 1;

			    % Add fragment data to data structure.
			    fragments_CNV(count).chr        = chr_num;
			    fragments_CNV(count).startbp    = bp_start;
			    fragments_CNV(count).endbp      = bp_end;
			    fragments_CNV(count).length     = fragment_length;
			    fragments_CNV(count).ave_reads  = ave_reads;
			    fragments_CNV(count).usable     = 1;
			end;
		end;
	end;
	fclose(data_RADseq);
	save([main_dir 'users/' user '/projects/' project '/fragment_CNV_data.mat'], 'fragments_CNV');
else
	load([main_dir 'users/' user '/projects/' project '/fragment_CNV_data.mat']);
end;
project_fragments_CNV = fragments_CNV;
clear fragments_CNV;


%%================================================================================================
% Load pre-processed ddRADseq fragment CNV data for parent project.
%-------------------------------------------------------------------------------------------------
if (strcmp(project,parentProject) == 1)
	parent_fragments_CNV = project_fragments_CNV;
else
	load([main_dir 'users/' parentUser '/projects/' parentProject '/fragment_CNV_data.mat']);
	parent_fragments_CNV = fragments_CNV;
	clear fragments_CNV;
end;


%%================================================================================================
% Load pre-processed ddRADseq fragment GC-bias data for genome.
%-------------------------------------------------------------------------------------------------
fprintf(['standard_bins_GC_ratios_file :\n\t' main_dir 'users/' genomeUser '/genomes/' genome '/' FastaName '.GC_ratios.MfeI_MboI.txt\n']);
GC_ratios_fid = fopen([main_dir 'users/' genomeUser '/genomes/' genome '/' FastaName '.GC_ratios.MfeI_MboI.txt'], 'r');
fprintf(['\t' num2str(GC_ratios_fid) '\n']);
fragID = 0;
while not (feof(GC_ratios_fid))
	dataLine = fgetl(GC_ratios_fid);
	if (length(dataLine) > 0)
		if (dataLine(1) ~= '#')
			% The number of valid lines found so far...  the number of usable restriction fragments with data so far.
			fragID              = fragID + 1;
			chr                 = str2num(sscanf(dataLine, '%s',1));
			fragment_start      = sscanf(dataLine, '%s',2);  for i = 1:size(sscanf(dataLine,'%s',1),2);      fragment_start(1) = []; end;    fragment_start = str2num(fragment_start);
			fragment_end        = sscanf(dataLine, '%s',3);  for i = 1:size(sscanf(dataLine,'%s',2),2);      fragment_end(1)   = []; end;    fragment_end   = str2num(fragment_end);
			GCratio             = sscanf(dataLine, '%s',4);  for i = 1:size(sscanf(dataLine,'%s',3),2);      GCratio(1)        = []; end;    GCratio        = str2num(GCratio);
			GCratioData(fragID) = GCratio;
		end;
	end;
end;
fclose(GC_ratios_fid);


%%================================================================================================
% Add reference and GC-bias data to common data structure : 'fragment_data'
%-------------------------------------------------------------------------------------------------
fragment_data = project_fragments_CNV;
numFragments  = length(fragment_data);
for fragID = 1:numFragments
	% Add GCratio data to the structure.
	fragment_data(fragID).GC_bias  = GCratioData(fragID);

	% Add parent project read depth to structure.
	fragment_data(fragID).ave_reads_parent = parent_fragments_CNV(fragID).ave_reads;

	% Initialize all data as being usable.
	fragment_data(fragID).usable = 1;

%	if (fragment_data(fragID).ave_reads <= 1)
%		fragment_data(fragID).usable = 0;
%	end;
%	if (ParentValid == 1)
%		if (fragment_data(fragID).ave_reads_parent <= 1)
%	        fragment_data(fragID).usable = 0;
%	    end;
%	end;
end;


%%================================================================================================
% Determine which fragments have usable data.
%-------------------------------------------------------------------------------------------------
fprintf('Gathering [fragment length], [CNV], & [GC bias] data for plotting bias and correction.\n');
count1 = 0;
count2 = 0;
count3 = 0;
count4 = 0;
for fragID = 1:numFragments
	if ((fragment_data(fragID).length == 0) || (fragment_data(fragID).ave_reads == 0))
		fragment_data(fragID).usable = 0;
		count1 = count1 + 1;
	else
		fragment_data(fragID).usable = 1;
		count2 = count2 + 1;
		count3 = count3 + 1;
		if (fragment_data(fragID).ave_reads_parent == 0)
			%fragment_data(fragID).usable = 0;
			count2 = count2 - 1;
			count1 = count1 + 1;
		else
			count4 = count4 + 1;
		end;
	end;
end;
fprintf(['unusable loci       = ' num2str(count1) '\n']);
fprintf(['usable loci         = ' num2str(count2) '\n']);
fprintf(['usable project loci = ' num2str(count3) '\n']);
fprintf(['usable parent loci  = ' num2str(count4) '\n']);


%%================================================================================================
% Prepare for LOWESS fitting processes : remove empty entries.
%-------------------------------------------------------------------------------------------------
X_length        = zeros(1,numFragments);
X_GCbias        = zeros(1,numFragments);
Y_reads_project = zeros(1,numFragments);
Y_reads_parent  = zeros(1,numFragments);
for fragID = 1:numFragments
    X_length(fragID)        = fragment_data(fragID).length; 
    X_GCbias(fragID)        = fragment_data(fragID).GC_bias;
    Y_reads_project(fragID) = fragment_data(fragID).ave_reads;
    Y_reads_parent(fragID)  = fragment_data(fragID).ave_reads_parent;
end;
fprintf(['\tnumFragments            = ' num2str(numFragments)            '\n']);
fprintf(['\tlength(fragment_data)   = ' num2str(length(fragment_data))   '\n']);
fprintf(['\tlength(X_length)        = ' num2str(length(X_length))        '\n']);
fprintf(['\tlength(Y_reads_project) = ' num2str(length(Y_reads_project)) '\n']);
fprintf(['\tlength(Y_reads_parent)  = ' num2str(length(Y_reads_parent))  '\n']);
fprintf(['\tlength(X_GCbias)        = ' num2str(length(X_GCbias))        '\n']);


%============================================================*
% Trim main dataset before initial LOWESS fitting.           |
%------------------------------------------------------------*
	% Remove fragments longer than 1000bp.
	Y_reads_project( X_length > 1000) = [];
	Y_reads_parent(  X_length > 1000) = [];
	X_GCbias(        X_length > 1000) = [];
	X_length(        X_length > 1000) = [];
%------------------------------------------------------------*
	% Remove fragments with less than 1.5 average reads.
	Y_reads_parent(  Y_reads_project < 1.5) = [];
	X_GCbias(        Y_reads_project < 1.5) = [];
	X_length(        Y_reads_project < 1.5) = [];
	Y_reads_project( Y_reads_project < 1.5) = [];
%------------------------------------------------------------*
%                                                            |
%============================================================*

fprintf('\n');
fprintf(['\tlength(X_length)        = ' num2str(length(X_length))        '\n']);
fprintf(['\tlength(Y_reads_project) = ' num2str(length(Y_reads_project)) '\n']);
fprintf(['\tlength(Y_reads_parent)  = ' num2str(length(Y_reads_parent))  '\n']);
fprintf(['\tlength(X_GCbias)        = ' num2str(length(X_GCbias))        '\n']);
fprintf('\n');

% subplot handle list.
sh=zeros(11,1);


%%================================================================================================
% Construct figures illustrating LOWESS fitting and normalization.
%-------------------------------------------------------------------------------------------------
fprintf('Generating figure with plot of length vs. coverage data.\n');
fig = figure(1);


%-------------------------------------------------------------------------------------------------
% Plotting average read depth vs. fragment length for project.
%-------------------------------------------------------------------------------------------------
fprintf('Subplot 1/11 : [EXPERIMENT] (Ave read depth) vs. (Fragment length).\n');
sh(1) = subplot(5,4,[1 2]);
fprintf('\tLOWESS fitting to project data.\n');
[newX1_project, newY1_project] = optimize_mylowess(X_length,Y_reads_project,10);
fprintf('\tLOWESS fitting to project data complete.\n');
plot(X_length,Y_reads_project,'.', 'color', [0.0, 0.66667, 0.33333], 'MarkerSize', 1);
hold on;
plot(newX1_project, newY1_project,'color','k','linestyle','-', 'linewidth',1);
hold off;
axis normal;
h = title('Examine bias between ddRADseq fragment length and read count. [EXPERIMENT]');   set(h, 'FontSize', 5);
h = ylabel('Average Reads');     set(h, 'FontSize', 5);
h = xlabel('Fragment Length');   set(h, 'FontSize', 5);
set(gca,'FontSize',5);
ylim([1 4*max(newY1_project)]);
xlim([0 1000]);

%-------------------------------------------------------------------------------------------------
% Plotting average read depth vs. fragment length for parent.
%-------------------------------------------------------------------------------------------------
fprintf('Subplot 2/11 : [REFERENCE] (Ave read count) vs. (Fragment length).\n');
sh(2) = subplot(5,4,[3 4]);
fprintf('\tLOWESS fitting to reference data.\n');
[newX1_parent, newY1_parent] = optimize_mylowess(X_length,Y_reads_parent,10);
fprintf('\tLOWESS fitting to referemce data complete.\n');
if (ParentValid == 1)
	plot(X_length,Y_reads_parent,'.', 'color', [0.0, 0.66667, 0.33333], 'MarkerSize', 1);
	hold on;
	plot(newX1_parent, newY1_parent,'color','k','linestyle','-', 'linewidth',1);
	hold off;
	axis normal;
	h = title('Examine bias between ddRADseq fragment length and read count. [REFERENCE]');   set(h, 'FontSize', 5);
	h = ylabel('Average Reads');     set(h, 'FontSize', 5);
	h = xlabel('Fragment Length');   set(h, 'FontSize', 5);
	set(gca,'FontSize',5);
	ylim([1 4*max(newY1_parent)]);
	xlim([0 1000]);
end;

%-------------------------------------------------------------------------------------------------
% Resetting data vectors
%-------------------------------------------------------------------------------------------------
X_length        = zeros(1,numFragments);
X_GCbias        = zeros(1,numFragments);
Y_reads_project = zeros(1,numFragments);
Y_reads_parent  = zeros(1,numFragments);
for frag = 1:numFragments
	X_length(frag)        = fragment_data(frag).length;
	X_GCbias(frag)        = fragment_data(frag).GC_bias;
	Y_reads_project(frag) = fragment_data(frag).ave_reads;
	Y_reads_parent(frag)  = fragment_data(frag).ave_reads_parent;
end;

%-------------------------------------------------------------------------------------------------
% Plotting corrected read depth vs. fragment length for project.
%-------------------------------------------------------------------------------------------------
fprintf('Subplot 3/11 : [EXPERIMENT] (Ave read count, with bias corrected) vs. (Fragment length).\n');
sh(3) = subplot(5,4,[5 6]);
hold on;
Y_target = 1;
% Calculate bias_corrected ave_read_copy data for plotting and later analysis.
Y_fitCurve1         = interp1(newX1_project,newY1_project,X_length,'spline');  
Y_reads_corrected_1 = Y_reads_project./Y_fitCurve1*Y_target;
fprintf(['\tlength(Y_fitCurve1)         = ' num2str(length(Y_fitCurve1))         '\n']);
fprintf(['\tlength(Y_reads_corrected_1) = ' num2str(length(Y_reads_corrected_1)) '\n']);
fprintf('\t');
for fragID = 1:numFragments
	% Add corrected ave_data to fragment data structure.
	fragment_data(fragID).reads_ave_corrected_1 = Y_reads_corrected_1(fragID);

	% Define data as no useful if correction fit term falls below 5 reads.
	if (Y_fitCurve1(fragID)         <= 1   );   fragment_data(fragID).usable = 0;   end;
	if (Y_reads_project(fragID)     <= 1   );   fragment_data(fragID).usable = 0;   end;
	if (Y_reads_corrected_1(fragID) <= 0   );   fragment_data(fragID).usable = 0;   end;
	if (X_length(fragID)            <  50  );   fragment_data(fragID).usable = 0;   end;
	if (X_length(fragID)            >  1000);   fragment_data(fragID).usable = 0;   end;

	X_datum = fragment_data(fragID).length;                  % X_Length
	Y_datum = fragment_data(fragID).reads_ave_corrected_1;   % Y_reads_corrected_1
    if (mod(fragID,1000) == 0)
        fprintf('.');
    end;
end;
for fragID = 1:numFragments
	% Plot each corrected data point, colored depending on if data is useful or not.
	if (fragment_data(fragID).usable == 1)
		plot(X_datum,Y_datum,'.', 'color', [0.0, 0.66667, 0.33333], 'MarkerSize',0.2);
	else
		plot(X_datum,Y_datum,'.', 'color', [1.0, 0.0, 0.0], 'MarkerSize',0.2);
	end;
end;
fprintf('\n');
plot(newX1_project,Y_target,'color','k','linestyle','-', 'linewidth',1);
hold off;
axis normal;
h = title('Corrected Read count vs. ddRADseq fragment length. [EXPERIMENT]');   set(h, 'FontSize', 5);
h = ylabel('Corrected Reads');   set(h, 'FontSize', 5);
h = xlabel('Fragment Length');   set(h, 'FontSize', 5);
set(gca,'FontSize',5);
ylim([0 4]);
xlim([1 1000]);

%-------------------------------------------------------------------------------------------------
% Plotting corrected read depth vs. fragment length for parent.
%-------------------------------------------------------------------------------------------------
fprintf('Subplot 4/11 : [REFERENCE] (Ave read count, with bias corrected) vs. (Fragment length).\n');
sh(4) = subplot(5,4,[7 8]);
hold on;
Y_target = 1;
% Calculate bias_corrected ave_read_copy data for plotting and later analysis.
Y_fitCurve1_parent         = interp1(newX1_parent,newY1_parent,X_length,'spline');   
Y_reads_corrected_1_parent = Y_reads_parent./Y_fitCurve1_parent*Y_target;
fprintf(['\tlength(Y_fitCurve1_parent)         = ' num2str(length(Y_fitCurve1_parent))         '\n']);
fprintf(['\tlength(Y_reads_corrected_1_parent) = ' num2str(length(Y_reads_corrected_1_parent)) '\n']);
fprintf('\t');
for fragID = 1:numFragments
	% Add corrected ave_data to fragment data structure.
	fragment_data(fragID).reads_ave_corrected_parent_1 = Y_reads_corrected_1_parent(fragID);

	% Define data as not useful if correction fit term falls below 5 reads.
	fragment_data(fragID).usable_parent = 1;
	if (Y_fitCurve1_parent(fragID)         <= 1   );   fragment_data(fragID).usable_parent = 0;   end;
	if (Y_reads_parent(fragID)             <= 1   );   fragment_data(fragID).usable_parent = 0;   end;
	if (Y_reads_corrected_1_parent(fragID) <= 0   );   fragment_data(fragID).usable_parent = 0;   end;
	if (X_length(fragID)                   <  50  );   fragment_data(fragID).usable_parent = 0;   end;
	if (X_length(fragID)                   >  1000);   fragment_data(fragID).usable_parent = 0;   end;

	X_datum = fragment_data(fragID).length;                         % X_Length
	Y_datum = fragment_data(fragID).reads_ave_corrected_parent_1;   % Y_reads_corrected_1_parent

	if (mod(fragID,1000) == 0)
		fprintf('.');
	end;
end;
fprintf('\n');
if (ParentValid == 1)
	for fragID = 1:numFragments
		% Plot each corrected data point, colored depending on if data is useful or not.
		if (fragment_data(fragID).usable_parent == 1)
			plot(X_datum,Y_datum,'.', 'color', [0.0, 0.66667, 0.33333], 'MarkerSize',0.2);
		else
			plot(X_datum,Y_datum,'.', 'color', [1.0, 0.0, 0.0], 'MarkerSize',0.2);
		end;
	end;
	plot(newX1_parent,Y_target,'color','k','linestyle','-', 'linewidth',1);
	hold off;
	axis normal;
	h = title('Corrected Read count vs. ddRADseq fragment length. [REFERENCE]');   set(h, 'FontSize', 5);
	h = ylabel('Corrected Reads');   set(h, 'FontSize', 5);
	h = xlabel('Fragment Length');   set(h, 'FontSize', 5);
	set(gca,'FontSize',5);
	ylim([0 4]);
	xlim([1 1000]);
end;


%-------------------------------------------------------------------------------------------------
% Gather usable sub-selections of data for LOWESS fitting of GC_bias vs read count (previously corrected for length bias), for project.
%-------------------------------------------------------------------------------------------------
fprintf('\tPreparing for LOWESS fitting : GC_ratio vs corrected_count_1.\n');
X_data   = zeros(1,numFragments);
Y_data   = zeros(1,numFragments);
tracking = zeros(1,numFragments);

fragIDs  = 1:numFragments;
X_data   = X_GCbias;
Y_data   = Y_reads_corrected_1;
for fragID = 1:numFragments;
	if (fragment_data(fragID).usable == 1)
		tracking(fragID) = 1;
	end;
end;
X_data(tracking == 0) = [];
Y_data(tracking == 0) = [];


%% Perform LOWESS fitting.
fprintf('\tLOWESS fitting to project data.\n');
[newX2_project, newY2_project] = optimize_mylowess(X_data,Y_data,10);
fprintf('\tLOWESS fitting to project data complete.\n');
Y_target = 1;
% Calculate repetitiveness_bias_corrected length_bia_corrected ave_read_count data for plotting and later analysis.
Y_fitCurve2         = interp1(newX2_project,newY2_project,X_GCbias,'spline');
fprintf(['& length(Y_reads_corrected_1) = ' num2str(length(Y_reads_corrected_1)) '\n']);    % length(Y_reads_corrected_1) = 17676
fprintf(['& length(Y_fitCurve2)         = ' num2str(length(Y_fitCurve2        )) '\n']);    % length(Y_fitCurve2)         = 10401
Y_reads_corrected_2 = Y_reads_corrected_1./Y_fitCurve2*Y_target;


%-------------------------------------------------------------------------------------------------
% Gather usable sub-selections of data for LOWESS fitting of GC_ratio vs read count (previously corrected for length bias), for parent.
%-------------------------------------------------------------------------------------------------
fprintf('\tPreparing for LOWESS fitting : GC_ratio vs corrected_count_1.\n');
tracking = zeros(1,numFragments);
fragIDs  = 1:numFragments;
X_data   = X_GCbias;
Y_data   = Y_reads_corrected_1_parent;
for fragID = 1:numFragments;
	if (fragment_data(fragID).usable_parent == 1)
		tracking(fragID) = 1;
	end;
end;
X_data(tracking == 0)  = [];
Y_data(tracking == 0)  = [];
%% Perform LOWESS fitting.
fprintf('\tLOWESS fitting to reference data.\n');
[newX2_parent, newY2_parent] = optimize_mylowess(X_data,Y_data,10);
fprintf('\tLOWESS fitting to reference data complete.\n');
Y_target = 1;
% Calculate repetitiveness_bias_corrected length_bia_corrected ave_read_count data for plotting and later analysis.
Y_fitCurve2_parent          = interp1(newX2_parent,newY2_parent,X_GCbias,'spline');
fprintf(['& length(Y_reads_corrected_1_parent) = ' num2str(length(Y_reads_corrected_1_parent)) '\n']);
fprintf(['& length(Y_fitCurve2_parent)         = ' num2str(length(Y_fitCurve2_parent        )) '\n']);
Y_reads_corrected_2_parent  = Y_reads_corrected_1_parent./Y_fitCurve2_parent*Y_target;


%-------------------------------------------------------------------------------------------------
% Add corrected read depths to fragment data structure.
%-------------------------------------------------------------------------------------------------
for fragID = 1:numFragments
	fragment_data(fragID).reads_ave_corrected_2        = Y_reads_corrected_2(fragID);
	fragment_data(fragID).reads_ave_corrected_2_parent = Y_reads_corrected_2_parent(fragID);
end;


%-------------------------------------------------------------------------------------------------
% Send status to log file.
%-------------------------------------------------------------------------------------------------
fprintf(['\tlength(X_length)                   = ' num2str(length(X_length))                   '\n']);
fprintf(['\tlength(Y_reads_project)            = ' num2str(length(Y_reads_project))            '\n']);
fprintf(['\tlength(Y_reads_parent)             = ' num2str(length(Y_reads_parent))             '\n']);
fprintf(['\tlength(Y_reads_corrected_1)        = ' num2str(length(Y_reads_corrected_1))        '\n']);
fprintf(['\tlength(Y_reads_corrected_1_parent) = ' num2str(length(Y_reads_corrected_1_parent)) '\n']);
fprintf(['\tlength(Y_reads_corrected_2)        = ' num2str(length(Y_reads_corrected_2))        '\n']);
fprintf(['\tlength(Y_reads_corrected_2_parent) = ' num2str(length(Y_reads_corrected_2_parent)) '\n']);
fprintf(['\tlength(X_GCbias)                   = ' num2str(length(X_GCbias))                   '\n']);
fprintf('\n');


%-------------------------------------------------------------------------------------------------
% Show GC bias in corrected' read depth for project.
%-------------------------------------------------------------------------------------------------
fprintf('Subplot 5/11 : [EXP] (Corrected reads 1) vs. (Fragment GC ratio).\n');
sh(5) = subplot(5,4,9);
hold on;
for frag = 1:numFragments
    if (fragment_data(frag).usable == 1)
        if (fragment_data(frag).reads_ave_corrected_1 < Y_fitCurve2(frag))
            plot(fragment_data(frag).GC_bias,fragment_data(frag).reads_ave_corrected_1,'.', 'color', [0.0, 0.66667, 0.33333], 'MarkerSize',1);
        else
            plot(fragment_data(frag).GC_bias,fragment_data(frag).reads_ave_corrected_1,'.', 'color', [0.0, 0.45, 0.55], 'MarkerSize',1);
        end;
    end;
end;
plot(newX2_project,newY2_project,'color','k','linestyle','-', 'linewidth',1);
hold off;
axis normal;
h = title('GC ratio vs. Corrected reads 1. [EXP]');   set(h, 'FontSize', 5);
h = ylabel('Corrected reads 1');  set(h, 'FontSize', 5);
h = xlabel('GC ratio');     set(h, 'FontSize', 5);
set(gca,'FontSize',5);
xlim([0 1]);
ylim([0 4]);

if (ParentValid == 1)
	%-------------------------------------------------------------------------------------------------
	% Show GC bias in corrected' read depth for parent.
	%-------------------------------------------------------------------------------------------------
	fprintf('Subplot 6/11 : [REF] (Corrected reads 1) vs. (Fragment GC ratio).\n');
	sh(7) = subplot(5,4,11);
	hold on;
	for frag = 1:numFragments
	    if (fragment_data(frag).usable_parent == 1)
	        if (fragment_data(frag).reads_ave_corrected_parent_1 < Y_fitCurve2_parent(frag))
	            plot(fragment_data(frag).GC_bias,fragment_data(frag).reads_ave_corrected_parent_1,'.', 'color', [0.0, 0.66667, 0.33333], 'MarkerSize',1);
	        else
	            plot(fragment_data(frag).GC_bias,fragment_data(frag).reads_ave_corrected_parent_1,'.', 'color', [0.0, 0.45, 0.55], 'MarkerSize',1);
	        end;
	    end;
	end;
	plot(newX2_parent,newY2_parent,'color','k','linestyle','-', 'linewidth',1);
	hold off;
	axis normal;
	h = title('GC ratio vs. Corrected reads 1. [REF]');   set(h, 'FontSize', 5);
	h = ylabel('Corrected reads 1');  set(h, 'FontSize', 5);
	h = xlabel('GC ratio');     set(h, 'FontSize', 5);
	set(gca,'FontSize',5);
	xlim([0 1]);
	ylim([0 4]);
end;


%-------------------------------------------------------------------------------------------------
% Show GC values in corrected'' read depth for project.
%-------------------------------------------------------------------------------------------------
fprintf('Subplot 7/11 : [EXP] (Corrected reads 2) vs. (Fragment GC ratio).\n');
sh(6) = subplot(5,4,10);
hold on;
for frag = 1:numFragments
    % Plot each corrected data point, colored depending on if data is useful or not.
    if (fragment_data(frag).usable == 1)
        if (fragment_data(frag).reads_ave_corrected_2 < 1)
            plot(fragment_data(frag).GC_bias,fragment_data(frag).reads_ave_corrected_2,'.', 'color', [0.0, 0.66667, 0.33333], 'MarkerSize',1);
        else
            plot(fragment_data(frag).GC_bias,fragment_data(frag).reads_ave_corrected_2,'.', 'color', [0.0, 0.45, 0.55], 'MarkerSize',1);
        end;
    end;
    if (mod(frag,1000) == 0)
        fprintf('.');
    end;
end;
fprintf('\n');
plot(newX2_project,Y_target,'color','k','linestyle','-', 'linewidth',1);
hold off;
axis normal;
h = title('GC ratio bias corrected.');   set(h, 'FontSize', 5);
h = ylabel('Corrected reads 2'); set(h, 'FontSize', 5);
h = xlabel('GC ratio');    set(h, 'FontSize', 5);
set(gca,'FontSize',5);
xlim([0 1]);
ylim([0 4]);

%-------------------------------------------------------------------------------------------------
% Show GC bias in corrected'' read depth for parent.
%-------------------------------------------------------------------------------------------------
if (ParentValid == 1)
	fprintf('Subplot 8/11 : [REF] (Corrected reads 2) vs. (Fragment GC ratio).\n');
	sh(8) = subplot(5,4,12);
	hold on;
	for fragID = 1:numFragments
	    % Plot each corrected data point, colored depending on if data is useful or not.
	    if (fragment_data(fragID).usable_parent == 1)
	        if (fragment_data(fragID).reads_ave_corrected_2_parent < 1)
	            plot(fragment_data(fragID).GC_bias,fragment_data(fragID).reads_ave_corrected_2_parent,'.', 'color', [0.0, 0.66667, 0.33333], 'MarkerSize',1);
	        else
	            plot(fragment_data(fragID).GC_bias,fragment_data(fragID).reads_ave_corrected_2_parent,'.', 'color', [0.0, 0.45, 0.55], 'MarkerSize',1);
	        end;
	    end;
	    if (mod(fragID,1000) == 0)
	        fprintf('.');
	    end;
	end;
	fprintf('\n');
	plot(newX2_parent,Y_target,'color','k','linestyle','-', 'linewidth',1);
	hold off;
	axis normal;
	h = title('GC ratio bias corrected.');   set(h, 'FontSize', 5);
	h = ylabel('Corrected reads 2'); set(h, 'FontSize', 5);
	h = xlabel('GC ratio');    set(h, 'FontSize', 5);
	set(gca,'FontSize',5);
	xlim([0 1]);
	ylim([0 4]);
end;

%-------------------------------------------------------------------------------------------------
% Show corrected'' read depth vs length for project.
%-------------------------------------------------------------------------------------------------
fprintf('Subplot 9/11 : [EXP] (Corrected reads 3) vs. (Fragment length).\n');
sh(9) = subplot(5,4,[13 14]);
hold on;
for fragID = 1:numFragments
    X_datum = fragment_data(fragID).length;
    Y_datum = fragment_data(fragID).reads_ave_corrected_2;

    % Plot each corrected data point, colored depending on if data is useful or not.
    if (fragment_data(fragID).usable == 1)
    plot(X_datum,Y_datum,'.', 'color', [0.0, 0.66667, 0.33333], 'MarkerSize',1);
    else
    plot(X_datum,Y_datum,'.', 'color', [1.0, 0.0, 0.0], 'MarkerSize',1);
    end;
end;
plot(newX1_project,1,'color','k','linestyle','-', 'linewidth',1);
hold off;
axis normal;
h = title('Corrected Read count vs. ddRADseq fragment length. [EXP]');   set(h, 'FontSize', 5);
h = ylabel('Corrected reads');   set(h, 'FontSize', 5);
h = xlabel('Fragment Length');   set(h, 'FontSize', 5);
set(gca,'FontSize',5);
ylim([0 4]);
xlim([1 1000]);

if (ParentValid == 1)
	%-------------------------------------------------------------------------------------------------
	% Show corrected'' read depth vs length for parent.
	%-------------------------------------------------------------------------------------------------
	fprintf('Subplot 10/10 : [REF] (Corrected reads 3) vs. (Fragment length).\n');
	sh(10) = subplot(5,4,[15 16]);
	hold on;
	for fragID = 1:numFragments
	    X_datum = fragment_data(fragID).length;
	    Y_datum = fragment_data(fragID).reads_ave_corrected_2_parent;
	
	    % Plot each corrected data point, colored depending on if data is useful or not.
	    if (fragment_data(fragID).usable_parent == 1)
	        plot(X_datum,Y_datum,'.', 'color', [0.0, 0.66667, 0.33333], 'MarkerSize',1);
	    else
	        plot(X_datum,Y_datum,'.', 'color', [1.0, 0.0, 0.0], 'MarkerSize',1);
	    end;
	end;
	plot(newX1_parent,1,'color','k','linestyle','-', 'linewidth',1);
	hold off;
	axis normal;
	h = title('Corrected Read count vs. ddRADseq fragment length. [REF]');   set(h, 'FontSize', 5);
	h = ylabel('Corrected reads');   set(h, 'FontSize', 5);
	h = xlabel('Fragment Length');   set(h, 'FontSize', 5);
	set(gca,'FontSize',5);
	ylim([0 4]);
	xlim([1 1000]);
end;

if (ParentValid == 1)
	%-------------------------------------------------------------------------------------------------
	% Subfigure showing corrected project data divided by corrected parent data.
	%-------------------------------------------------------------------------------------------------
	fprintf('Subplot 11/11 : [REF] (Corrected reads 3) vs. (Fragment length).\n');
	sh(11) = subplot(5,4,[18 19]);
	hold on;
	for fragID = 1:numFragments
		X_datum  = fragment_data(fragID).length;
		Y_datum1 = fragment_data(fragID).reads_ave_corrected_2;
		Y_datum2 = fragment_data(fragID).reads_ave_corrected_2_parent;

		% Plot each corrected data point, colored depending on if data is useful or not.
		if (fragment_data(fragID).usable_parent == 1)
			plot(X_datum,Y_datum1./Y_datum2,'.', 'color', [0.0, 0.66667, 0.33333], 'MarkerSize',1);
		else
			plot(X_datum,Y_datum1./Y_datum2,'.', 'color', [1.0, 0.0, 0.0], 'MarkerSize',1);
		end;
	end;
	hold off;
	axis normal;
	h = title('Final Corrected Read count vs. ddRADseq fragment length.');   set(h, 'FontSize', 5);
	h = ylabel('Corrected reads');   set(h, 'FontSize', 5);
	h = xlabel('Fragment Length');   set(h, 'FontSize', 5);
	set(gca,'FontSize',5);
	ylim([0 4]);
	xlim([1 1000]);
end;


%-------------------------------------------------------------------------------------------------
% Saving figure.
%-------------------------------------------------------------------------------------------------
fprintf('Saving figure.\n');
set(fig,'PaperPosition',[0 0 8 6]*2);
saveas(fig, [main_dir 'users/' user '/projects/' project '/fig.examine_bias.eps'], 'epsc');
saveas(fig, [main_dir 'users/' user '/projects/' project '/fig.examine_bias.png'], 'png');
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
save([main_dir 'users/' user '/projects/' project '/fragment_CNV_data.mat'], 'fragments_CNV');

%%=========================================================================
%%= Analyze corrected CNV data for RADseq datasets. =======================
%%=========================================================================
fprintf(['\nGenerating CNV figure from ''' project ''' sequence data corrected by removing restriction fragment length bias.\n']);
% Initializes vectors used to hold copy number data.
for chr = 1:num_chrs   % number of chrs.
	if (chr_in_use(chr) == 1)
		for j = 1:2       % 2 categories tracked : total read counts in bin and number of data entries in region (size of bin in base-pairs).
			chr_CNVdata_RADseq_project{chr,j} = zeros(1,ceil(chr_size(chr)/bases_per_bin));
			chr_CNVdata_RADseq_parent{chr,j}  = zeros(1,ceil(chr_size(chr)/bases_per_bin));
		end;
	end;
end;
% Output chromosome lengths to log file.
for i = 1:length(chr_name)
	fprintf(['\nchr' num2str(i) ' = ''' chr_name{i} '''.\tCGHlength = ' num2str(length(chr_CNVdata_RADseq_project{i,1}))]);
end;

chr_CNVdata_RADseq = chr_CNVdata_RADseq_project;



if (exist([main_dir 'users/' user '/projects/' project '/corrected_CNV.mat'],'file') == 0)
	fprintf('\nMAT file containing project CNV information not found, regenerating from prior data files.\n');
	% Convert bias corrected ave_copy_number per restriction digest fragment into CNV data for plotting
	fprintf('\n# Adding fragment corrected_ave_read_copy values to map bins.');
	for fragID = 1:numFragments
		if (fragment_data(fragID).usable == 1)
			% Load important data from fragments data structure.
			chrID      = fragment_data(fragID).chr;
			posStart   = fragment_data(fragID).startbp;
			posEnd     = fragment_data(fragID).endbp;
			count      = fragment_data(fragID).reads_ave_corrected_2;
			fragLength = fragment_data(fragID).length;

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

	save([main_dir 'users/' user '/projects/' project '/corrected_CNV.mat'],'chr_CNVdata_RADseq');
else
	fprintf('\nProject CNV MAT file found, loading.\n');
	load([main_dir 'users/' user '/projects/' project '/corrected_CNV.mat']);
end;
chr_CNVdata_RADseq_project = chr_CNVdata_RADseq;


chr_CNVdata_RADseq = chr_CNVdata_RADseq_parent;
fprintf('\nRefrence CNV MAT file found, loading.\n');
load([main_dir 'users/' parentUser '/projects/' parentProject '/corrected_CNV.mat']);
chr_CNVdata_RADseq_parent = chr_CNVdata_RADseq;

        
% basic plot parameters not defined per genome.
TickSize        = -0.005;  %negative for outside, percentage of longest chr figure.
maxY            = ploidyBase*2;

%% -----------------------------------------------------------------------------------------
% Make figures
%-------------------------------------------------------------------------------------------
fig = figure(1);
set(gcf, 'Position', [0 70 1024 600]);
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		for pos = 1:length(chr_CNVdata_RADseq{chr,1})
			% Plot the sum of the data in each region, divided by the number of data points in each region.
			% Then divided by this value calculated for SC5314 data.
			if (chr_CNVdata_RADseq_project{chr,2}(pos) == 0) || (chr_CNVdata_RADseq_parent{chr,2}(pos) == 0)
				% No data elements => null value is plotted.
				CNVplot{chr}(pos) = 0;
			else
				% Sum of data elements is divided by the number of data elements.
				project_factor = chr_CNVdata_RADseq_project{chr,1}(pos)/chr_CNVdata_RADseq_project{chr,2}(pos);
				parent_factor  = chr_CNVdata_RADseq_parent{chr,1}(pos)/chr_CNVdata_RADseq_parent{chr,2}(pos)
				if (ParentValid == 1)
					CNVplot{chr}(pos) = project_factor/parent_factor;
				else
					CNVplot{chr}(pos) = project_factor;
				end;
			end;
		end;
		% chr_max(chr) = max(CNVplot{chr});
		% chr_med(chr) = median(CNVplot{chr});
	end;
end;
% max_count     = max(chr_max);
% median_count  = sum(chr_med)/length(chr_med);
for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
		% CNVplot2{chr} = CNVplot{chr}/median_count;
		CNVplot2{chr} = CNVplot{chr};
	end;
end;

% Save presented CNV data in a file format common across data types being processed.
fprintf('\nSaving "Common_CNV" data file.');
genome_CNV = genome;
save([main_dir 'users/' user '/projects/' project '/Common_CNV.mat'], 'CNVplot2','genome_CNV');

ploidy = str2num(ploidyEstimate);
fprintf(['Ploidy string = "' ploidyEstimate '"\n']);
[chr_breaks, chrCopyNum, ploidyAdjust] = FindChrSizes_4(Aneuploidy,CNVplot2,ploidy,num_chrs,chr_in_use)
largestChr = find(chr_width == max(chr_width));


%% -----------------------------------------------------------------------------------------
% Output corrected CNV data to pileup-type file.
%-------------------------------------------------------------------------------------------
%% Output corrected reads, multiplied by ploidyEstimate, to : [main_dir 'users/' user '/projects/' project '/RADseq_digest_analysis_CNV_corrected.txt']
%  Including : [chrNum, bpStart, bpEnd, correctedReads*ploidy, fragmentLength]
%   1       9638    10115   4796    478
%   1       10116   10123   1283    8
%   1       13170   13841   55175   672
fprintf('\nOutputting fragment read copy after correction.\n');
outfile     = [main_dir 'users/' user '/projects/' project '/RADseq_digest_analysis_CNV_corrected.txt'];
outfile_fid = fopen(outfile,'w');
for fragID = 1:numFragments
	string1 = num2str(fragment_data(fragID).chr);
	string2 = num2str(fragment_data(fragID).startbp);
	string3 = num2str(fragment_data(fragID).endbp);
	if (Low_quality_ploidy_estimate == true)
		if (ParentValid == 1)
			string4 = num2str(fragment_data(fragID).reads_ave_corrected_2/fragment_data(fragID).reads_ave_corrected_2_parent*ploidy*ploidyAdjust);
		else
			string4 = num2str(fragment_data(fragID).reads_ave_corrected_2*ploidy*ploidyAdjust);
		end;
	else
		if (ParentValid == 1)
			string4 = num2str(fragment_data(fragID).reads_ave_corrected_2/fragment_data(fragID).reads_ave_corrected_2_parent*ploidy);
		else
			string4 = num2str(fragment_data(fragID).reads_ave_corrected_2*ploidy);
		end;
	end;
	string5 = num2str(fragment_data(fragID).length);
	if (fragment_data(fragID).usable == 1)
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
	Linear_fig           = figure(2);
	Linear_genome_size   = sum(chr_size);
	Linear_Chr_max_width = 0.85;
	Linear_left_start    = 0.07;
	Linear_left_chr_gap  = 0.01*7/(num_chrs-1);
	Linear_height        = 0.6;
	Linear_base          = 0.1;
	Linear_TickSize      = -0.01;  %negative for outside, percentage of longest chr figure.
	maxY                 = ploidyBase*2;
	Linear_left          = Linear_left_start;
end;


%% -----------------------------------------------------------------------------------------
% Median normalize CNV data before figure generation.
%-------------------------------------------------------------------------------------------
% Gather CGH data for LOWESS fitting.
CNVdata_all = [];
for chr = 1:num_chrs
    if (chr_in_use(chr) == 1)
        CNVdata_all = [CNVdata_all   CNVplot2{chr}];
    end;
end;
CNVdata_all(CNVdata_all == 0) = [];
medianCNV = median(CNVdata_all)
for chr = 1:num_chrs
    if (chr_in_use(chr) == 1)
        CNVplot2{chr} = CNVplot2{chr}/medianCNV;
    end;
end;


%% -----------------------------------------------------------------------------------------
% Make figures
%-------------------------------------------------------------------------------------------
% initialize string to contain chromosome copy number estimates for later output.
stringChrCNVs = '';

for chr = 1:num_chrs
	if (chr_in_use(chr) == 1)
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
		x2 = chr_size(chr)/bases_per_bin;
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
        
		% This section sets the Y-axis labelling.
		switch ploidyBase
			case 1
				set(gca,'YTick',[0 maxY/2 maxY]);
				set(gca,'YTickLabel',{'','',''});
				text(-50000/bases_per_bin, maxY/2,   '1','HorizontalAlignment','right','Fontsize',5);
				text(-50000/bases_per_bin, maxY,     '2','HorizontalAlignment','right','Fontsize',5);
			case 2
				set(gca,'YTick',[0 maxY/4 maxY/2 maxY/4*3 maxY]);
				set(gca,'YTickLabel',{'','','','',''});
				text(-50000/bases_per_bin, maxY/4,   '1','HorizontalAlignment','right','Fontsize',5);
				text(-50000/bases_per_bin, maxY/2,   '2','HorizontalAlignment','right','Fontsize',5);
				text(-50000/bases_per_bin, maxY/4*3, '3','HorizontalAlignment','right','Fontsize',5);
				text(-50000/bases_per_bin, maxY,     '4','HorizontalAlignment','right','Fontsize',5);
			case 3
				set(gca,'YTick',[0 maxY/6 maxY/3 maxY/2 maxY/3*2 maxY/6*5 maxY]);
				set(gca,'YTickLabel',{'','','','','','',''});
				text(-50000/bases_per_bin, maxY/2,   '3','HorizontalAlignment','right','Fontsize',5);
				text(-50000/bases_per_bin, maxY,     '6','HorizontalAlignment','right','Fontsize',5);
			case 4
				set(gca,'YTick',[0 maxY/8 maxY/4 maxY/8*3 maxY/2 maxY/8*5 maxY/4*3 maxY/8*7 maxY]);
				set(gca,'YTickLabel',{'','','','','','','','',''});
				text(-50000/bases_per_bin, maxY/4,   '2','HorizontalAlignment','right','Fontsize',5);
				text(-50000/bases_per_bin, maxY/2,   '4','HorizontalAlignment','right','Fontsize',5);
				text(-50000/bases_per_bin, maxY/4*3, '6','HorizontalAlignment','right','Fontsize',5);
				text(-50000/bases_per_bin, maxY,     '8','HorizontalAlignment','right','Fontsize',5);
		end;

		set(gca,'FontSize',6);
		if (chr == find(chr_posY == max(chr_posY)))
			title([ project ' CNV map'],'Interpreter','none','FontSize',12);
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
		x1 = cen_start(chr)/bases_per_bin;
		x2 = cen_end(chr)/bases_per_bin;
		leftEnd  = 0.5*(5000/bases_per_bin);
		rightEnd = chr_size(chr)/bases_per_bin-0.5*(5000/bases_per_bin);
		if (Centromere_format == 0)
			% standard chromosome cartoons in a way which will not cause segfaults when running via commandline.
			dx = cen_tel_Xindent;
			dy = cen_tel_Yindent;
			% draw white triangles at corners and centromere locations.
			fill([leftEnd   leftEnd   leftEnd+dx ],       [maxY-dy   maxY      maxY],         [1.0 1.0 1.0], 'EdgeColor',[1.0 1.0 1.0]);    % top left corner.
			fill([leftEnd   leftEnd   leftEnd+dx ],       [dy        0         0   ],         [1.0 1.0 1.0], 'EdgeColor',[1.0 1.0 1.0]);    % bottom left corner.
			fill([rightEnd  rightEnd  rightEnd-dx],       [maxY-dy   maxY      maxY],         [1.0 1.0 1.0], 'EdgeColor',[1.0 1.0 1.0]);    % top right corner.
			fill([rightEnd  rightEnd  rightEnd-dx],       [dy        0         0   ],         [1.0 1.0 1.0], 'EdgeColor',[1.0 1.0 1.0]);    % bottom right corner.
			fill([x1-dx     x1        x2           x2+dx],[maxY      maxY-dy   maxY-dy  maxY],[1.0 1.0 1.0], 'EdgeColor',[1.0 1.0 1.0]);    % top centromere.
			fill([x1-dx     x1        x2           x2+dx],[0         dy        dy       0   ],[1.0 1.0 1.0], 'EdgeColor',[1.0 1.0 1.0]);    % bottom centromere.
			% draw outlines of chromosome cartoon.   (drawn after horizontal lines to that cartoon edges are not interrupted by horiz lines.
			plot([leftEnd   leftEnd   leftEnd+dx   x1-dx   x1        x2        x2+dx   rightEnd-dx   rightEnd   rightEnd   rightEnd-dx x2+dx   x2   x1   x1-dx   leftEnd+dx   leftEnd],...
			     [dy        maxY-dy   maxY         maxY    maxY-dy   maxY-dy   maxY    maxY          maxY-dy    dy         0           0       dy   dy   0       0            dy     ],...
			      'Color',[0 0 0]);
		elseif (Centromere_format == 1)
			leftEnd  = 0;
			rightEnd = chr_size(chr)/bases_per_bin;

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
			else
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
			x2 = chr_size(chr)/bases_per_bin;
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
			x1 = cen_start(chr)/bases_per_bin;
			x2 = cen_end(chr)/bases_per_bin;
			leftEnd  = 0.5*(5000/bases_per_bin);
			rightEnd = chr_size(chr)/bases_per_bin-0.5*(5000/bases_per_bin);
			if (Centromere_format == 0)
				% standard chromosome cartoons in a way which will not cause segfaults when running via commandline.
				dx = cen_tel_Xindent;
				dy = cen_tel_Yindent;
				% draw white triangles at corners and centromere locations.
				fill([leftEnd   leftEnd   leftEnd+dx ],       [maxY-dy   maxY      maxY],         [1.0 1.0 1.0], 'EdgeColor',[1.0 1.0 1.0]);    % top left corner.
				fill([leftEnd   leftEnd   leftEnd+dx ],       [dy        0         0   ],         [1.0 1.0 1.0], 'EdgeColor',[1.0 1.0 1.0]);    % bottom left corner.
				fill([rightEnd  rightEnd  rightEnd-dx],       [maxY-dy   maxY      maxY],         [1.0 1.0 1.0], 'EdgeColor',[1.0 1.0 1.0]);    % top right corner.
				fill([rightEnd  rightEnd  rightEnd-dx],       [dy        0         0   ],         [1.0 1.0 1.0], 'EdgeColor',[1.0 1.0 1.0]);    % bottom right corner.
				fill([x1-dx     x1        x2           x2+dx],[maxY      maxY-dy   maxY-dy  maxY],[1.0 1.0 1.0], 'EdgeColor',[1.0 1.0 1.0]);    % top centromere.
				fill([x1-dx     x1        x2           x2+dx],[0         dy        dy       0   ],[1.0 1.0 1.0], 'EdgeColor',[1.0 1.0 1.0]);    % bottom centromere.
				% draw outlines of chromosome cartoon.   (drawn after horizontal lines to that cartoon edges are not interrupted by horiz lines.
				plot([leftEnd   leftEnd   leftEnd+dx   x1-dx   x1        x2        x2+dx   rightEnd-dx   rightEnd   rightEnd   rightEnd-dx   x2+dx   x2   x1   x1-dx   leftEnd+dx  leftEnd],...
				     [dy        maxY-dy   maxY         maxY    maxY-dy   maxY-dy   maxY    maxY          maxY-dy    dy         0             0       dy   dy   0       0           dy     ],...
				     'Color',[0 0 0]);
			elseif (Centromere_format == 1)
				leftEnd  = 0;
				rightEnd = chr_size(chr)/bases_per_bin;
        
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
			if (chr == 1)
				ylabel(project, 'Rotation', 0, 'HorizontalAlign', 'right', 'VerticalAlign', 'bottom','Interpreter','none','FontSize',5);

				% This section sets the Y-axis labelling.
				switch ploidyBase
					case 1
						set(gca,'YTick',[0 maxY/2 maxY]);
						set(gca,'YTickLabel',{'','',''});
						text(-50000/bases_per_bin, maxY/2,   '1','HorizontalAlignment','right','Fontsize',5);
						text(-50000/bases_per_bin, maxY,     '2','HorizontalAlignment','right','Fontsize',5);
					case 2
						set(gca,'YTick',[0 maxY/4 maxY/2 maxY/4*3 maxY]);
						set(gca,'YTickLabel',{'','','','',''});
						text(-50000/bases_per_bin, maxY/4,   '1','HorizontalAlignment','right','Fontsize',5);
						text(-50000/bases_per_bin, maxY/2,   '2','HorizontalAlignment','right','Fontsize',5);
						text(-50000/bases_per_bin, maxY/4*3, '3','HorizontalAlignment','right','Fontsize',5);
						text(-50000/bases_per_bin, maxY,     '4','HorizontalAlignment','right','Fontsize',5);
					case 3
						set(gca,'YTick',[0 maxY/6 maxY/3 maxY/2 maxY/3*2 maxY/6*5 maxY]);
						set(gca,'YTickLabel',{'','','','','','',''});
						text(-50000/bases_per_bin, maxY/2,   '3','HorizontalAlignment','right','Fontsize',5);
						text(-50000/bases_per_bin, maxY,     '6','HorizontalAlignment','right','Fontsize',5);
					case 4
						set(gca,'YTick',[0 maxY/8 maxY/4 maxY/8*3 maxY/2 maxY/8*5 maxY/4*3 maxY/8*7 maxY]);
						set(gca,'YTickLabel',{'','','','','','','','',''});
						text(-50000/bases_per_bin, maxY/4,   '2','HorizontalAlignment','right','Fontsize',5);
						text(-50000/bases_per_bin, maxY/2,   '4','HorizontalAlignment','right','Fontsize',5);
						text(-50000/bases_per_bin, maxY/4*3, '6','HorizontalAlignment','right','Fontsize',5);
						text(-50000/bases_per_bin, maxY,     '8','HorizontalAlignment','right','Fontsize',5);
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
end;
            


% Make figure output have transparent background.
set(gcf, 'color', 'none',...
         'inverthardcopy', 'off');
            
% Save original arrangement of chromosomes.
set(fig,'PaperPosition',[0 0 8 6]*2);
saveas(fig, [main_dir 'users/' user '/projects/' project '/fig.CNV-map.1.eps'], 'epsc');
saveas(fig, [main_dir 'users/' user '/projects/' project '/fig.CNV-map.1.png'], 'png');
delete(fig);
    
% Save horizontally arranged chromosomes.
set(Linear_fig,'PaperPosition',[0 0 8 0.62222222]*2);
saveas(Linear_fig, [main_dir 'users/' user '/projects/' project '/fig.CNV-map.2.eps'], 'epsc');
saveas(Linear_fig, [main_dir 'users/' user '/projects/' project '/fig.CNV-map.2.png'], 'png');
delete(Linear_fig);

% Output chromosome copy number estimates.
textFileName = [main_dir 'users/' user '/projects/' project '/CNV-map.3.txt'];
fprintf(['Text output of CNVs : "' textFileName '"\n']);
textFileID = fopen(textFileName,'w');
fprintf(textFileID,stringChrCNVs);
fclose(textFileID);

end
