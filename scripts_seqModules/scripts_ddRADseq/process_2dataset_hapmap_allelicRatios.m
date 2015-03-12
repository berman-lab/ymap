function [] = process_2dataset_hapmap_allelicRatios(project1dir, project2dir, hapmapDir, chr_size, chr_name, chr_in_use, SNP_verString);


%%============================================================================================================
% Define files for processing.
%-------------------------------------------------------------------------------------------------------------
%	% Child data at loci where parent has allelic ratio on range [0.25-0.75].
%	C_datafile = [project1dir 'trimmed_SNPs_v4.txt'       ];


% Parent data at loci where allelic ratio is on range [0.25-0.75].
P_datafile = [project1dir 'trimmed_SNPs_v4.parent.txt'];

% Child data at loci from hapmap.
C_datafile = [project1dir 'trimmed_SNPs_v5.txt'       ];

% Hapmap dataset.
H_datafile = [hapmapDir   'SNPdata_parent.txt'        ];


%%============================================================================================================
% Preallocate data vectors the length of each chromosome.
%-------------------------------------------------------------------------------------------------------------
fprintf(['process_2dataset_hapmap_allelicRatios.m: Preallocate data vectors.\n']);
fprintf('Preallocate data vectors for each chromosome.\n');
C_chr_SNP_data_positions = cell(length(chr_size),1);   % coordinate of SNP.
C_chr_SNP_data_ratios    = cell(length(chr_size),1);   % allelic ratio of SNP.
C_chr_baseCall           = cell(length(chr_size),1);   % majority basecall of SNP.
C_chr_count              = cell(length(chr_size),1);   % number of reads at SNP coordinate.
C_chr_SNP_homologA       = cell(length(chr_size),1);   % hapmap homolog a basecall.
C_chr_SNP_homologB       = cell(length(chr_size),1);   % hapmap homolog b basecall.
C_chr_SNP_flipHomologs   = cell(length(chr_size),1);   % does hapmap entry need flipped?
C_chr_SNP_keep           = cell(length(chr_size),1);   % is SNP coordinate found in hapmap?

P_chr_SNP_data_positions = cell(length(chr_size),1);
P_chr_SNP_data_ratios    = cell(length(chr_size),1);
P_chr_baseCall           = cell(length(chr_size),1);
P_chr_count              = cell(length(chr_size),1);
P_chr_SNP_homologA       = cell(length(chr_size),1);
P_chr_SNP_homologB       = cell(length(chr_size),1);
P_chr_SNP_flipHomologs   = cell(length(chr_size),1);
P_chr_SNP_keep           = cell(length(chr_size),1);

H_chr_SNP_data_positions = cell(length(chr_size),1);
H_chr_SNP_alleleA        = cell(length(chr_size),1);
H_chr_SNP_alleleB        = cell(length(chr_size),1);
H_chr_SNP_hapmapEntry    = cell(length(chr_size),1);

for chrID = 1:length(chr_size)
	if (chr_in_use(chrID) == 1)
		C_chr_SNP_data_positions{chrID} = zeros(chr_size(chrID),1);
		C_chr_SNP_data_ratios{   chrID} = zeros(chr_size(chrID),1);
		C_chr_count{             chrID} = zeros(chr_size(chrID),1);
		C_chr_baseCall{          chrID} = cell( chr_size(chrID),1);
		C_chr_SNP_homologA{      chrID} = cell( chr_size(chrID),1);
		C_chr_SNP_homologB{      chrID} = cell( chr_size(chrID),1);
		C_chr_SNP_flipHomologs{  chrID} = zeros(chr_size(chrID),1);
		C_chr_SNP_keep{          chrID} = ones( chr_size(chrID),1);
		C_chr_lines_analyzed(    chrID) = 0;

		P_chr_SNP_data_positions{chrID} = zeros(chr_size(chrID),1);
		P_chr_SNP_data_ratios{   chrID} = zeros(chr_size(chrID),1);
		P_chr_count{             chrID} = zeros(chr_size(chrID),1);
		P_chr_baseCall{          chrID} = cell( chr_size(chrID),1);
		P_chr_SNP_homologA{      chrID} = cell( chr_size(chrID),1);
		P_chr_SNP_homologB{      chrID} = cell( chr_size(chrID),1);
		P_chr_SNP_flipHomologs{  chrID} = zeros(chr_size(chrID),1);
		P_chr_SNP_keep{          chrID} = ones( chr_size(chrID),1);
		P_chr_lines_analyzed(    chrID) = 0;

		H_chr_SNP_data_positions{chrID} = zeros(chr_size(chrID),1);
		H_chr_SNP_alleleA{       chrID} = cell( chr_size(chrID),1);
		H_chr_SNP_alleleB{       chrID} = cell( chr_size(chrID),1);
		H_chr_SNP_hapmapEntry{   chrID} = zeros(chr_size(chrID),1);
		H_chr_lines_analyzed(    chrID) = 0;
	end;
end;


%%============================================================================================================
% Process parent project dataset.
%-------------------------------------------------------------------------------------------------------------
fprintf(['process_2dataset_hapmap_allelicRatios.m: Process parent project dataset (trimmed_SNPs_v4.parent.txt).\n']);
P_data      = fopen(P_datafile, 'r');
allele_list = ['A' 'T' 'G' 'C'];
old_chr     = 0;
while not (feof(P_data))
	P_dataLine = fgetl(P_data);
	if (P_dataLine(1) ~= '#') 
		if (length(P_dataLine) > 0)
			% process the loaded line into data channels.
			values           = strsplit(strtrim(P_dataLine),'	');
			P_SNP_chr_name   = values{1};
			P_SNP_coordinate = values{2};
			P_SNP_countA     = values{3};
			P_SNP_countT     = values{4};
			P_SNP_countG     = values{5};
			P_SNP_countC     = values{6};
			P_chr_num        = find(strcmp(P_SNP_chr_name, chr_name));
			if (length(P_chr_num) > 0)
				if (P_chr_num ~= old_chr)
					fprintf(['\tchr = ' num2str(P_chr_num) '\n']);
				end;
				P_SNP_countA       = str2num(P_SNP_countA);
				P_SNP_countT       = str2num(P_SNP_countT);
				P_SNP_countG       = str2num(P_SNP_countG);
				P_SNP_countC       = str2num(P_SNP_countC);
				P_count_vector1    = [P_SNP_countA P_SNP_countT P_SNP_countG P_SNP_countC];
				P_chr_read_max1    = max(P_count_vector1);
				P_SNP_coordinate   = str2num(P_SNP_coordinate);
				P_chr_lines_analyzed(P_chr_num) = P_chr_lines_analyzed(P_chr_num)+1;
				P_chr_SNP_data_positions{P_chr_num}(P_chr_lines_analyzed(P_chr_num)) = P_SNP_coordinate;
				P_chr_SNP_data_ratios{   P_chr_num}(P_chr_lines_analyzed(P_chr_num)) = P_chr_read_max1/sum(P_count_vector1);
				P_chr_count{             P_chr_num}(P_chr_lines_analyzed(P_chr_num)) = sum(P_count_vector1);
				allele_call_id = find(P_count_vector1==max(P_count_vector1));
				if (length(allele_call_id) > 1)
					P_chr_read_id = 'N';
				else
					P_chr_read_id = allele_list(allele_call_id);
				end;
				P_chr_baseCall{          P_chr_num}{P_chr_lines_analyzed(P_chr_num)} = P_chr_read_id;
				old_chr = P_chr_num;
			else
				old_chr = 0;
			end;
		end;
	end;
end;
fclose(P_data);


%%============================================================================================================
% Process child project dataset.
%-------------------------------------------------------------------------------------------------------------
fprintf(['process_2dataset_hapmap_allelicRatios.m: Process child project dataset (trimmed_SNPs_v5.txt).\n']);
C_data      = fopen(C_datafile, 'r');
allele_list = ['A' 'T' 'G' 'C'];
old_chr     = 0;
while not (feof(C_data))
	C_dataLine = fgetl(C_data);
	if (C_dataLine(1) ~= '#')
		if (length(C_dataLine) > 0)
			% process the loaded line into data channels.
			values           = strsplit(strtrim(C_dataLine),'	');
			C_SNP_chr_name   = values{1};
			C_SNP_coordinate = values{2};
			C_SNP_countA     = values{3};
			C_SNP_countT     = values{4};
			C_SNP_countG     = values{5};
			C_SNP_countC     = values{6};
			C_chr_num        = find(strcmp(C_SNP_chr_name, chr_name));
			if (length(C_chr_num) > 0)
				if (C_chr_num ~= old_chr);   fprintf(['\tchr = ' num2str(C_chr_num) '\n']);   end;
				C_SNP_countA                                                         = str2num(C_SNP_countA);
				C_SNP_countT                                                         = str2num(C_SNP_countT);
				C_SNP_countG                                                         = str2num(C_SNP_countG);
				C_SNP_countC                                                         = str2num(C_SNP_countC);
				C_count_vector1                                                      = [C_SNP_countA C_SNP_countT C_SNP_countG C_SNP_countC];
				C_chr_read_max1                                                      = max(C_count_vector1);
				C_SNP_coordinate                                                     = str2num(C_SNP_coordinate);
				C_chr_lines_analyzed(C_chr_num)                                      = C_chr_lines_analyzed(C_chr_num)+1;
				C_chr_SNP_data_positions{C_chr_num}(C_chr_lines_analyzed(C_chr_num)) = C_SNP_coordinate;
				C_chr_SNP_data_ratios{   C_chr_num}(C_chr_lines_analyzed(C_chr_num)) = C_chr_read_max1/sum(C_count_vector1);
				C_chr_count{             C_chr_num}(C_chr_lines_analyzed(C_chr_num)) = sum(C_count_vector1);
				allele_call_id                                                       = find(C_count_vector1==max(C_count_vector1));
				if (length(allele_call_id) > 1)
					C_chr_read_id                                                = 'N';
				else
					C_chr_read_id                                                = allele_list(allele_call_id);
				end;
				C_chr_baseCall{          C_chr_num}{C_chr_lines_analyzed(C_chr_num)} = C_chr_read_id;
				old_chr = C_chr_num;
				
			else
				old_chr = 0;
			end;
		end;
	end;
end;
fclose(C_data);


%%============================================================================================================
% Process hapmap dataset.
%-------------------------------------------------------------------------------------------------------------
fprintf(['process_2dataset_hapmap_allelicRatios.m: Process hapmap dataset ([Hapmap_dir]/SNPdata_parent.txt).\n']);
H_data      = fopen(H_datafile, 'r');
old_chr     = 0;
while not (feof(H_data))
	H_dataLine = fgetl(H_data);
	if (H_dataLine(1) ~= '#')
		if (length(H_dataLine) > 0)
			% process the loaded line into data channels.
			values                   = strsplit(strtrim(H_dataLine),'	');
			H_SNP_chr_name           = values{1};
			H_SNP_coordinate         = str2num(values{2});
			H_SNP_alleleA            = values{3};
			H_SNP_alleleB            = values{4};
			H_chr_num                = find(strcmp(H_SNP_chr_name, chr_name));
			% darren: Make hapmapEntry consensus determination.
			H_SNP_hapmapEntry_values = [];
			for i = 5:length(values)
				H_SNP_hapmapEntry_values = [H_SNP_hapmapEntry_values str2num(values{i})];
			end;
			count_0     = length(find(H_SNP_hapmapEntry_values == 0 ));   % 0  : correct phase.
			count_1     = length(find(H_SNP_hapmapEntry_values == 1 ));   % 1  : incorrect phase.
			count_10    = length(find(H_SNP_hapmapEntry_values == 10));   % 10 : no phase information, due to no data at coodinate in child dataset.
			count_11    = length(find(H_SNP_hapmapEntry_values == 11));   % 11 : no phase information, due to coordinate not matching a LOH fragment in child dataset.
			count_12    = length(find(H_SNP_hapmapEntry_values == 12));   % 12 : no phase information, due to surprise allele in child dataset.
			count_error = count_10 + count_11 + count_12;
			if (count_0 > count_1)       H_SNP_hapmapEntry_consensus = 0;               % More evidence for correct phasing.
			elseif (count_1 > count_0)   H_SNP_hapmapEntry_consensus = 1;               % More evidence for incorrect phasing.
			elseif (count_0 == 0)        H_SNP_hapmapEntry_consensus = 10;              % No phasing information. All error codes will be treated the same.
			else                         H_SNP_hapmapEntry_consensus = round(rand());   % Equal non-zero evidence for correct and incorrect phasing, so choose 0 or 1 at random.
			end;
			if (length(H_chr_num) > 0)
				if (H_chr_num ~= old_chr)
					fprintf(['\tchr = ' num2str(H_chr_num) '\n']);
				end;
				H_chr_lines_analyzed(    H_chr_num)                                  = H_chr_lines_analyzed(H_chr_num)+1;
				H_chr_SNP_data_positions{H_chr_num}(H_chr_lines_analyzed(H_chr_num)) = H_SNP_coordinate;
				H_chr_SNP_alleleA{       H_chr_num}{H_chr_lines_analyzed(H_chr_num)} = H_SNP_alleleA;
				H_chr_SNP_alleleB{       H_chr_num}{H_chr_lines_analyzed(H_chr_num)} = H_SNP_alleleB;
				H_chr_SNP_hapmapEntry{   H_chr_num}(H_chr_lines_analyzed(H_chr_num)) = H_SNP_hapmapEntry_consensus;
				old_chr = H_chr_num;
			else
				old_chr = 0;
			end;
		end;
	end;
end;
fclose(H_data);


%%============================================================================================================
% Clean up data vectors.
%-------------------------------------------------------------------------------------------------------------
fprintf(['process_2dataset_hapmap_allelicRatios.m: clean up data.\n']);
for chrID = 1:length(chr_size)
	if (chr_in_use(chrID) == 1)
		C_chr_SNP_data_ratios{   chrID}(C_chr_SNP_data_positions{chrID} == 0)  = [];
		C_chr_count{             chrID}(C_chr_SNP_data_positions{chrID} == 0)  = [];
		C_chr_baseCall{          chrID}(C_chr_SNP_data_positions{chrID} == 0)  = [];
		C_chr_SNP_homologA{      chrID}(C_chr_SNP_data_positions{chrID} == 0)  = [];
		C_chr_SNP_homologB{      chrID}(C_chr_SNP_data_positions{chrID} == 0)  = [];
		C_chr_SNP_flipHomologs{  chrID}(C_chr_SNP_data_positions{chrID} == 0)  = [];
		C_chr_SNP_keep{          chrID}(C_chr_SNP_data_positions{chrID} == 0)  = [];
		C_chr_SNP_data_positions{chrID}(C_chr_SNP_data_positions{chrID} == 0)  = [];

		P_chr_SNP_data_ratios{   chrID}(P_chr_SNP_data_positions{chrID} == 0)  = [];
		P_chr_count{             chrID}(P_chr_SNP_data_positions{chrID} == 0)  = [];
		P_chr_baseCall{          chrID}(P_chr_SNP_data_positions{chrID} == 0)  = [];
		P_chr_SNP_homologA{      chrID}(P_chr_SNP_data_positions{chrID} == 0)  = [];
		P_chr_SNP_homologB{      chrID}(P_chr_SNP_data_positions{chrID} == 0)  = [];
		P_chr_SNP_flipHomologs{  chrID}(P_chr_SNP_data_positions{chrID} == 0)  = [];
		P_chr_SNP_keep{          chrID}(P_chr_SNP_data_positions{chrID} == 0)  = [];
		P_chr_SNP_data_positions{chrID}(P_chr_SNP_data_positions{chrID} == 0)  = [];

		H_chr_SNP_alleleA{       chrID}(H_chr_SNP_data_positions{chrID} == 0)  = [];
		H_chr_SNP_alleleB{       chrID}(H_chr_SNP_data_positions{chrID} == 0)  = [];
		H_chr_SNP_hapmapEntry{   chrID}(H_chr_SNP_data_positions{chrID} == 0)  = [];
		H_chr_SNP_data_positions{chrID}(H_chr_SNP_data_positions{chrID} == 0)  = [];
	end;
end;


%%============================================================================================================
% Determine child dataset values at hapmap loci.
%-------------------------------------------------------------------------------------------------------------
fprintf(['process_2dataset_hapmap_allelicRatios.m: Determine child values at hapmap loci.\n']);
start = 1;
for chrID = 1:length(chr_size)
	if (chr_in_use(chrID) == 1)
		fprintf(['\tchr = ' num2str(chrID) '\n']);
		for projectDatumID = 1:length(C_chr_SNP_data_positions{chrID})
			pos   = C_chr_SNP_data_positions{chrID}(projectDatumID);
			found = false;
			for hapmapDatumID = start:length(H_chr_SNP_data_positions{chrID})
				hapmap_pos = H_chr_SNP_data_positions{chrID}(hapmapDatumID);
				if (pos == hapmap_pos)
					found = true;
					break;
				end;
			end;
			if (found == true)
				C_chr_SNP_homologA{    chrID}{projectDatumID} = H_chr_SNP_alleleA{       chrID}{hapmapDatumID};;
				C_chr_SNP_homologB{    chrID}{projectDatumID} = H_chr_SNP_alleleB{       chrID}{hapmapDatumID};;
				C_chr_SNP_flipHomologs{chrID}(projectDatumID) = H_chr_SNP_hapmapEntry{   chrID}(hapmapDatumID);;
				C_chr_SNP_keep{        chrID}(projectDatumID) = 1;
				start = projectDatumID;
			else
				C_chr_SNP_keep{        chrID}(projectDatumID) = 0;
				start = 1;
			end;
		end;

		% Clean up data for chromosome.
		C_chr_SNP_data_ratios{   chrID}(C_chr_SNP_keep{chrID} == 0) = [];
		C_chr_SNP_data_positions{chrID}(C_chr_SNP_keep{chrID} == 0) = [];
		C_chr_baseCall{          chrID}(C_chr_SNP_keep{chrID} == 0) = [];
		C_chr_SNP_homologA{      chrID}(C_chr_SNP_keep{chrID} == 0) = [];
		C_chr_SNP_homologB{      chrID}(C_chr_SNP_keep{chrID} == 0) = [];
		C_chr_SNP_flipHomologs{  chrID}(C_chr_SNP_keep{chrID} == 0) = [];
		C_chr_count{             chrID}(C_chr_SNP_keep{chrID} == 0) = [];
		C_chr_SNP_keep{          chrID}(C_chr_SNP_keep{chrID} == 0) = [];
	end;
end;


%%============================================================================================================
% Determine parent dataset values at hapmap loci.
%-------------------------------------------------------------------------------------------------------------
fprintf(['process_2dataset_hapmap_allelicRatios.m: Determine parent values at hapmap loci.\n']);
start = 1;
for chrID = 1:length(chr_size)
	if (chr_in_use(chrID) == 1)
		fprintf(['\tchr = ' num2str(chrID) '\n']);
		for projectDatumID = 1:length(P_chr_SNP_data_positions{chrID})
			pos   = P_chr_SNP_data_positions{chrID}(projectDatumID);
			found = false;
			for hapmapDatumID = start:length(H_chr_SNP_data_positions{chrID})
				hapmap_pos = H_chr_SNP_data_positions{chrID}(hapmapDatumID);
				if (pos == hapmap_pos)
					found = true;
					break;
				end;
			end;
			if (found == true)
				P_chr_SNP_homologA{    chrID}{projectDatumID} = H_chr_SNP_alleleA{       chrID}{hapmapDatumID};;
				P_chr_SNP_homologB{    chrID}{projectDatumID} = H_chr_SNP_alleleB{       chrID}{hapmapDatumID};;
				P_chr_SNP_flipHomologs{chrID}(projectDatumID) = H_chr_SNP_hapmapEntry{   chrID}(hapmapDatumID);;
				P_chr_SNP_keep{        chrID}(projectDatumID) = 1;
				start = projectDatumID;
			else
				P_chr_SNP_keep{        chrID}(projectDatumID) = 0;
				start = 1;
			end;
		end;

		% Clean up data for chromosome.
		P_chr_SNP_data_ratios{   chrID}(P_chr_SNP_keep{chrID} == 0) = [];
		P_chr_SNP_data_positions{chrID}(P_chr_SNP_keep{chrID} == 0) = [];
		P_chr_baseCall{          chrID}(P_chr_SNP_keep{chrID} == 0) = [];
		P_chr_SNP_homologA{      chrID}(P_chr_SNP_keep{chrID} == 0) = [];
		P_chr_SNP_homologB{      chrID}(P_chr_SNP_keep{chrID} == 0) = [];
		P_chr_SNP_flipHomologs{  chrID}(P_chr_SNP_keep{chrID} == 0) = [];
		P_chr_count{             chrID}(P_chr_SNP_keep{chrID} == 0) = [];
		P_chr_SNP_keep{          chrID}(P_chr_SNP_keep{chrID} == 0) = [];
	end;
end;


%%============================================================================================================
% Save processed data file.
%-------------------------------------------------------------------------------------------------------------
fprintf(['process_2dataset_hapmap_allelicRatios.m: Save processed data.\n']);
save([project1dir 'SNP_' SNP_verString '.all3.mat'],'C_chr_SNP_data_positions','C_chr_SNP_data_ratios','C_chr_count','C_chr_baseCall','C_chr_SNP_homologA','C_chr_SNP_homologB','C_chr_SNP_flipHomologs', ...
                                                    'P_chr_SNP_data_positions','P_chr_SNP_data_ratios','P_chr_count','P_chr_baseCall','P_chr_SNP_homologA','P_chr_SNP_homologB','P_chr_SNP_flipHomologs');
