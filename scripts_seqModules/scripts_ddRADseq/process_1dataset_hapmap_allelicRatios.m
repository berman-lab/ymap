function [] = process_1dataset_hapmap_allelicRatios(project1dir, project2dir, chr_size, chr_name, chr_in_use, SNP_verString);


%%============================================================================================================
% Define file names for processing.
%-------------------------------------------------------------------------------------------------------------
fprintf('Define file names for processing.\n');
C_datafile = [project1dir 'putative_SNPs_v4.txt'];
P_datafile = [project2dir 'SNPdata_parent.txt'];


%%============================================================================================================
% Preallocate data vectors the length of each chromosome.
%-------------------------------------------------------------------------------------------------------------
fprintf('Preallocate data vectors for each chromosome.\n');
C_chr_SNP_data_positions = cell(length(chr_size),1);
C_chr_SNP_data_ratios    = cell(length(chr_size),1);
C_chr_baseCall           = cell(length(chr_size),1);
C_chr_count              = cell(length(chr_size),1);
C_chr_SNP_homologA       = cell(length(chr_size),1);
C_chr_SNP_homologB       = cell(length(chr_size),1);
C_chr_SNP_flipHomologs   = cell(length(chr_size),1);
C_chr_SNP_keep           = cell(length(chr_size),1);

P_chr_SNP_data_positions = cell(length(chr_size),1);
P_chr_SNP_alleleA        = cell(length(chr_size),1);
P_chr_SNP_alleleB        = cell(length(chr_size),1);
P_chr_SNP_needToFlip     = cell(length(chr_size),1);

for chrID = 1:length(chr_size)
	if (chr_in_use(chrID) == 1)
		C_chr_SNP_data_positions{chrID} = zeros(chr_size(chrID),1);
		C_chr_SNP_data_ratios{   chrID} = zeros(chr_size(chrID),1);
		C_chr_count{             chrID} = zeros(chr_size(chrID),1);
		C_chr_baseCall{          chrID} = cell( chr_size(chrID),1);
		C_chr_SNP_homologA{      chrID} = cell( chr_size(chrID),1);
		C_chr_SNP_homologB{      chrID} = cell( chr_size(chrID),1);
		C_chr_SNP_flipHomologs{  chrID} = zeros(chr_size(chrID),1);
		C_chr_SNP_keep{          chrID} = ones(chr_size(chrID),1);
		C_chr_lines_analyzed(    chrID) = 0;
		P_chr_SNP_data_positions{chrID} = zeros(chr_size(chrID),1);
		P_chr_SNP_alleleA{       chrID} = cell( chr_size(chrID),1);
		P_chr_SNP_alleleB{       chrID} = cell( chr_size(chrID),1);
		P_chr_SNP_needToFlip{    chrID} = zeros(chr_size(chrID),1);
		P_chr_lines_analyzed(    chrID) = 0;
	end;
end;


%%============================================================================================================
% Load project information.
%-------------------------------------------------------------------------------------------------------------
fprintf('Load project information.\n');
C_data      = fopen(C_datafile, 'r');
allele_list = ['A' 'T' 'G' 'C'];
old_chr     = 0;
while not (feof(C_data))
	C_dataLine = fgetl(C_data);
	if (length(C_dataLine) > 0)
		% process the loaded line into data channels.
		C_SNP_chr_name   = sscanf(C_dataLine, '%s',1);
		C_SNP_coordinate = sscanf(C_dataLine, '%s',2);   for i = 1:size(sscanf(C_dataLine,'%s',1),2);   C_SNP_coordinate(1) = [];   end;
		C_SNP_reference  = sscanf(C_dataLine, '%s',3);   for i = 1:size(sscanf(C_dataLine,'%s',2),2);   C_SNP_reference(1)  = [];   end;
		C_SNP_countA     = sscanf(C_dataLine, '%s',4);   for i = 1:size(sscanf(C_dataLine,'%s',3),2);   C_SNP_countA(1)     = [];   end;
		C_SNP_countT     = sscanf(C_dataLine, '%s',5);   for i = 1:size(sscanf(C_dataLine,'%s',4),2);   C_SNP_countT(1)     = [];   end;
		C_SNP_countG     = sscanf(C_dataLine, '%s',6);   for i = 1:size(sscanf(C_dataLine,'%s',5),2);   C_SNP_countG(1)     = [];   end;
		C_SNP_countC     = sscanf(C_dataLine, '%s',7);   for i = 1:size(sscanf(C_dataLine,'%s',6),2);   C_SNP_countC(1)     = [];   end;
		C_chr_num        = find(strcmp(C_SNP_chr_name, chr_name));
		if (length(C_chr_num) > 0)
			if (C_chr_num ~= old_chr)
				fprintf(['\tchr = ' num2str(C_chr_num) '\n']);
			end;
			C_SNP_countA       = str2num(C_SNP_countA);
			C_SNP_countT       = str2num(C_SNP_countT);
			C_SNP_countG       = str2num(C_SNP_countG);
			C_SNP_countC       = str2num(C_SNP_countC);
			C_count_vector1    = [C_SNP_countA C_SNP_countT C_SNP_countG C_SNP_countC];
			C_chr_read_max1    = max(C_count_vector1);
			C_SNP_coordinate   = str2num(C_SNP_coordinate);
			C_chr_lines_analyzed(C_chr_num) = C_chr_lines_analyzed(C_chr_num)+1;
			C_chr_SNP_data_positions{C_chr_num}(C_chr_lines_analyzed(C_chr_num)) = C_SNP_coordinate;
			C_chr_SNP_data_ratios{   C_chr_num}(C_chr_lines_analyzed(C_chr_num)) = C_chr_read_max1/sum(C_count_vector1);
			C_chr_count{             C_chr_num}(C_chr_lines_analyzed(C_chr_num)) = sum(C_count_vector1);
			allele_call_id = find(C_count_vector1==max(C_count_vector1));
			if (length(allele_call_id) > 1)
				C_chr_read_id = 'N';
			else
				C_chr_read_id = allele_list(allele_call_id);
			end;
			C_chr_baseCall{          C_chr_num}{C_chr_lines_analyzed(C_chr_num)} = C_chr_read_id;
			old_chr = C_chr_num;
		else
			old_chr = 0;
		end;
	end;
end;
fclose(C_data);


%%============================================================================================================
% Load hapmap information.
%-------------------------------------------------------------------------------------------------------------
fprintf('Load hapmap information.\n');
P_data      = fopen(P_datafile, 'r');
old_chr     = 0;
while not (feof(P_data))
	P_dataLine = fgetl(P_data);
	if (length(P_dataLine) > 0)
		% process the loaded line into data channels.
		P_SNP_chr_name   = sscanf(P_dataLine, '%s',1);
		P_SNP_coordinate = sscanf(P_dataLine, '%s',2);   for i = 1:size(sscanf(P_dataLine,'%s',1),2);   P_SNP_coordinate(1) = [];   end;
		P_SNP_alleleA    = sscanf(P_dataLine, '%s',3);   for i = 1:size(sscanf(P_dataLine,'%s',2),2);   P_SNP_alleleA(1)    = [];   end;
		P_SNP_alleleB    = sscanf(P_dataLine, '%s',4);   for i = 1:size(sscanf(P_dataLine,'%s',3),2);   P_SNP_alleleB(1)    = [];   end;
		P_SNP_needToFlip = sscanf(P_dataLine, '%s',5);   for i = 1:size(sscanf(P_dataLine,'%s',4),2);   P_SNP_needToFlip(1) = [];   end;
		P_chr_num        = find(strcmp(P_SNP_chr_name, chr_name));
		if (length(P_chr_num) > 0)
			if (P_chr_num ~= old_chr)
				fprintf(['\tchr = ' num2str(P_chr_num) '\n']);
			end;
			P_SNP_coordinate                                                     = str2num(P_SNP_coordinate);
			P_SNP_needToFlip                                                     = str2num(P_SNP_needToFlip);
			P_chr_lines_analyzed(P_chr_num)                                      = P_chr_lines_analyzed(P_chr_num)+1;
			P_chr_SNP_data_positions{P_chr_num}(P_chr_lines_analyzed(P_chr_num)) = P_SNP_coordinate;
			P_chr_SNP_alleleA{       P_chr_num}{P_chr_lines_analyzed(P_chr_num)} = P_SNP_alleleA;
			P_chr_SNP_alleleB{       P_chr_num}{P_chr_lines_analyzed(P_chr_num)} = P_SNP_alleleB;
			P_chr_SNP_needToFlip{    P_chr_num}(P_chr_lines_analyzed(P_chr_num)) = P_SNP_needToFlip;
			old_chr = P_chr_num;
		else
			old_chr = 0;
		end;
	end;
end;
fclose(P_data);


%%============================================================================================================
% Clean up data vectors.
%-------------------------------------------------------------------------------------------------------------
fprintf('Clean up data vectors.\n');
for chrID = 1:length(chr_size)
	if (chr_in_use(chrID) == 1)
		fprintf(['\tchr = ' num2str(chrID) '\n']);
		C_chr_SNP_data_ratios{   chrID}(C_chr_SNP_data_positions{chrID} == 0)  = [];
		C_chr_count{             chrID}(C_chr_SNP_data_positions{chrID} == 0)  = [];
		C_chr_baseCall{          chrID}(C_chr_SNP_data_positions{chrID} == 0)  = []; 
		C_chr_SNP_homologA{      chrID}(C_chr_SNP_data_positions{chrID} == 0)  = [];
		C_chr_SNP_homologB{      chrID}(C_chr_SNP_data_positions{chrID} == 0)  = [];
		C_chr_SNP_flipHomologs{  chrID}(C_chr_SNP_data_positions{chrID} == 0)  = [];
		C_chr_SNP_keep{          chrID}(C_chr_SNP_data_positions{chrID} == 0)  = [];
		C_chr_SNP_data_positions{chrID}(C_chr_SNP_data_positions{chrID} == 0)  = [];
		C_chr_SNP_data_ratios{   chrID}(C_chr_count{             chrID} <= 20) = [];
		C_chr_SNP_data_positions{chrID}(C_chr_count{             chrID} <= 20) = [];
		C_chr_baseCall{          chrID}(C_chr_count{             chrID} <= 20) = [];
		C_chr_SNP_homologA{      chrID}(C_chr_count{             chrID} <= 20) = []; 
		C_chr_SNP_homologB{      chrID}(C_chr_count{             chrID} <= 20) = []; 
		C_chr_SNP_flipHomologs{  chrID}(C_chr_count{             chrID} <= 20) = []; 
		C_chr_SNP_keep{          chrID}(C_chr_count{             chrID} <= 20) = []; 
		C_chr_count{             chrID}(C_chr_count{             chrID} <= 20) = [];

		P_chr_SNP_alleleA{       chrID}(P_chr_SNP_data_positions{chrID} == 0)  = [];
		P_chr_SNP_alleleB{       chrID}(P_chr_SNP_data_positions{chrID} == 0)  = [];
		P_chr_SNP_needToFlip{    chrID}(P_chr_SNP_data_positions{chrID} == 0)  = [];
		P_chr_SNP_data_positions{chrID}(P_chr_SNP_data_positions{chrID} == 0)  = [];
	end;
end;


%%============================================================================================================
% Determine hapmap information for dataset values.
%-------------------------------------------------------------------------------------------------------------
fprintf('Determine hapmap information for dataset values.\n');
start = 1;
for chrID = 1:length(chr_size)
	if (chr_in_use(chrID) == 1)
		fprintf(['\tchr = ' num2str(chrID) '\n']);
		for projectDatumID = 1:length(C_chr_SNP_data_positions{chrID})
			pos   = C_chr_SNP_data_positions{chrID}(projectDatumID);
			found = false;
			for hapmapDatumID = start:length(P_chr_SNP_data_positions{chrID})
				hapmap_pos = P_chr_SNP_data_positions{chrID}(hapmapDatumID);
				if (pos == hapmap_pos)
					found = true;
					break;
				end;
			end;
			if (found == true)
				fprintf(['\tchr' num2str(chrID) ':' num2str(pos) ' +\n']);
				C_chr_SNP_homologA{    chrID}{projectDatumID} = P_chr_SNP_alleleA{       chrID}{hapmapDatumID};;
				C_chr_SNP_homologB{    chrID}{projectDatumID} = P_chr_SNP_alleleB{       chrID}{hapmapDatumID};;
				C_chr_SNP_flipHomologs{chrID}(projectDatumID) = P_chr_SNP_needToFlip{    chrID}(hapmapDatumID);;
				C_chr_SNP_keep{        chrID}(projectDatumID) = 1;
				start = projectDatumID;
			else
				% fprintf(['\tchr' num2str(chrID) ':' num2str(pos) '\n']);
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
% Save processed data file.
%-------------------------------------------------------------------------------------------------------------
fprintf('Save output data file.\n');
save([project1dir 'SNP_' SNP_verString '.all2.mat'],'C_chr_SNP_data_positions','C_chr_SNP_data_ratios','C_chr_count','C_chr_baseCall','C_chr_SNP_homologA','C_chr_SNP_homologB','C_chr_SNP_flipHomologs');
