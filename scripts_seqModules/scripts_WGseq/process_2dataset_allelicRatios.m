function [] = process_2dataset_allelicRatios(project1dir, project2dir, chr_size, chr_name, chr_in_use, SNP_verString);


%%============================================================================================================
% Define files for processing.
%-------------------------------------------------------------------------------------------------------------
C_datafile = [project1dir 'putative_SNPs_v4.txt'];
P_datafile = [project2dir 'putative_SNPs_v4.txt'];


%%============================================================================================================
% Preallocate data vectors the length of each chromosome.
%-------------------------------------------------------------------------------------------------------------
fprintf(['process_2dataset_allelicRatios.m: Preallocate data vectors.\n']);
C_chr_SNP_data_positions = cell(length(chr_size),1);
C_chr_SNP_data_ratios    = cell(length(chr_size),1);
C_chr_count              = cell(length(chr_size),1);
P_chr_SNP_data_positions = cell(length(chr_size),1);
P_chr_SNP_data_ratios    = cell(length(chr_size),1);
P_chr_count              = cell(length(chr_size),1);
for chrID = 1:length(chr_size)
	if (chr_in_use(chrID) == 1)
		C_chr_SNP_data_positions{chrID} = zeros(chr_size(chrID),1);
		C_chr_SNP_data_ratios{   chrID} = zeros(chr_size(chrID),1);
		C_chr_count{             chrID} = zeros(chr_size(chrID),1);
		C_chr_lines_analyzed(    chrID) = 0;
		P_chr_SNP_data_positions{chrID} = zeros(chr_size(chrID),1);
	 	P_chr_SNP_data_ratios{   chrID} = zeros(chr_size(chrID),1);
		P_chr_count{             chrID} = zeros(chr_size(chrID),1);
		P_chr_lines_analyzed(    chrID) = 0;
	end;
end;


%%============================================================================================================
% Process project 1 dataset.
%-------------------------------------------------------------------------------------------------------------
fprintf(['process_2dataset_allelicRatios.m: Process project 1 dataset.\n']);
C_data     = fopen(C_datafile, 'r');
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
			C_SNP_countA                                                         = str2num(C_SNP_countA);
			C_SNP_countT                                                         = str2num(C_SNP_countT);
			C_SNP_countG                                                         = str2num(C_SNP_countG);
			C_SNP_countC                                                         = str2num(C_SNP_countC);
			C_count_vector                                                       = [C_SNP_countA C_SNP_countT C_SNP_countG C_SNP_countC];
			C_SNP_coordinate                                                     = str2num(C_SNP_coordinate);
			C_chr_lines_analyzed(C_chr_num)                                      = C_chr_lines_analyzed(C_chr_num)+1;
			C_chr_SNP_data_positions{C_chr_num}(C_chr_lines_analyzed(C_chr_num)) = C_SNP_coordinate;
			C_chr_SNP_data_ratios   {C_chr_num}(C_chr_lines_analyzed(C_chr_num)) = max(C_count_vector)/sum(C_count_vector);
			C_chr_count             {C_chr_num}(C_chr_lines_analyzed(C_chr_num)) = sum(C_count_vector);
		end;
	end;
end;
fclose(C_data);


%%============================================================================================================
% Process project 2 dataset.
%-------------------------------------------------------------------------------------------------------------
fprintf(['process_2dataset_allelicRatios.m: Process project 2 dataset.\n']);
P_data     = fopen(P_datafile, 'r');
while not (feof(P_data))
	P_dataLine = fgetl(P_data);
	if (length(P_dataLine) > 0)
		% process the loaded line into data channels.
		P_SNP_chr_name   = sscanf(P_dataLine, '%s',1);
		P_SNP_coordinate = sscanf(P_dataLine, '%s',2);   for i = 1:size(sscanf(P_dataLine,'%s',1),2);   P_SNP_coordinate(1) = [];   end;
		P_SNP_reference  = sscanf(P_dataLine, '%s',3);   for i = 1:size(sscanf(P_dataLine,'%s',2),2);   P_SNP_reference(1)  = [];   end;
		P_SNP_countA     = sscanf(P_dataLine, '%s',4);   for i = 1:size(sscanf(P_dataLine,'%s',3),2);   P_SNP_countA(1)     = [];   end;
		P_SNP_countT     = sscanf(P_dataLine, '%s',5);   for i = 1:size(sscanf(P_dataLine,'%s',4),2);   P_SNP_countT(1)     = [];   end;
		P_SNP_countG     = sscanf(P_dataLine, '%s',6);   for i = 1:size(sscanf(P_dataLine,'%s',5),2);   P_SNP_countG(1)     = [];   end;
		P_SNP_countC     = sscanf(P_dataLine, '%s',7);   for i = 1:size(sscanf(P_dataLine,'%s',6),2);   P_SNP_countC(1)     = [];   end;
		P_chr_num        = find(strcmp(P_SNP_chr_name, chr_name));
		if (length(P_chr_num) > 0)
			P_SNP_countA                                                         = str2num(P_SNP_countA);
			P_SNP_countT                                                         = str2num(P_SNP_countT);
			P_SNP_countG                                                         = str2num(P_SNP_countG);
			P_SNP_countC                                                         = str2num(P_SNP_countC);
			P_count_vector                                                       = [P_SNP_countA P_SNP_countT P_SNP_countG P_SNP_countC];
			P_SNP_coordinate                                                     = str2num(P_SNP_coordinate);
			P_chr_lines_analyzed(P_chr_num)                                      = P_chr_lines_analyzed(P_chr_num)+1;
			P_chr_SNP_data_positions{P_chr_num}(P_chr_lines_analyzed(P_chr_num)) = P_SNP_coordinate;
			P_chr_SNP_data_ratios   {P_chr_num}(P_chr_lines_analyzed(P_chr_num)) = max(P_count_vector)/sum(P_count_vector);
			P_chr_count             {P_chr_num}(P_chr_lines_analyzed(P_chr_num)) = sum(P_count_vector);
		end;
	end;
end;
fclose(P_data);


%%============================================================================================================
% Clean up data vectors.
%-------------------------------------------------------------------------------------------------------------
fprintf(['process_2dataset_allelicRatios.m: Clean up data.\n']);
for chrID = 1:length(chr_size)
	if (chr_in_use(chrID) == 1)
		C_chr_SNP_data_ratios{   chrID}(C_chr_SNP_data_positions{chrID} == 0) = [];
		C_chr_count{             chrID}(C_chr_SNP_data_positions{chrID} == 0) = [];
		C_chr_SNP_data_positions{chrID}(C_chr_SNP_data_positions{chrID} == 0) = [];
		C_chr_SNP_data_ratios{   chrID}(C_chr_count{chrID} < 20)              = [];
		C_chr_SNP_data_positions{chrID}(C_chr_count{chrID} < 20)              = [];
		C_chr_count{             chrID}(C_chr_count{chrID} < 20)              = [];

		P_chr_SNP_data_ratios{   chrID}(P_chr_SNP_data_positions{chrID} == 0) = [];
		P_chr_count{             chrID}(P_chr_SNP_data_positions{chrID} == 0) = [];
		P_chr_SNP_data_positions{chrID}(P_chr_SNP_data_positions{chrID} == 0) = [];
		P_chr_SNP_data_ratios{   chrID}(P_chr_count{chrID} < 20)              = [];
		P_chr_SNP_data_positions{chrID}(P_chr_count{chrID} < 20)              = [];
		P_chr_count{             chrID}(P_chr_count{chrID} < 20)              = [];
	end;
end;


%%============================================================================================================
% Save processed data file.
%-------------------------------------------------------------------------------------------------------------
fprintf(['process_2dataset_allelicRatios.m: Save data.\n']);
save([project1dir 'SNP_' SNP_verString '.all1.mat'],'C_chr_SNP_data_ratios','C_chr_SNP_data_positions','P_chr_SNP_data_ratios','P_chr_SNP_data_positions','C_chr_count','P_chr_count');
