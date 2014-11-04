%% ========================================================================
% Reads in Tab-delimited text format output from Jason Funt's script
% processing SAM data for potential SNPs.
%
% Heterozygous SNPs are identified by having an allelic fraction between
% 0.25 and 0.75.
%
% This isolates het SNPs in any strain and has been used to identify LOHs in
% later strains in an isolate series.
%
% This version attempts to introduce speed improvements:
%    * Preprocessing only determines total SNP count, cleans up empty cells
%    later.
%    * This version attempts to deal with issues surrounding very large datasets
%    and file system limitations.
%==========================================================================
function [] = preprocess_SNP_data_files_1(projectName,SNP_verString,projectDir)
if (exist([projectDir 'SNP_' SNP_verString '.het.mat'],'file') == 0) && ...
	(exist([projectDir 'SNP_' SNP_verString '.hom1.mat'],'file') == 0) && ...
		(exist([projectDir 'SNP_' SNP_verString '.hom2.mat'],'file') == 0)    % Control variables for SNP and CGH probe lists.

	datafile = [projectDir 'putative_SNPs_' SNP_verString '.txt'];

	%% This section determines how many lines are in the SNP data file.   The number of lines is used to assign
	% large, memory-intensive, arrays for storing the actual data about those SNPs.   Multiple arrays were
	% needed to avoid the problem of running out of memory from placing everything into one data structure on
	% desktop system.   This requirement has not been tested in the server system.
	j     = 0;
	k     = 0;
	fprintf(['\tPreprocessing raw SNP data: ' datafile]);
	[status, result] = system(['wc -l < ' datafile]);
	tot = str2num(result)-1;
	fprintf(['\n\t\t' num2str(tot) ' lines.']);

	%% This section runs through the data file and stores data about the SNPs into data structures for later
	% analysis.   A separate data structure is made for SNPs identified as heterozygous and those identified
	% as homozygous.
	j     = 0;
	k     = 0;
	data  = fopen(datafile);
	null  = fgetl(data);
	fprintf(['\n\tLoading raw SNP data: ' datafile '\n\t\t[']);
	dataset_het           = [];
	dataset_hom1          = [];
	dataset_hom2          = [];
	dataset_het{tot}.chr  = 'null';
	dataset_hom1{tot}.chr = 'null';
	dataset_hom2{tot}.chr = 'null';
	het_set               = zeros(1,tot);
	hom_set1              = zeros(1,tot);
	hom_set2              = zeros(1,tot);
	counts                = zeros(1,tot);
	hets                  = 0;
	homs1                 = 0;
	homs2                 = 0;

	while ~feof(data)
		j = j+1;
		% load line from data file and tokenize into interesting collumns.
		line_of_text = fgetl(data);
		[chr      ,remain] = strtok(line_of_text);
		[position ,remain] = strtok(remain);
		[reference,remain] = strtok(remain);
		[A        ,remain] = strtok(remain);
		[T        ,remain] = strtok(remain);
		[G        ,remain] = strtok(remain);
		[C        ,remain] = strtok(remain);
		A = str2double(A);
		T = str2double(T);
		G = str2double(G);
		C = str2double(C);
		counts(j) = A+T+G+C;

		% determine the highest and second highest base read count.
		list = [A T G C];
		v1 = max(list);
		list(find(list == v1,1)) = [];
		v2 = max(list);
        
		% randomly choose 1' vs 2' to assign potential SNPs randomly to two halves.
		list_ = [A T G C];
		if (rand() < 0.5)
			num1	= v1;
			num2	= v2;
		else
			num1	= v2;
			num2	= v1;
		end;

		% determine which bases are the 1' and 2'.
		a1 = find(list_ == v1);
		if (length(a1) > 1)
			% if more than one base has the number 'v1', choose randomly.
			a1 = a1(ceil(rand(1)*length(a1)));
		end;
		a2 = find(list_ == v2);
		if (length(a2) > 1)
			% if more than one base has the number 'v2', choose randomly.
			a2 = a2(ceil(rand(1)*length(a2)));
		end;
		% determine letter associated with base.
		if (a1 == 1)		allele1 = 'A';
		elseif(a1 == 2)		allele1 = 'T';
		elseif(a1 == 3)		allele1 = 'G';
		elseif(a1 == 4)		allele1 = 'C';
		end;
		if (a2 == 1)		allele2 = 'A';
		elseif(a2 == 2)		allele2 = 'T';
		elseif(a2 == 3)		allele2 = 'G';
		elseif(a2 == 4)		allele2 = 'C';
		end;
		% calculate the total count of reads.
		total = A+T+G+C;

		% calculate the allelic fraction; here the random 1' vs. 2' choice has impact.
		allelicFraction = num1/(total);
		%angle = atan2(num1,num2);

		% Build data structures containing heterozygous and homozygous data, defined by allelic ratios.
		% [0.00,0.25) : Homozygous aa.
		% [0.25,0.75] : Heterozygous ab.
		% (0.75,1.00] : Homozygous bb.
		if (allelicFraction > 0.25) && (allelicFraction < 0.75)
			% The alleleic ratios closer to 1:1 are considered heterozygous, at least initially.
			hets = hets+1;
			dataset_het{hets}.chr      = chr;
			dataset_het{hets}.position = position;
			dataset_het{hets}.a        = allele1;
			dataset_het{hets}.b        = allele2;
			dataset_het{hets}.a_count  = num1;
			dataset_het{hets}.b_count  = num2;
			dataset_het{hets}.total    = total;
			het_set(hets)              = 1;
		else
			% The allelic ratios closer to 1:0 and 0:1 are considered homozygous, at least initially.   These are broken into two
			% datasets to avoid file system and memory size limits on single structures.
			if (allelicFraction < 0.5)
				homs1                        = homs1+1;
				dataset_hom1{homs1}.chr      = chr;
				dataset_hom1{homs1}.position = position;
				dataset_hom1{homs1}.a        = allele1;
				dataset_hom1{homs1}.b        = allele2;
				dataset_hom1{homs1}.a_count  = num1;
				dataset_hom1{homs1}.b_count  = num2;
				dataset_hom1{homs1}.total    = total;
				hom_set1(homs1)              = 1;
			else
				homs2                        = homs2+1;
				dataset_hom2{homs2}.chr      = chr;
				dataset_hom2{homs2}.position = position;
				dataset_hom2{homs2}.a        = allele1;
				dataset_hom2{homs2}.b        = allele2;
				dataset_hom2{homs2}.a_count  = num1;
				dataset_hom2{homs2}.b_count  = num2;
				dataset_hom2{homs2}.total    = total;
				hom_set2(homs2)              = 1;
			end;
		end;

		% update text output to let user know script is progressing.
		if (mod(j,100000) == 0)
			fprintf('.');
		end;
        if (mod(j,5800000) == 0)
			k = k+1;
			fprintf([num2str(6000*k) 'K\n\t\t ']);
		end;
	end;
	fprintf(']');
	fclose(data);

	% Cleanup of data structures: removes empty end of data structures.
	fprintf('\n\tCleaning het/hom SNP datasets. ');
	fprintf('[1');
	dataset_het(hets+1:end) = [];
	fprintf(',2');
	dataset_hom1(homs1+1:end) = [];
	fprintf(',3');
	dataset_hom2(homs2+1:end) = [];
	fprintf(']');

	counts(counts == 0) = [];
	ave_read_count = median(counts);

	%% Save files to scratch space.
	fprintf('\n\tSaving het SNP dataset.');
	save([projectDir 'SNP_' SNP_verString '.het.mat'],'dataset_het','ave_read_count','-v7.3');
	fprintf('\n\tSaving hom SNP dataset. [1/2]');
	save([projectDir 'SNP_' SNP_verString '.hom1.mat'],'dataset_hom1','-v7.3');
	fprintf('\n\tSaving hom SNP dataset. [2/2]');
	save([projectDir 'SNP_' SNP_verString '.hom2.mat'],'dataset_hom2','-v7.3');

	fprintf('\n');
else
	fprintf(['\t' projectName ' already processed.\n']);
end;
end
