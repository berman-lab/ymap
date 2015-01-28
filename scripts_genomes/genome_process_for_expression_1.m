function results = genome_process_for_expression_1(user,genome)
% This script attempts to reformat a FASTA genome file into a new FASTA genome fragmented by restriction digest sites.
%     for use in expression analysis.

% log file start, for in-process analysis.
fprintf(['Genome : [[[' genome '[[[\n']);

workingDir             = ['../../users/' user '/genomes/' genome '/'];
nmer_length            = 10;

[centromeres, chr_sizes, figure_details, annotations, ploidy_default] = Load_genome_information_1(workingDir,genome);

%% Determine number of chromosomes from figure_details.
num_chrs   = 0;
chr_labels = [];
inUse      = zeros(1,length(figure_details));
allNames   = {};
for i = 1:length(figure_details)
	allNames{i}  = figure_details(i).name;
    if (figure_details(i).chr ~= 0)
		inUse(i) = 1;
        num_chrs = num_chrs+1;
        chr_names{figure_details(i).chr}  = figure_details(i).name;
        length_nts(figure_details(i).chr) = chr_sizes(figure_details(i).chr).size;
    end;
end;

%% Initialize the cell arrays needed to contain chromosome sequences.
sequences         = cell(1,length(figure_details));
rev_com_sequences = cell(1,length(figure_details));

%% Determine reference genome FASTA file in use.
reference_file = [workingDir '/reference.txt'];
reference_fid  = fopen(reference_file, 'r');
refFASTA       = fgetl(reference_fid);
fclose(reference_fid);
FastaName      = strrep(refFASTA,'.fasta','');

if (exist([workingDir FastaName '.expression.fasta'],'file') == 0)
    %% ====================================================================
    % Parse sequences from FASTQ file, or load MAT files if already generated for this genome.
    % ---------------------------------------------------------------------
    if (exist([workingDir FastaName 'seq.mat'],'file') == 0) && (exist([workingDir FastaName '.seqRevCom.mat'],'file') == 0)
		fprintf('\nParsing sequences from FASTQ file.');
		SequenceData = fastaread([workingDir refFASTA ],'Blockread',[1,inf],'TrimHeaders', true);
		% 'Blockread'   : read in specific entries from the FASTA file.
		% 'TrimHeaders' : trims the header lines for each FASTA entry to the first white-space delimited string.

		fprintf(['\n\tLoading FASTA file: \"' refFASTA '\"\n\n']);
		for i = 1:length(SequenceData)
			for j = 1:length(chr_names)
				if (strcmp(chr_names{j},SequenceData(i).Header) == 1)
					length_nts(j)        = length(SequenceData(i).Sequence);
					sequences{j}         = SequenceData(i).Sequence;
					rev_com_sequences{j} = rev_com(sequences{j});
				end;
			end;
		end;
		save([workingDir FastaName '.seq.mat']      ,'sequences');
		save([workingDir FastaName '.seqRevCom.mat'],'rev_com_sequences');
	else
		fprintf('\nLoading pre-determined sequences.');
		load([workingDir FastaName '.seq.mat']);
		load([workingDir FastaName '.seqRevCom.mat']);
	end;


%	%% ====================================================================
%	% Load description of "chromosome_features.txt" file from "expression.txt" file.
%	% ---------------------------------------------------------------------
%	expression_file = [workingDir '/expression.txt'];
%	expression_fid  = fopen(expression_file, 'r');
%	header_count    = str2double(fgetl(expression_fid));
%	col_chrName     = str2double(fgetl(expression_fid));
%	col_startBP     = str2double(fgetl(expression_fid));
%	col_endBP       = str2double(fgetl(expression_fid));
%	fclose(expression_fid);


	%% ====================================================================
	% Load lines from "chromosome_features_2.txt" file and parse into new FASTA entries.
	% ---------------------------------------------------------------------
	NewSequenceData   = [];
	fragment          = 1;
	fprintf('\nFragmenting genome into expression units.');
	features_file     = [workingDir 'chromosome_features_2.txt']
	features_fid      = fopen(features_file, 'r');
%	for i = 1:header_count
%		discard_line  = fgetl(features_fid);
%	end;
	discard_single_header_line = fgetl(features_fid);

	for i = 1:length(chr_names)
		test = chr_names{i}
	end;

	while not (feof(features_fid))
		loaded_line   = fgetl(features_fid);
		parts_of_line = regexp(loaded_line,'\t','split');
		chrName       = parts_of_line{1};
		startBP       = str2double(parts_of_line{2});
		endBP         = str2double(parts_of_line{3});

		% Swap coordinates so that start is lowest coordinate.
		if (startBP > endBP)
			temp      = endBP;
			endBP     = startBP;
			startBP   = temp;
		end;

		% Identify the chrID associated with the chrName.
		chrID = 0;
		for chr = 1:length(chr_names)
			if (length(chr_names{chr}) > 0)
				if strcmp(chr_names{chr},chrName)
					test  = chr_names{chr};
					chrID = chr;
				end;
			end;
		end;

		if (chrID ~= 0)
			header_string                      = ['>' genome '.chr' num2str(chrID) ' (' num2str(startBP) '..' num2str(endBP) ')'];
			sequence_string                    = sequences{chrID}(startBP:endBP);
			NewSequenceData(fragment).Header   = header_string;
			NewSequenceData(fragment).Sequence = sequence_string;
			fragment                           = fragment+1;
		end;
	end;
	fclose(features_fid);

	% Annotate fragment headers with fragment length.
	fprintf('\nAdding fragment lengths to fragment headers.');
	for fragment = 1:(length(NewSequenceData))
		header_length = length(NewSequenceData(fragment).Sequence);
		NewSequenceData(fragment).Header = [NewSequenceData(fragment).Header ' (' num2str(length(NewSequenceData(fragment).Sequence)) 'bp)'];
	end;

	% Annotate fragments headers with '[*]' so downstream scripts will know entry is valid.
	for fragment = 1:(length(NewSequenceData))
		NewSequenceData(fragment).Header = [NewSequenceData(fragment).Header ' [*]'];
	end;

	outputFile = [workingDir FastaName '.expression.fasta'];

	% Clear pre-existing processed genome file if it already exists.
	if (exist(outputFile,'file') == 0)
		fprintf('\nPre-existing expression processed genome file not found.');
	else
		fprintf('\nDeleting pre-existing expression processed genome file.');
		delete(outputFile);
	end;

	% Write out the fragmented chromosome as a genome FASTA file.
	fprintf('\nWriting out fragmented chromosomes structure as new FASTA file.\n\n');
	fastawrite(outputFile, NewSequenceData);
else
	outputFile = [workingDir FastaName '.expression.fasta'];
	fprintf(['\nReference genome was already digested.\nDigested genome is found at :\n\t' outputFile '\n\n']);
end;

end
