function results = genome_process_for_RADseq_1_BamHI_BclI(user,genome)
% This script attempts to reformat a FASTA genome file into a new FASTA genome fragmented by restriction digest sites.
%     for use in ddRADseq analysis.

% log file start, for in-process analysis.
fprintf(['Genome : [[[' genome '[[[\n']);

workingDir             = ['/heap/hapmap/bermanlab/users/' user '/genomes/' genome '/'];
figureDir              = '~/';
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
testVar = inUse

%% Initialize the cell arrays needed to contain chromosome sequences.
sequences         = cell(1,length(figure_details));
rev_com_sequences = cell(1,length(figure_details));

%% Determine reference genome FASTA file in use.
reference_file = [workingDir '/reference.txt'];
refernce_fid   = fopen(reference_file, 'r');
refFASTA       = fgetl(refernce_fid);
fclose(refernce_fid);
FastaName      = strrep(refFASTA,'.fasta','');

if (exist([workingDir FastaName '.BamHI_BclI.fasta'],'file') == 0)
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

    %% ====================================================================
    % Fragment genome by restriction sites into new FASTA file.
    % ---------------------------------------------------------------------
    % restriction digest site 1: BamHI [G:GATCC]
    % restriction digest site 2: BclI  [T:GATCA]

    NewSequenceData = [];
    fragment = 1;
    fprintf('\nFragmenting genome by BamHI & BclI');
    for chr = 1:length(figure_details)
		if (inUse(chr) == 1)
			start_coordinate = 1;
			fprintf(['\n\tFragmenting : ' allNames{chr} ]);
			for bp = start_coordinate:length(sequences{chr})
				if (bp < length(sequences{chr})-5)
					% fragment chromosome by BamHI [G:GATCC].
					test_sequence = sequences{chr}(bp:(bp+5));
					if (strcmp(upper(test_sequence),'GGATCC') == 1)
						header_string                      = ['>' genome '.chr' num2str(chr) ' (' num2str(start_coordinate) '..' num2str(bp) ')'];
						sequence_string                    = sequences{chr}(start_coordinate:bp);
						NewSequenceData(fragment).Header   = header_string;
						NewSequenceData(fragment).Sequence = sequence_string;
						start_coordinate                   = bp+1;
						fragment                           = fragment+1;
					end;
					% fragment chromosome by BclI [T:GATCA].
					test_sequence = sequences{chr}(bp:(bp+5));
					if (strcmp(upper(test_sequence),'TGATCA') == 1)
						header_string                      = ['>' genome '.chr' num2str(chr) ' (' num2str(start_coordinate) '..' num2str(bp-1) ')'];
						sequence_string                    = sequences{chr}(start_coordinate:bp);
						NewSequenceData(fragment).Header   = header_string;
						NewSequenceData(fragment).Sequence = sequence_string;
						start_coordinate                   = bp+1;
						fragment                           = fragment+1;
					end;
				end;
				if (bp == length(sequences{chr}))
					header_string                      = ['>' genome '.chr' num2str(chr) ' (' num2str(start_coordinate) '..' num2str(length(sequences{chr})) ')'];
					sequence_string                    = sequences{chr}(start_coordinate:end);
					NewSequenceData(fragment).Header   = header_string;
					NewSequenceData(fragment).Sequence = sequence_string;
					fragment                           = fragment+1;
				end;
			end;
		else
			fprintf(['\n\tSkipping : ' allNames{chr}]);
		end;
	end;

	% Annotate fragment headers with fragment length.
	fprintf('\nAdding fragment lengths to fragment headers.');
	for fragment = 1:(length(NewSequenceData))
		header_length = length(NewSequenceData(fragment).Sequence);
		NewSequenceData(fragment).Header = [NewSequenceData(fragment).Header ' (' num2str(length(NewSequenceData(fragment).Sequence)) 'bp)'];
	end;

	% Annotate fragments headers with '[*]' if fragment is bounded by one restriction site at one end and the second restrion site at the other end.
	fprintf('\nChecking which fragments will be amplified during ddRADseq amplification.');
	fragment_test = zeros(1,length(NewSequenceData));
	fragment_test(end) = 0;
	for fragment = 1:(length(NewSequenceData)-1)
		if (length(NewSequenceData(fragment).Sequence) > 5) && (length(NewSequenceData(fragment+1).Sequence) > 5)
			% [G:GATCC]..[T:GATCA]
			test1 = upper(NewSequenceData(fragment  ).Sequence(1:5));
			test2 = upper(NewSequenceData(fragment  ).Sequence(end));
			test3 = upper(NewSequenceData(fragment+1).Sequence(1:5));
			% [T:GATCA]..[G:GATCC]
			test4 = upper(NewSequenceData(fragment  ).Sequence(1:5));
			test5 = upper(NewSequenceData(fragment  ).Sequence(end));
			test6 = upper(NewSequenceData(fragment+1).Sequence(1:4));

			if ((strcmp(test1,'GATCC') == 1) && (strcmp(test2,'T') == 1) && (strcmp(test3,'GATCA') == 1)) || ...
			   ((strcmp(test4,'GATCA') == 1) && (strcmp(test5,'G') == 1) && (strcmp(test6,'GATCC') == 1))
				NewSequenceData(fragment).Header = [NewSequenceData(fragment).Header ' [*]'];
				fragment_test(fragment) = 1;
			else
				% testing all fragments per MZAnderson protocol.
				NewSequenceData(fragment).Header = [NewSequenceData(fragment).Header ' [*]'];
				fragment_test(fragment) = 1;
			end;
		end;
	end;

	fprintf('\nClearing unusable restriction fragments from digested genome.');
	NewSequenceData(fragment_test == 0) = [];

	outputFile = [workingDir FastaName '.BamHI_BclI.fasta'];

	% Clear pre-existing processed genome file if it already exists.
	if (exist(outputFile,'file') == 0)
		fprintf('\nPre-existing ddRADseq processed genome file not found.');
	else
		fprintf('\nDeleting pre-existing ddRADseq processed genome file.');
		delete(outputFile);
	end;

	% Write out the fragmented chromosome as a genome FASTA file.
	fprintf('\nWriting out fragmented chromosomes structure as new FASTA file.\n\n');
	fastawrite(outputFile, NewSequenceData);
else
	outputFile = [workingDir FastaName '.BamHI_BclI.fasta'];
	fprintf(['\nReference genome was already digested.\nDigested genome is found at :\n\t' outputFile '\n\n']);
end;

end
