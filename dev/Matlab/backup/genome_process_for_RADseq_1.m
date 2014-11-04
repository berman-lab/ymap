function results = genome_process_for_RADseq_1(user,genome)
% This script attempts to reformat a FASTA genome file into a new FASTA genome fragmented by restriction digest sites.
%     for use in ddRADseq analysis.

% log file start, for in-process analysis.
fprintf(['Genome : [[[' genome '[[[\n']);

workingDir             = ['/heap/hapmap/bermanlab/users/' user '/genomes/' genome '/'];
figureDir              = '~/';
nmer_length            = 10;

[centromeres, chr_sizes, figure_details, annotations, ploidy_default] = Load_genome_information_1(workingDir,figureDir,user,genome);

%% Determine number of chromosomes from figure_details.
num_chr    = 0;
chr_labels = [];
for i = 1:length(figure_details)
    if (figure_details(i).chr ~= 0)
        num_chr = num_chr+1;
        chr_names{figure_details(i).chr}  = figure_details(i).name;
        length_nts(figure_details(i).chr) = chr_sizes(figure_details(i).chr).size;
    end;
end;

%% Initialize the cell arrays needed to contain chromosome sequences.
sequences         = cell(1,num_chr);
rev_com_sequences = cell(1,num_chr);

%% Determine reference genome FASTA file in use.
reference_file = [workingDir '/reference.txt'];
refernce_fid   = fopen(reference_file, 'r');
refFASTA       = fgetl(refernce_fid);
fclose(refernce_fid);

if (exist([workingDir strrep(refFASTA,'fasta','MfeI_MboI.fasta')],'file') == 0)
    %% ====================================================================
    % Parse sequences from FASTQ file, or load MAT files if already generated for this genome.
    % ---------------------------------------------------------------------
    if (exist([workingDir genome '_seq.mat'],'file') == 0) && (exist([workingDir genome '_seqRevCom.mat'],'file') == 0)
		fprintf('\nParsing sequences from FASTQ file.');
		SequenceData = fastaread([workingDir refFASTA ],'Blockread',[1,num_chr],'TrimHeaders', true);
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
		save([workingDir genome '_seq.mat']      ,'sequences');
		save([workingDir genome '_seqRevCom.mat'],'rev_com_sequences');
	else
		fprintf('\nLoading pre-determined sequences.');
		load([workingDir genome '_seq.mat']);
		load([workingDir genome '_seqRevCom.mat']);
	end;

    %% ====================================================================
    % Fragment genome by restriction sites into new FASTA file.
    % ---------------------------------------------------------------------
    % restriction digest site 1: MfeI [C:AATTG]
    % restriction digest site 2: MboI [:GATC]

    NewSequenceData = [];
    fragment = 1;
    fprintf('\nFragmenting genome by MfeI & MboI');
    for chr = 1:num_chr
		start_coordinate = 1;
		fprintf(['\n\tFragmenting chr' num2str(chr)]);
		for bp = start_coordinate:length(sequences{chr})
			if (bp < length(sequences{chr})-5)
				% fragment chromosome by MfeI [C:AATTG].
				test_sequence = sequences{chr}(bp:(bp+5));
				if (strcmp(upper(test_sequence),'CAATTG') == 1)
					header_string                      = ['>' genome '.chr' num2str(chr) ' (' num2str(start_coordinate) '..' num2str(bp) ')'];
					sequence_string                    = sequences{chr}(start_coordinate:bp);
					NewSequenceData(fragment).Header   = header_string;
					NewSequenceData(fragment).Sequence = sequence_string;
					start_coordinate                   = bp+1;
					fragment                           = fragment+1;
				end;

				% fragment chromosome by MboI [:GATC].
				test_sequence = sequences{chr}(bp:(bp+3));
				if (strcmp(upper(test_sequence),'GATC') == 1)
					header_string                      = ['>' genome '.chr' num2str(chr) ' (' num2str(start_coordinate) '..' num2str(bp-1) ')'];
					sequence_string                    = sequences{chr}(start_coordinate:(bp-1));
					NewSequenceData(fragment).Header   = header_string;
					NewSequenceData(fragment).Sequence = sequence_string;
					start_coordinate                   = bp;
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
		if (length(NewSequenceData(fragment).Sequence) > 4) && (length(NewSequenceData(fragment+1).Sequence) > 4)
			% [:GATC]..[C:AATTG]
			test1 = upper(NewSequenceData(fragment  ).Sequence(1:4));
			test2 = upper(NewSequenceData(fragment  ).Sequence(end));
			test3 = upper(NewSequenceData(fragment+1).Sequence(1:5));
			% [C:AATTG]..[:GATC]
			test4 = upper(NewSequenceData(fragment  ).Sequence(1:5));
			test5 = upper(NewSequenceData(fragment+1).Sequence(1:4));

			if ((strcmp(test1,'GATC') == 1) && (strcmp(test2,'C') == 1) && (strcmp(test3,'AATTG') == 1)) || ...
			   ((strcmp(test4,'AATTG') == 1) && (strcmp(test5,'GATC') == 1))
				NewSequenceData(fragment).Header = [NewSequenceData(fragment).Header ' [*]'];
				fragment_test(fragment) = 1;
			end;
		end;
	end;

	fprintf('\nClearing unusable restriction fragments from digested genome.');
	NewSequenceData(fragment_test == 0) = [];

	outputFile = [workingDir strrep(refFASTA,'fasta','MfeI_MboI.fasta')];

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
	outputFile = [workingDir strrep(refFASTA,'fasta','MfeI_MboI.fasta')];
	fprintf(['\nReference genome was already digested.\nDigested genome is found at :\n\t' outputFile '\n\n']);
end;

end
