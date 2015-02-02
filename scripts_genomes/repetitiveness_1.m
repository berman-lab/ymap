function results = repetitiveness_1(user,genome)
% This script attempts to analyze genomic sequences for repetitive elements.
% Loads multiple FASTA-formatted sequence files defining a genome, then processes it using the repetitiveness algorithm.

% log file start, for in-process analysis.
fprintf(['** Genome : [[[' genome '[[[\n']);

workingDir             = ['../users/' user '/genomes/' genome '/'];
nmer_length            = 10;

[centromeres, chr_sizes, figure_details, annotations, ploidy_default] = Load_genome_information_1(workingDir,genome);

%% Determine number of chromosomes from figure_details.
num_chr    = 0;
chr_labels = [];
for i = 1:length(figure_details)
    if (figure_details(i).chr ~= 0)
        num_chr = num_chr+1;
        chr_names{figure_details(i).chr}  = figure_details(i).name;
    end;
end;

%% Determine reference genome FASTA file in use.
%  Read in and parse : "links_dir/main_script_dir/genome_specific/[genome]/reference.txt"
reference_file = [workingDir '/reference.txt'];
refernce_fid   = fopen(reference_file, 'r');
refFASTA       = fgetl(refernce_fid);
fclose(refernce_fid);
FastaName      = strrep(refFASTA,'.fasta','');

fprintf(['\n** number of chrs in figure = ' num2str(num_chr)]);

%% ====================================================================
% Parse sequences from FASTQ file, or load MAT files if already generated for this genome.
% ---------------------------------------------------------------------
if (exist([workingDir FastaName '.seq.mat'],'file') == 0) && (exist([workingDir FastaName '.seqRevCom.mat'],'file') == 0)
    fprintf('\n** Parsing sequences from FASTQ file.');
    SequenceData = fastaread([workingDir refFASTA ],'TrimHeaders', true);
    num_chr = length(SequenceData);

    %% Initialize the cell arrays needed to contain chromosome sequences.
    sequences         = cell(1,num_chr);
    rev_com_sequences = cell(1,num_chr);

    fprintf(['\n**\tLoading FASTA file: \"' refFASTA '\"\n']);
	for i = 1:length(SequenceData)
		sequences{i}         = SequenceData(i).Sequence;
		name                 = SequenceData(i).Header;
		fprintf([num2str(i) ' : ' name '\n']);
		name_strings         = textscan(name,'%s').';
		names{i}             = name_strings{1};
		rev_com_sequences{i} = rev_com(sequences{i});
    end;
    save([workingDir FastaName '.seq.mat']      ,'sequences', 'names');
    save([workingDir FastaName '.seqRevCom.mat'],'rev_com_sequences');
else
    fprintf('\n** Loading pre-determined sequences.');
    load([workingDir FastaName '.seq.mat']);
    load([workingDir FastaName '.seqRevCom.mat']);
end;

for i = 1:length(sequences)
    length_nts(i) = length(sequences{i});
end;

%% ====================================================================
% Tally incidence of nmers across genome.
% ---------------------------------------------------------------------
num_chr = length(sequences);
fprintf(['\n**\t' num2str(num_chr) ' chromosomes in genome sequence file.']);
if (exist([workingDir FastaName '.nmerCounts_' num2str(nmer_length) '.mat'],'file') == 0)
    fprintf('\n** Determining incidence of nmers across across genome.\n');
    nmer_counts = zeros(1,4^nmer_length,'int32');
    for chr = 1:num_chr
        Query_length = nmer_length;

	fprintf(['\n**\tAnalyzing ' num2str(nmer_length) 'mer repeats in : ' names{chr}{1}]);
%	disp(names{chr});

	fprintf(['\n**\t\tseq length = ' num2str(length_nts(chr))]);

	for Query_offset = 1:(length_nts(chr)-Query_length+1)
	    Query         = sequences{chr}(Query_offset:(Query_offset+Query_length-1));
	    rev_com_Query = rev_com_sequences{chr}((length_nts(chr)-Query_offset-Query_length+2):(length_nts(chr)-Query_offset+1));

	    % Add to counts of query sequence and reverse complement of query sequence.
	    [forward,err1] = find_nmer(Query);
	    [reverse,err2] = find_nmer(rev_com_Query);
	    if (err1 == false)
		nmer_counts(forward) = nmer_counts(forward)+1;
	    end;
	    if (err2 == false)
		nmer_counts(reverse) = nmer_counts(reverse)+1;
	    end;
	end;
    end;
    % full genome analysis for 5mers took:   27,253.760544 s = 454.229342 min = 7.570489 hr
    % and then I wrote 'find_nmers.m'...
    % full genome analysis for 10mers took:   1,093.824293 s =  18.230405 min = 0.303840 hr (~18 minutes)

    fprintf('\n** Saving incidence of nmers across genome.');
    save([workingDir FastaName '.nmerCounts_' num2str(nmer_length) '.mat'],'nmer_counts');
    clear counter forward reverse Query rev_com_Query Query_length;
    clear Query_offset err1 err2;
else
    fprintf('\n** Loading incidence of nmers across genome.');
    load([workingDir FastaName '.nmerCounts_' num2str(nmer_length) '.mat']);
end;

%% ====================================================================
% Determine distribution of repeats across genome.
% ---------------------------------------------------------------------
if (exist([workingDir FastaName '.nmerDists_' num2str(nmer_length) '.mat'],'file') == 0)
    fprintf('\n** Determining distribution of repeats across genome.');
    % Make a vectors to track incidence of repetitive sequences.
    clear repeats;
    repeats = cell(1,num_chr);
    for chr = 1:num_chr
        repeats{chr}(1:length_nts(chr)) = 0;
    end;
    
    % Examine sequence...
    for chr = 1:num_chr
        counter = 0;
        Query_length = nmer_length;
        fprintf(['\n**\tAnalyzing ' num2str(nmer_length) 'mer repeats in : ' names{chr}{1}]);
%		disp(names{chr});
        for Query_offset = 1:(length_nts(chr)-Query_length+1)
            counter = counter+1;
            Query = sequences{chr}(Query_offset:(Query_offset+Query_length-1));
            [forward, err1] = find_nmer(Query);
            if (err1 == false)
                if (~isempty(forward))
                    for j = Query_offset:(Query_offset+Query_length-1)
                        repeats{chr}(j) = repeats{chr}(j) + nmer_counts(forward);
                    end;
                end;
            end;
        end;
    end;
    fprintf('\n** Saving distribution of repeats across genome.');
    save([workingDir FastaName '.nmerDists_' num2str(nmer_length) '.mat'],'repeats');
    clear i j counter forward reverse Query Query_length Query_offset;
    clear rev_com_Query err1;
else
    fprintf('\n** Loading distribution of repeats across genome.');
    load([workingDir FastaName '.nmerDists_' num2str(nmer_length) '.mat']);
    for chr = 1:num_chr
        length_nts(chr) = length(repeats{chr});
    end;
    clear input_files sequences rev_com_sequences chr_names i;
end;


%% ====================================================================
% Output raw repetitiveness score per base pair.
% ---------------------------------------------------------------------
fprintf('\n** Outputing repetitiveness score file for genome.');

outFile = [workingDir FastaName '.repetitiveness.txt'];
fid = fopen(outFile,'w');
fprintf(fid,'## Repetitiveness score per bp location.\n');
fprintf(fid,'##\n');
fprintf(fid,'## columns = [chrName, bpCoordinate, repetitivenessScore]\n');
for chr = 1:num_chr
	fprintf(['\n**\t' names{chr}{1}]);
%	disp(names{chr}{1});
	for bp = 1:length_nts(chr)
		repetitiveness_value = repeats{chr}(bp);
		fprintf(fid, names{chr}{1} );
		fprintf(fid,['\t' num2str(bp) '\t' num2str(repetitiveness_value) '\n']);
 	end;
end;
fclose(fid);
fprintf('\n**\n');


%% ====================================================================
% Smooth repetitiveness score & output
% ---------------------------------------------------------------------
%for smooth_level = 1:20
%    sigma = smooth_level*1000;
%    fprintf(['\nSmoothing raw repetitiveness plot (sigma = ' num2str(sigma) ') [']);
%    for chr = 1:num_chr
%	fprintf(num2str(chr));
%	% smooth_gaussian(data,sigma,size)
%	%    data  : input vector with raw data.
%	%    sigma : standard deviation of the gaussian distribution used in the smoothing.
%	%    size  : size of vector over which smoothing function is applied.   (2-3 sigmas is usually good.)
%	s1{chr} = smooth_gaussian(repeats{chr},sigma,sigma*2);
%    end;
%    fprintf(']');
%
%    outFile = [workingDir FastaName '.repetitiveness.sigma_' num2str(sigma) '.txt'];
%    fid = fopen(outFile,'w');
%    fprintf(fid,'## Repetitiveness score per bp location.\n');
%    fprintf(fid,['## Raw scores are Gaussian smoothed with sigma = ' num2str(sigma) '.\n']);
%    fprintf(fid,'##\n');
%    fprintf(fid,'## columns = [chrName, bpCoordinate, repetitivenessScore]\n');
%    for chr = 1:num_chr
%	fprintf(['\n**\t' names{chr}{1}]);
%	for bp = 1:length_nts(chr)
%	    repetitiveness_value = repeats{chr}(bp);
%	    fprintf(fid,[names{chr} '\t' num2str(bp) '\t' num2str(repetitiveness_value) '\n']);
%	end;
%    end;
%    fclose(fid);
%    fprintf('\n**\n');
%end;

end
