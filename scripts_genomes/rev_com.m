% DNA sequence tool.
%   Returns reverse complement of input DNA sequence.
function [rev_com_seq] = rev_com(seq)
   % reverses input sequence.
   rev_seq = fliplr(seq);

   % initialize rev_com_seq to rev_seq.   this allows un-resolved
   % characters in the sequences to remain unresolved as space-fillers.
   rev_com_seq = rev_seq;

   % complements input sequence.
   rev_com_seq(rev_seq == 'A') = 'T';
   rev_com_seq(rev_seq == 'T') = 'A';
   rev_com_seq(rev_seq == 'C') = 'G';
   rev_com_seq(rev_seq == 'G') = 'C';
   rev_com_seq(rev_seq == 'a') = 'T';
   rev_com_seq(rev_seq == 't') = 'A';
   rev_com_seq(rev_seq == 'c') = 'G';
   rev_com_seq(rev_seq == 'g') = 'C';
end
