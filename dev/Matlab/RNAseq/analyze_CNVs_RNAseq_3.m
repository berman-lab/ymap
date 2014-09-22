function [] = analyze_CNVs_RNAseq_3(main_dir,user,genomeUser,project,parent,hapmap,genome,ploidyEstimate,ploidyBase)
% A componant of the sequence analysis pipeline, analyzing CNVs in RNAseq data only.
%    <main_dir>       : base directory of web pipeline.
%    <user>           : name of user.
%    <genomeUser>     : name of user for genome.
%    <project>        : name of project.
%    <projectParent>  : name of parent project.
%    <genome>         : name of genome in use.
%    <ploidyEstimate> : numerical estimate of strain ploidy.
%    <ploidyBase>     : baseline ploidy to draw figures from.

% log file start, for in-process analysis.
fprintf(['project : [[[' project '[[[\n']);

CNV_verString   = 'v1';
INDEL_verString = 'v1';
SNP_verString   = 'v4';
rDNA_verString  = 'v1';
displayBREAKS   = true;
referenceCHR    = 1;

%%% Attempting to resolve error : final normalization vs. reference dataset doesn't work as in MSI pipeline.
%CNV_v6_fragmentLengthCorrected_6(main_dir,user,genome,genomeUser,project,projectParent,projectParentUser,ploidyEstimate,ploidyBase, ...
%                                 CNV_verString,displayBREAKS);

%% ORF length bias correction.
%% GC bias correction.
%% Repetitiveness bias correction returned.
%% Chromosome end bias correction added.
CNV_v6_fragmentLengthCorrected_8(main_dir,user,genomeUser,project,parent,genome,ploidyEstimate,ploidyBase, ...
                                 CNV_verString,displayBREAKS);


fprintf('*--- End of ''analyze_CNVs_RNAseq_1.m'' was reached ---*\n');
% log file end, for in-process analysis.
fprintf(['project : ]]]' project ']]]\n']);
end
