function  analyze_CNVs_1(main_dir, user, genomeUser, project, genome, ploidyEstimateString,ploidyBaseString)
% analyze_CNVS(<Project Name>,<Project File>,<Ploidy Estimate)
% A componant of the sequence analysis pipeline, analyzing CNVs only.
%    <main_dir>             : the main pipeline directory.
%    <user>                 : the user account owner of the project.
%    <genomeUser>           : the user account owner of the genome used in this project.
%    <project>              : the name of the project.
%    <genome>               : the name of the genome used in this project.
%    <ploidyEstimateString> : ploidy determined elsewhere for current dataset.
%    <ploidyBaseString>     : exepected ploidy for species.

% log file start, for in-process analysis.
fprintf(['project : [[[' project '[[[\n']);

projectDirectory       = [main_dir '/users/' user '/projects/' project '/'];
genomeDirectory        = [main_dir '/users/' genomeUser '/genomes/' genome '/'];
CNV_verString          = 'v1';
INDEL_verString        = 'v1';
SNP_verString          = 'v4';
rDNA_verString         = 'v1';
displayBREAKS          = true;
referenceCHR           = 1;


CNV_v6_6(main_dir,user,genomeUser,project,genome,ploidyEstimateString,ploidyBaseString,CNV_verString,rDNA_verString,displayBREAKS, referenceCHR);

CNV_v6_6_highTop(main_dir,user,genomeUser,project,genome,ploidyEstimateString,ploidyBaseString,CNV_verString,rDNA_verString,displayBREAKS, referenceCHR);


fprintf('*--- End of ''analyze_CNVs_1.m'' was reached ---*\n');
end
