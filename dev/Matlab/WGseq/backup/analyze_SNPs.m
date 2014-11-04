function [] = analyze_SNPs(main_dir, user, genomeUser, projectChild, projectParent, genome, ploidyEstimateString,ploidyBaseString)
% analyze_CNVS(<Project Name>,<Project File>,<Ploidy Estimate)
% A componant of the sequence analysis pipeline, analyzing CNVs only.
%    <Project Name>    : the name of the project.
%    <Project File>    : the name of the project file, including directory location.
%    <Ploidy Estimate> : a numerical value for ploidy from other work.

% log file start, for in-process analysis.
fprintf(['project : [[[' projectChild '[[[\n']);

workingDir             = [main_dir 'users/' user '/projects/' projectChild '/'];
figureDir              = workingDir;
CNV_verString          = 'v1';
INDEL_verString        = 'v1';
SNP_verString          = 'v4';
rDNA_verString         = 'v1';
LOH_verString          = 'v2';
displayBREAKS          = false;
referenceCHR           = 1;

%% Grab the first column of the first line of the putative_SNP pileup file.
datafile       = [workingDir 'putative_SNPs_' SNP_verString '.txt'];
data           = fopen(datafile);
fprintf(['datafile : ' datafile]);
line           = fgetl(data);           % grab the first line of the putative_CNV datafile.
fprintf(['\n' datafile '::' num2str(data) '::' line '\n']);
exampleChrName = sscanf(line, '%s',1);  % grab the first column, ex : 'ChrA_C_glabrata_CBS138';
fclose(data);

%% Generate SNP map from allele ratios across genome.
full_data_threshold = 50;

% LOH_v2_1(project,project,genome,full_data_threshold,SNP_verString,LOH_verString,workingDir,figureDir,displayBREAKS);

% Same as above, with better annotation display.
% LOH_v2_2(project,project,genome,full_data_threshold,SNP_verString,LOH_verString,workingDir,figureDir,displayBREAKS);

% Same as above, but tuned to online pipeline.
% LOH_v2_3(main_dir,user,genomeUser,projectChild,projectParent,genome,ploidyEstimateString,ploidyBaseString, ...
%          SNP_verString,LOH_verString,CNV_verString,displayBREAKS);

% Same as above, but tuned to online pipeline and preprocessing via python script : dataset_process_for_SNP_analysis.py
LOH_v2_4(main_dir,user,genomeUser,projectChild,projectParent,genome,ploidyEstimateString,ploidyBaseString, ...
         SNP_verString,LOH_verString,CNV_verString,displayBREAKS);

fprintf('*--- End of ''analyze_SNPs.m'' was reached ---*\n');
end
