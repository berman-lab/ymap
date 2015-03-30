function [] = analyze_SNPs_hapmap(main_dir, user, genomeUser, project, parent_or_hapmap, genome, ploidyEstimateString,ploidyBaseString)
% analyze_CNVS(<Project Name>,<Project File>,<Ploidy Estimate)
% A componant of the sequence analysis pipeline, analyzing SNPs only.
%    <Project Name>    : the name of the project.
%    <Project File>    : the name of the project file, including directory location.
%    <Ploidy Estimate> : a numerical value for ploidy from other work.

% log file start, for in-process analysis.
fprintf(['project : [[[' project '[[[\n']);

workingDir             = [main_dir 'users/' user '/projects/' project '/'];
figureDir              = workingDir;
CNV_verString          = 'v1';
INDEL_verString        = 'v1';
SNP_verString          = 'v4';
rDNA_verString         = 'v1';
LOH_verString          = 'v2';
displayBREAKS          = true;
referenceCHR           = 1;

%% Grab the first column of the first line of the putative_SNP pileup file.
datafile       = [workingDir 'putative_SNPs_' SNP_verString '.txt'];
data           = fopen(datafile);
fprintf(['datafile : ' datafile]);
line           = fgetl(data);           % grab the first line of the putative_CNV datafile.
fprintf(['\n' datafile '::' num2str(data) '::' line '\n']);
exampleChrName = sscanf(line, '%s',1);  % grab the first column, ex : 'ChrA_C_glabrata_CBS138';
fclose(data);


LOH_hapmap_v4(main_dir,user,genomeUser,project,parent_or_hapmap,genome,ploidyEstimateString,ploidyBaseString,SNP_verString,LOH_verString,CNV_verString,displayBREAKS);

parent = parent_or_hapmap;
hapmap = parent_or_hapmap;
allelic_ratios_WGseq(main_dir,user,genomeUser,project,parent,hapmap,genome,ploidyEstimateString,ploidyBaseString,SNP_verString,LOH_verString,CNV_verString,displayBREAKS);


fprintf('*--- End of ''analyze_SNPs.m'' was reached. ---*\n');
end
