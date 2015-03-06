function [] = analyze_SNPs_RADseq(main_dir, user, genomeUser, project, parent, hapmap, genome, ploidyEstimateString,ploidyBaseString)
% analyze_CNVS(<Project Name>,<Project File>,<Ploidy Estimate)
% A componant of the sequence analysis pipeline, analyzing CNVs only.

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

%% Generate SNP map from allele ratios across genome.
% LOH analysis using haplotype map generated by website subsystems.
% LOH_hapmap_v1(main_dir,user,genomeUser,project,hapmap,genome,ploidyEstimateString,ploidyBaseString, SNP_verString,LOH_verString,CNV_verString,displayBREAKS);
% processes phased and unphased data.
% LOH_hapmap_v2(main_dir,user,genomeUser,project,hapmap,genome,ploidyEstimateString,ploidyBaseString, SNP_verString,LOH_verString,CNV_verString,displayBREAKS);


% Adds processing of preprocessed_SNPs file columns for data coordinates.
LOH_hapmap_v3_ddRADseq(main_dir,user,genomeUser,project,hapmap,genome,ploidyEstimateString,ploidyBaseString, SNP_verString,LOH_verString,CNV_verString,displayBREAKS);


% Alternate presentation of SNP data, showing allelic ratios of all data : fire-plot
allelic_ratios_ddRADseq_A(main_dir,user,genomeUser,project,parent,hapmap,genome,ploidyEstimateString,ploidyBaseString, SNP_verString,LOH_verString,CNV_verString,displayBREAKS);


% Alternate presentation of SNP data, showing allelic ratios of all data : ratio-histogram-plot and standard-plot
allelic_ratios_ddRADseq_B(main_dir,user,genomeUser,project,parent,hapmap,genome,ploidyEstimateString,ploidyBaseString, SNP_verString,LOH_verString,CNV_verString,displayBREAKS);


fprintf('*--- End of ''analyze_SNPs.m'' was reached ---*\n');
end
