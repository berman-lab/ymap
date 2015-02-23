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

% workingDir = projectDirectory;
% figureDir  = projectDirectory;

%% Generate CNV map from read count across genome.
% CNV_v6_1(project,genome,ploidyEstimateString, ...
%          CNV_verString,rDNA_verString,workingDir,figureDir,displayBREAKS, referenceCHR);
% CNV_v6_2(project,genome,ploidyEstimateString,ploidyBaseString, ...
%          CNV_verString,rDNA_verString,workingDir,figureDir,displayBREAKS, referenceCHR);

% % Same as above, but improved annotation display.
% CNV_v6_3(project,genome,ploidyEstimateString,ploidyBaseString, ...
%          CNV_verString,rDNA_verString,workingDir,figureDir,displayBREAKS, referenceCHR);

% % Same as above, but adjusted for online pipeline.
% CNV_v6_4(main_dir,user,genomeUser,project,genome,ploidyEstimateString,ploidyBaseString, ...
%          CNV_verString,rDNA_verString,displayBREAKS, referenceCHR);

% % Same as above, but adjusted for online pipeline, including using a python preprocessed putative_CNVs file to speed things up.
% CNV_v6_5(main_dir,user,genomeUser,project,genome,ploidyEstimateString,ploidyBaseString, ...
%          CNV_verString,rDNA_verString,displayBREAKS, referenceCHR);

% CNV_v6_6         : Repetitiveness bias correction and Chromosome end bias correction added.
% CNV_v6_6_highTop : re-renders output from CNV_v6_6.m, with a higher max 3x the height of the chromosome.
% CNV_v6_7         : not used, includes repetitiveness correction.
CNV_v6_6(main_dir,user,genomeUser,project,genome,ploidyEstimateString,ploidyBaseString, ...
         CNV_verString,rDNA_verString,displayBREAKS, referenceCHR);
CNV_v6_6_highTop(main_dir,user,genomeUser,project,genome,ploidyEstimateString,ploidyBaseString, ...
                 CNV_verString,rDNA_verString,displayBREAKS, referenceCHR);

%% rDNA copy estimation.
%if (strcmp(genome,'Ca_a') == 1)
%	load([projectDirectory 'data.rDNA-CNV_' rDNA_verString '.mat']);
%	rDNA = countRDNA/sizeRDNA;
%	ref  = countREF/sizeREF;
%	fid  = fopen([projectDirectory 'data.rDNA-CNV_analysis.txt'],'w');
%	fprintf(fid,['Project ''' project ''' has an estimated [' num2str(rDNA/ref) '] rDNA copies per genome copy.\n']);
%	fprintf(fid, '--------\n');
%	fprintf(fid, 'rDNA copy number is calculated by comparing the average read count in the rDNA region to the average read count across chr1.\n');
%	fprintf(fid, 'If chr1 represents an aneuploidy, then rDNA copy number estimates will be need to be adjsted.\n');
%	fprintf(fid, 'This analysis is performed as a subset of general CNV analysis.\n');
%	fprintf(fid, 'Currently this analysis only supports Candida albicans.');
%	fclose(fid);
%end;

fprintf('*--- End of ''analyze_CNVs_1.m'' was reached ---*\n');
end
