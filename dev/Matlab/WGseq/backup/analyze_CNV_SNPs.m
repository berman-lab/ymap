function [] = analyze_CNV_SNPs(main_dir,user,genomeUser,projectChild,projectParent,genome,ploidyEstimateString,ploidyBaseString)

% analyze_CNV_LOHs(<Parent>,<child>,<parentSNP>,<childSNP>,<childCNV>,<Ploidy Estimate)
% A componant of the sequence analysis pipeline, analyzing CNVs only.
%    <Project Name>    : the name of the project.
%    <Project File>    : the name of the project file, including directory location.
%    <Ploidy Estimate> : a numerical value for ploidy from other work.

workingDir      = [main_dir 'users/' user '/projects/' projectChild '/'];
figureDir       = workingDir;
CNV_verString   = 'v1';
INDEL_verString = 'v1';
SNP_verString   = 'v4';
LOH_verString   = 'v2';
displayBREAKS   = false;
referenceCHR    = 1;

fprintf('\n');
fprintf(['[analyze_CNV_LOHs.m data file inputs:\n']);
fprintf(['    projectParent                 : ' projectParent '\n']);
fprintf(['    projectChild                  : ' projectChild '\n']);
fprintf(['    Ploidy                        : ' ploidyEstimateString '\n']);
fprintf('\n');

%% Generate SNP map from allele ratios across genome.
ave_copy_num = 34;
% CNV_SNP_v1(projectParent,projectChild,genome1,ploidyEstimateString,ploidyBaseString,...
%            ave_copy_num,CNV_verString,SNP_verString,LOH_verString,workingDir,figureDir,displayBREAKS);

% Same as above, but tweaked for online pipeline.
CNV_SNP_v2(main_dir,user,genomeUser,projectChild,projectParent,genome,ploidyEstimateString,ploidyBaseString, ...
           CNV_verString,SNP_verString,LOH_verString,displayBREAKS, referenceCHR);

fprintf('*--- End of ''analyze_CNV_LOHs.m'' was reached ---*\n');
% log file end, for in-process analysis
fprintf(['ProjectName : ]]]' projectParent '|||' projectChild ']]]\n']);
end
