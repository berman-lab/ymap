function [] = analyze_CNV_SNPs_hapmap(main_dir, user, genomeUser, project, hapmap, genome, ploidyEstimateString,ploidyBaseString)

testPloidy     = ploidyEstimateString
testPloidyBase = ploidyBaseString

workingDir      = [main_dir 'users/' user '/projects/' project '/'];
figureDir       = workingDir;
INDEL_verString = 'v1';
SNP_verString   = 'v4';
LOH_verString   = 'v2';
CNV_verString   = 'v1';
displayBREAKS   = true;
referenceCHR    = 1;

fprintf('\n');
fprintf(['[analyze_CNV_SNPs_hapmap.m data inputs:\n']);
fprintf(['    hapmap  : ' hapmap '\n']);
fprintf(['    project : ' project '\n']);
fprintf(['    Ploidy  : ' ploidyEstimateString '\n']);
fprintf('\n');


CNV_SNP_hapmap_v4(main_dir,user,genomeUser,project,hapmap,genome,ploidyEstimateString,ploidyBaseString,SNP_verString,LOH_verString,CNV_verString,displayBREAKS);

CNV_SNP_hapmap_v4_RedGreen(main_dir,user,genomeUser,project,hapmap,genome,ploidyEstimateString,ploidyBaseString,SNP_verString,LOH_verString,CNV_verString,displayBREAKS);

CNV_manualLOH_v1( main_dir,user,genomeUser,project,hapmap,genome,ploidyEstimateString,ploidyBaseString,SNP_verString,LOH_verString,CNV_verString,displayBREAKS);


fprintf('*--- End of ''analyze_CNV_SNPs_hapmap.m'' was reached ---*\n');
% log file end, for in-process analysis
fprintf(['ProjectName : ]]]project''' project '''|||hapmap''' hapmap ''']]]\n']);
end
