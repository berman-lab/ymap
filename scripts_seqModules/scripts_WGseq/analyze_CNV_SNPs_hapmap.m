function [] = analyze_CNV_SNPs_hapmap(main_dir, user, genomeUser, project, hapmap, genome, ploidyEstimateString,ploidyBaseString)

testPloidy     = ploidyEstimateString
testPloidyBase = ploidyBaseString

workingDir      = [main_dir 'users/' user '/projects/' project '/'];
figureDir       = workingDir;
CNV_verString   = 'v1';
INDEL_verString = 'v1';
SNP_verString   = 'v4';
LOH_verString   = 'v2';
displayBREAKS   = true;
referenceCHR    = 1;

fprintf('\n');
fprintf(['[analyze_CNV_SNPs_hapmap.m data inputs:\n']);
fprintf(['    hapmap  : ' hapmap '\n']);
fprintf(['    project : ' project '\n']);
fprintf(['    Ploidy  : ' ploidyEstimateString '\n']);
fprintf('\n');

%% Generate CNV/SNP/LOH map from allele ratios across genome, using haplotype map information.
% Doesn't display SNP information for locations not included in the hapmap.
% CNV_SNP_hapmap_v1(main_dir,user,genomeUser,project,hapmap,genome,ploidyEstimateString,ploidyBaseString, ...
%                   SNP_verString,LOH_verString,CNV_verString,displayBREAKS);

%% Display hapmap data in colors, and data not included in hapmap in alternate default system.
% CNV_SNP_hapmap_v2(main_dir,user,genomeUser,project,hapmap,genome,ploidyEstimateString,ploidyBaseString, ...
%                   SNP_verString,LOH_verString,CNV_verString,displayBREAKS);

% Adds CGD GBrowse output for SNP data.
% CNV_SNP_hapmap_v3          : old version.
% CNV_SNP_hapmap_v4          : When no hapmap is available, draws SNP ratios as shades of red (with red meaning homozygous).
% CNV_SNP_hapmap_v4_RedGreen : When no hapmap is available, Draws SNP ratios as red:homozygou, grey:heterozygouss and green:non-het/non-hom/unusual-ratios.
% CNV_manualLOH_v1           : Draws CNV plot with manual LOH annotations, if available.

% CNV_SNP_hapmap_v3(main_dir,user,genomeUser,project,hapmap,genome,ploidyEstimateString,ploidyBaseString, ...
%                   SNP_verString,LOH_verString,CNV_verString,displayBREAKS);

CNV_SNP_hapmap_v4(main_dir,user,genomeUser,project,hapmap,genome,ploidyEstimateString,ploidyBaseString, ...
                  SNP_verString,LOH_verString,CNV_verString,displayBREAKS);

CNV_SNP_hapmap_v4_RedGreen(main_dir,user,genomeUser,project,hapmap,genome,ploidyEstimateString,ploidyBaseString, ...
                           SNP_verString,LOH_verString,CNV_verString,displayBREAKS);

CNV_manualLOH_v1( main_dir,user,genomeUser,project,hapmap,genome,ploidyEstimateString,ploidyBaseString, ...
                  SNP_verString,LOH_verString,CNV_verString,displayBREAKS);

fprintf('*--- End of ''analyze_CNV_SNPs_hapmap.m'' was reached ---*\n');
% log file end, for in-process analysis
fprintf(['ProjectName : ]]]project''' project '''|||hapmap''' hapmap ''']]]\n']);
end
