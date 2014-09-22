function [] = analyze_CNV_SNPs_hapmap(main_dir, user, genomeUser, project, hapmap, genome, ploidyEstimateString,ploidyBaseString)

testPloidy     = ploidyEstimateString
testPloidyBase = ploidyBaseString

workingDir      = [main_dir 'users/' user '/projects/' project '/'];
figureDir       = workingDir;
CNV_verString   = 'v1';
INDEL_verString = 'v1';
SNP_verString   = 'v4';
LOH_verString   = 'v2';
displayBREAKS   = false;
referenceCHR    = 1;

fprintf('\n');
fprintf(['[analyze_CNV_SNPs_hapmap.m data inputs:\n']);
fprintf(['    hapmap  : ' hapmap '\n']);
fprintf(['    project : ' project '\n']);
fprintf(['    Ploidy  : ' ploidyEstimateString '\n']);
fprintf('\n');

%%	Generate CNV/SNP/LOH map from allele ratios across genome, using haplotype map information.
%	Doesn't display SNP information for locations not included in the hapmap.
% CNV_SNP_hapmap_v1(main_dir,user,genomeUser,project,hapmap,genome,ploidyEstimateString,ploidyBaseString, ...
%                   SNP_verString,LOH_verString,CNV_verString,displayBREAKS);

%%	Display hapmap data in colors, and data not included in hapmap in alternate default system.
%	Second script draws SNP data with red/green color scheme to visualize LOH vs. ratio shift.
% CNV_SNP_hapmap_v2(main_dir,user,genomeUser,project,hapmap,genome,ploidyEstimateString,ploidyBaseString, ...
%                   SNP_verString,LOH_verString,CNV_verString,displayBREAKS);

%%	Adds CGD GBrowse output for SNP data.
%	Third script draws CNV plot with manualLOH annotations, if available.
% CNV_SNP_hapmap_v3(main_dir,user,genomeUser,project,hapmap,genome,ploidyEstimateString,ploidyBaseString, ...
%                   SNP_verString,LOH_verString,CNV_verString,displayBREAKS);

%%	4th
CNV_SNP_hapmap_v4(main_dir,user,genomeUser,project,hapmap,genome,ploidyEstimateString,ploidyBaseString, ...
                  SNP_verString,LOH_verString,CNV_verString,displayBREAKS);

%%	4th oldview changes color scheme for non-homozygous regions from red to green.
%	Attempts to mimic original red/green color scheme used in paper figure.
CNV_SNP_hapmap_v4_oldview(main_dir,user,genomeUser,project,hapmap,genome,ploidyEstimateString,ploidyBaseString, ...
                          SNP_verString,LOH_verString,CNV_verString,displayBREAKS);
CNV_manualLOH_v1( main_dir,user,genomeUser,project,hapmap,genome,ploidyEstimateString,ploidyBaseString, ...
                  SNP_verString,LOH_verString,CNV_verString,displayBREAKS);

fprintf('*--- End of ''analyze_CNV_SNPs_hapmap.m'' was reached ---*\n');
% log file end, for in-process analysis
fprintf(['ProjectName : ]]]project''' project '''|||hapmap''' hapmap ''']]]\n']);
end
