% profile clear;
% profile on;
% analyze_CNVs_RADseq_3('/var/www/html/Ymap/','darren1','default','test_ddRADseq_11461_vs_SC5314_and_novelHapmap_2','test_ddRADseq_SC5314','test_hapmap_2','C_albicans_SC5314_vA21-s02-m09-r07','2.0','2.0');
% analyze_SNPs_RADseq('/var/www/html/Ymap/','darren1','default','test_ddRADseq_11461_vs_SC5314_and_novelHapmap_2','test_ddRADseq_SC5314','test_hapmap_2','C_albicans_SC5314_vA21-s02-m09-r07','2.0','2.0');
analyze_CNV_SNPs_RADseq('/var/www/html/Ymap/','darren1','default','test_ddRADseq_11461_vs_SC5314_and_novelHapmap_2','test_ddRADseq_SC5314','test_hapmap_2','C_albicans_SC5314_vA21-s02-m09-r07','2.0','2.0');
% profile off;
% profile viewer;