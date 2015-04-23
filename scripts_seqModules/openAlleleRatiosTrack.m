function [alleleRatiosFid] = openAlleleRatiosTrack(projectDir, projectName)

alleleRatiosFid = fopen(fullfile(projectDir, ['allele_ratios.' projectName  '.bed']), 'w');
fprintf(alleleRatiosFid, ['track name=' projectName 'AlleleRatios description="' projectName ' allele ratios" useScore=0 itemRGB=On\n']);

end