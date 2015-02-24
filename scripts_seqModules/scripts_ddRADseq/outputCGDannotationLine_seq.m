function outputCGDannotationLine_seq(CGDid, chrName, coordinate, allele1, allele2, Output_CGD_annotations, color, copyNumberEstimate)
%outputCGDuntitled2 outputs a BED-formated line to a text file.
%   BED format is a genome annotation format for uploading to CGD.
if (Output_CGD_annotations == true)
	% fprintf(['chr = ' chrName '; bp = ' num2str(coordinate) '; copyNum = ' num2str(copyNumberEstimate) '; color = ' num2str(color) '\n']);
	fprintf(CGDid,[chrName ' ' num2str(coordinate) ' ' num2str(coordinate) ' [' allele1 '/' allele2 '] 0 + ' num2str(coordinate) ' ' num2str(coordinate) ' ' num2str(round(color(1)*255)) ',' num2str(round(color(2)*255)) ',' num2str(round(color(3)*255)) '\n']);
end;
end

