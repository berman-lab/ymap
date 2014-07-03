function outputCGDannotationLine(CGDid,probeset1,i,Output_CGD_annotations, color)
%outputCGDuntitled2 outputs a BED-formated line to a text file.
%   BED format is a genome annotation format for uploading to CGD.
if (Output_CGD_annotations == true)
    %fprintf(CGDid,[name '\t"[' probeset1(i).probe_ID(length(probeset1(i).probe_ID)) '/' probeset1(i+1).probe_ID(length(probeset1(i+1).probe_ID)) ']"\t']);
    %if (probeset1(i).probe_chromosome == 8)
    %    fprintf(CGDid,'Ca21chrR');
    %else
    %    fprintf(CGDid,['Ca21chr' num2str(probeset1(i).probe_chromosome)]);
    %end;
    %fprintf(CGDid,['_C_albicans_SC5314:' num2str(probeset1(i).probe_location) '..' num2str(probeset1(i).probe_location) '\n']);
    if (probeset1(i).probe_chromosome == 8)
        fprintf(CGDid,'Ca21chrR');
    else
        fprintf(CGDid,['Ca21chr' num2str(probeset1(i).probe_chromosome)]);
    end;
    fprintf(CGDid,['_C_albicans_SC5314 ' num2str(probeset1(i).probe_location) ' ' num2str(probeset1(i).probe_location) ' [' probeset1(i+1).probe_ID(length(probeset1(i+1).probe_ID)) '/' probeset1(i).probe_ID(length(probeset1(i).probe_ID)) '] 0 + ' num2str(probeset1(i).probe_location) ' ' num2str(probeset1(i).probe_location) ' ' num2str(round(color(1)*255)) ',' num2str(round(color(2)*255)) ',' num2str(round(color(3)*255)) '\n']);
end;
end

