function [UsedData] = calculateValue(probeset1,i,DataTypeToUse)
if (isfield(probeset1(1),'probe_polarity') == 1)
    % Calculate allelic fraction for SNP pair.
    if (probeset1(i).probe_polarity == 1)
        AllelicFraction = probeset1(i+1).probe_Ratio/(probeset1(i).probe_Ratio+probeset1(i+1).probe_Ratio);
    elseif (probeset1(i).probe_polarity == 2)
        AllelicFraction = probeset1(i).probe_Ratio/(probeset1(i).probe_Ratio+probeset1(i+1).probe_Ratio);
    else
        AllelicFraction = 0;
    end;
    % Calculate angle for SNP pair.
    if (probeset1(i).probe_polarity == 1)
        Angle = atan2(probeset1(i+1).probe_Ratio,probeset1(i).probe_Ratio)/(pi/2);
    elseif (probeset1(i).probe_polarity == 2)
        Angle = atan2(probeset1(i).probe_Ratio,probeset1(i+1).probe_Ratio)/(pi/2);
    else
        Angle = 0;
    end;
else
    AllelicFraction = probeset1(i+1).probe_Ratio/(probeset1(i).probe_Ratio+probeset1(i+1).probe_Ratio);
    Angle = atan2(probeset1(i+1).probe_Ratio,probeset1(i).probe_Ratio)/(pi/2);
end;
% Determine which data type to use.
if (DataTypeToUse == 1)   % use AllelicFractions.
    UsedData = AllelicFraction;
else % (DataTypeToUse == 2)   % use Angles.
    UsedData = Angle;
end;
end

