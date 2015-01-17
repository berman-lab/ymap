function [SNP_count, Ratio_sum] = Collect_SNP_ratios(chr, chr_size, probeset1, projectName);
%% =========================================================================================
% Analyze SNP information for allelic ratios, then gather total and sum for each standard
%    bin across chromosomes.
%===========================================================================================
fprintf(['\nCalculating SNP ratios from project "' projectName '" SNP data.\n']);
SNP_probeset_length = length(probeset1);
bases_per_bin       = max(chr_size)/700;
DataTypeToUse       = 1;   % (1)AllelicFraction; (2)Angle.
show_unnassigned    = false;
SNP_count           = zeros(1,ceil(chr_size(chr)/bases_per_bin));   % SNP count.
Ratio_sum           = zeros(1,ceil(chr_size(chr)/bases_per_bin));   % ratio sum.

for i = 1:2:SNP_probeset_length
	probe_chromosome = probeset1(i).probe_chromosome;
	if (probe_chromosome == chr)
		probe_location   = ceil(probeset1(i).probe_location/bases_per_bin);
		if (length(probeset1(i).probe_Ratio) > 0) && (length(probeset1(i+1).probe_Ratio) > 0)   % Both SNP probes have a valid Red/Green ratio, we can examine this location.
			% Calculate value of SNP probe pair.
			[UsedData] = calculateValue(probeset1,i,DataTypeToUse);

			if (isfield(probeset1(1),'probe_polarity') == 1)
				%%%% Probe entry has a field for 'probe_polarity', meaning the data has been phased vs. haplotype map.
				if (probeset1(i).probe_polarity == 0)
					%%%% probe pair doesn't have probe_polarity assigned in haplotype map.
					if (show_unnassigned == true)
						SNP_count(probe_location) = SNP_count(probe_location) + 1;
						Ratio_sum(probe_location) = Ratio_sum(probe_location) + UsedData;
					end;
				elseif (probeset1(i).probe_polarity == 4)
					%%%% null action for probe design error (probes are identical).
				else
					%%%% Data is haplotype mapped and usable.
					SNP_count(probe_location) = SNP_count(probe_location) + 1;
					Ratio_sum(probe_location) = Ratio_sum(probe_location) + UsedData;
				end;
			else
				%%%% Data has not been phased vs. haplotype map => all data is usable.
				SNP_count(probe_location) = SNP_count(probe_location) + 1;
				Ratio_sum(probe_location) = Ratio_sum(probe_location) + UsedData;
			end;
		end;
	end;
end;

end
