function [which] = FindHighestGaussian_2(G)
% FindHighestGaussian determines which of called Gaussians is the highest
%    across a 1:200 vector.   This is used to determine which Gaussians
%    were meaningful in the multi-Gaussian fit.
	curve_count = length(G);
	which = [];
	% finds highest Gaussian for each of [0..200].
	for i = 1:200
		fit = [];
		num_zeros = 0;
		if (curve_count > 0)
			for c = 1:curve_count
				if (c <= curve_count/2)
					% right-skewed curves.
					if (i <= G{c}.b)
						fit(c) = G{c}.a*exp(-0.5*((i-G{c}.b)./G{c}.c).^2);
					else
						fit(c) = G{c}.a*exp(-0.5*((i-G{c}.b)./G{c}.c/(G{c}.d/(100-abs(100-G{c}.b)))).^2);
					end;
				elseif (c >= curve_count-curve_count/2+1)
					% left-skewed curves.
					if (i < G{c}.b)
						fit(c) = G{c}.a*exp(-0.5*((i-G{c}.b)./G{c}.c/(G{c}.d/(100-abs(100-G{c}.b)))).^2);
					else
						fit(c) = G{c}.a*exp(-0.5*((i-G{c}.b)./G{c}.c).^2);
					end;
				else
					fit(c) = G{c}.a*exp(-0.5*((i-G{c}.b)./G{c}.c).^2);
				end;
				if (fit(c) == 0)
					num_zeros = num_zeros+1;
				end;
			end;
			if (num_zeros == curve_count)
				which(i) = 0;
			else
				fit(fit~=max(fit)) = 0;
				which(i) = find(fit);
			end;
		end;
	end;
	% if the first [which] is zero, finds the nearest available [which] of not zero.
	if (which(1) == 0)
		last_valid = 0;
		for i = 200:-1:1
			if (which(i) ~= 0)
				last_valid = which(i);
			end;
		end;
		which(i) = last_valid;
	end;
	% if the last [which] is zero, finds the nearest available [which] of not zero.
	if (which(200) == 0)
		last_valid = 0;
		for i = 1:1:200
			if (which(i) ~= 0)
				last_valid = which(i);
			end;
		end;
		which(200) = last_valid;
	end;
	% fills [which] gaps of zeros.
	zero_count = 1;
	while (zero_count ~= 0)
		which2 = which;
		zero_count = 0;
		for j = 2:199
			if (which(j) == 0)
				zero_count = zero_count+1;
			end;
			if (which(j) == 0) && (which(j-1) ~= 0) && (which(j+1) == 0)
				which2(j) = which(j-1);
			end;
			if (which(j) == 0) && (which(j-1) == 0) && (which(j+1) ~= 0)
				which2(j) = which(j+1);
			end;
			if (which(j) == 0) && (which(j-1) ~= 0) && (which(j+1) ~= 0)
				%if (random('uniform',0,1) >= 0.5)
					which2(j) = which(j-1);
				%else
				%    which2(j) = which(j+1);
				%end;
			end;
		end;
		which = which2;
		if (zero_count == 198)
			zero_count = 0;
		end;
	end;
end

