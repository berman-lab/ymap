function [crossoverPoint] = FindGaussianCrossover_2(G1,G2, cross)

% FindGaussianCutoffs_2 finds the equal probability point between two skew
%     gaussians.   Equations found by manually solving for the intersection
%     of componant Gaussians of skew Gaussians.
% [a = height; b = location; c = width; d = skew]

method = 2;
	%% 0 : Very simplistic guess.
	%% 1 : Valid only for when Gaussian widths are identical.
	%% 2 : A more general solution for Gaussians with different widths.
	%% 3 : A solution incorporating skew gaussians.

if (method == 0)
	%% Very simplistic guess.
	crossoverPoint = cross+0.5;
elseif (method == 1)
	%% Valid only for when Gaussian widths are identical.
	crossoverPoint = (G1.b+G2.b) - (log(G1.a/G2.a)*4*G1.c^2*G2.c^2 + 2*G1.b^2*G2.c^2 - 2*G2.b^2*G1.c^2)/(4*G1.b*G2.c^2 - 4*G2.b*G1.c^2);
elseif (method == 2)
	%% A more general solution for Gaussians with different widths.
	test_G1a = G1.a
	test_G1b = G1.b
	test_G1c = G1.c
	test_G1d = G1.d
	test_G2a = G2.a
	test_G2b = G2.b
	test_G2c = G2.c
	test_G2d = G2.d

	A = 2*(G2.c^2-G1.c^2);
	B = 4*(G2.b*G1.c^2 - G1.b*G2.c^2);
	C = 2*(G1.b^2*G2.c^2 - G2.b^2*G1.c^2) - log(G1.a/G2.a)*4*G1.c^2*G2.c^2;
	if (A == 0)   % only 1 real answer.
		crossoverPoint = -C/B;
	else   % the closer solution to the wider peak is the most likely interesting solution.
		X1 = (B+(B^2-4*A*C)^(0.5))/(-2*A);
		X2 = (B-(B^2-4*A*C)^(0.5))/(-2*A);
		if (G1.c > G2.c)   %G1 is wider.
			if (abs(G2.b-X1) < abs(G2.b-X2))
				crossoverPoint = X1;
			else
				crossoverPoint = X2;
			end;
		else   % G2 is wider.
			if (abs(G2.b-X1) < abs(G2.b-X2))
				crossoverPoint = X1;
			else
				crossoverPoint = X2;
			end;
		end;
	end;
elseif (method == 3)
	%% A solution incorporating skew gaussians.
	% determine if G1 is skewed at crossover point.
	if (cross > 100)
		cross = cross+1;
	end;
	if (G1.b < 100)
		if (cross < G1.b)
			G1_skew = 0;
		elseif (cross > G1.b)
			G1_skew = 1;
		end;
	elseif (G1.b == 100)
		G1_skew = 0;
	elseif (G1.b > 100)
		if (cross > G1.b)
			G1_skew = 0;
		elseif (cross < G1.b)
			G1_skew = 1;
		end;
	end;
	% determine if G2 is skewed at crossover point.
	if (G2.b < 100)
		if (cross < G2.b)
			G2_skew = 0;
		elseif (cross > G2.b)
			G2_skew = 1;
		end;
	elseif (G2.b == 100)
		G2_skew = 0;
	elseif (G2.b > 100)
		if (cross > G2.b)
			G2_skew = 0;
		elseif (cross < G2.b)
			G2_skew = 1;
		end;
	end;
	% calculate width corrections.
	if (G1_skew == 0) && (G2_skew == 0)
		G1_c = G1.c;
		G2_c = G2.c;
	elseif (G1_skew == 1) && (G2_skew == 0)
		G1_c = G1.c*G1.d/(100.5-abs(100.5-G1.b));
		G2_c = G2.c;
	elseif (G1_skew == 0) && (G2_skew == 1)
		G1_c = G1.c;
		G2_c = G2.c*G2.d/(100.5-abs(100.5-G2.b));
	elseif (G1_skew == 1) && (G2_skew == 1)
		G1_c = G1.c*G1.d/(100.5-abs(100.5-G1.b));
		G2_c = G2.c*G2.d/(100.5-abs(100.5-G2.b));
	end;

	% calculate exact crossover point.
	A = 2*(G2_c^2-G1_c^2);
	B = 4*(G2.b*G1_c^2 - G1.b*G2_c^2);
	C = 2*(G1.b^2*G2_c^2 - G2.b^2*G1_c^2) - log(G1.a/G2.a)*4*G1_c^2*G2_c^2;
	if (A == 0)   % only 1 real answer.
		crossoverPoint = -C/B;
	else   % the closer solution to the 'cross' is the most likely interesting solution.
		X1 = (B+(B^2-4*A*C)^(0.5))/(-2*A);
		X2 = (B-(B^2-4*A*C)^(0.5))/(-2*A);
		if (abs(cross-X1) < abs(cross-X2))
			crossoverPoint = X1;
		else
			crossoverPoint = X2;
		end;
	end;
	crossoverPoint = real(crossoverPoint);
end;

end

