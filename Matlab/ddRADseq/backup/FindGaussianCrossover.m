function [crossoverPoint] = FindGaussianCrossover(G1,G2)
% FindGaussianCutoffs finds the equal probability point between two
%     gaussians.   Equation was found my manually solving for the intersection
%     of two gaussians.   [a = height; b = location; c = width]
% Valid only for when Gaussian widths are identical.
crossoverPoint1 = (G1.b+G2.b) - (log(G1.a/G2.a)*4*G1.c^2*G2.c^2 + 2*G1.b^2*G2.c^2 - 2*G2.b^2*G1.c^2)/(4*G1.b*G2.c^2 - 4*G2.b*G1.c^2);
% A more general solution for Gaussians with different widths.
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

end


