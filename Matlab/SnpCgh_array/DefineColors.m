function [ colorA,colorB,...
           colorAB,...
           colorAAB,colorABB,...
           colorAAAB,colorABBB,...
           colorAAAAB,colorAAABB,colorAABBB,colorABBBB,...
           colorAAAAAB,colorABBBBB,...
           colorPeak,colorCutoff ] = DefineColors()
% DefineColors defines colors used in the main script.
colorAB     = [0.667 0.667 0.667]; % heterozygous.
colorA      = [1.0   0.0   1.0  ]; % homozygous a:magenta.
colorB      = [0.0   1.0   1.0  ]; % homozygous b:cyan.
colorAAB    = [0.667 0.333 1.0  ]; % heterozygous trisomy, aab.
colorABB    = [0.333 0.667 1.0  ]; % heterozygous trisomy, abb.
colorAAAB   = [0.75   0.25   1.0]; % heterozygous tetrasomy, aaab.
colorABBB   = [0.25   0.75   1.0]; % heterozygous tetrasomy, abbb.
colorAAAAB  = [0.8    0.2    1.0]; % heterozygous pentasomy, aaaab.
colorAAABB  = [0.6    0.4    1.0]; % heterozygous pentasomy, aaabb.
colorAABBB  = [0.4    0.6    1.0]; % heterozygous pentasomy, aabbb.
colorABBBB  = [0.2    0.8    1.0]; % heterozygous pentasomy, abbbb.
colorAAAAAB = [0.8333 0.1667 1.0]; % heterozygous hexasomy, aaaaab.
colorABBBBB = [0.1667 0.8333 1.0]; % heterozygous hexasomy, abbbbb.
colorPeak   = [0.5   0.5   0.5  ]; % peak locations.
colorCutoff = [1.0   0.0   0.0  ]; % cutoff locations.
end