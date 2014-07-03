function [monosomy_cutoff,disomy_cutoff,trisomy_cutoff,tetrasomy_cutoff,pentasomy_cutoff,hexasomy_cutoff ] = ...
    FindCutoffs( monosomy_peak,disomy_peak,trisomy_peak,tetrasomy_peak,pentasomy_peak,hexasomy_peak )
% FindPeaks determines peak locations when given the homozygous peak location.

%% used for drawing locations of cutoffs in SNP histograms.
monosomy_cutoff(1)  = (monosomy_peak(1)+monosomy_peak(2))/2;   % 'a':'b'
disomy_cutoff(1)    = (disomy_peak(1)+disomy_peak(2))/2;       % 'aa':'ab'
disomy_cutoff(2)    = (disomy_peak(2)+disomy_peak(3))/2;       % 'ab':'bb'
trisomy_cutoff(1)   = (trisomy_peak(1)+trisomy_peak(2))/2;     % 'aaa':'aab'
trisomy_cutoff(2)   = (trisomy_peak(2)+trisomy_peak(3))/2;     % 'aab':'abb'
trisomy_cutoff(3)   = (trisomy_peak(3)+trisomy_peak(4))/2;     % 'abb':'bbb'
tetrasomy_cutoff(1) = (tetrasomy_peak(1)+tetrasomy_peak(2))/2; % 'aaaa':'aaab'
tetrasomy_cutoff(2) = (tetrasomy_peak(2)+tetrasomy_peak(3))/2; % 'aaab':'aabb'
tetrasomy_cutoff(3) = (tetrasomy_peak(3)+tetrasomy_peak(4))/2; % 'aabb':'abbb'
tetrasomy_cutoff(4) = (tetrasomy_peak(4)+tetrasomy_peak(5))/2; % 'abbb':'bbbb'
pentasomy_cutoff(1) = (pentasomy_peak(1)+pentasomy_peak(2))/2; % 'aaaaa':'aaaab'
pentasomy_cutoff(2) = (pentasomy_peak(2)+pentasomy_peak(3))/2; % 'aaaab':'aaabb'
pentasomy_cutoff(3) = (pentasomy_peak(3)+pentasomy_peak(4))/2; % 'aaabb':'aabbb'
pentasomy_cutoff(4) = (pentasomy_peak(4)+pentasomy_peak(5))/2; % 'aabbb':'abbbb'
pentasomy_cutoff(5) = (pentasomy_peak(5)+pentasomy_peak(6))/2; % 'abbbb':'bbbbb'
hexasomy_cutoff(1)  = (hexasomy_peak(1)+hexasomy_peak(2))/2;   % 'aaaaaa':'aaaaab'
hexasomy_cutoff(2)  = (hexasomy_peak(2)+hexasomy_peak(3))/2;   % 'aaaaab':'aaaabb'
hexasomy_cutoff(3)  = (hexasomy_peak(3)+hexasomy_peak(4))/2;   % 'aaaabb':'aaabbb'
hexasomy_cutoff(4)  = (hexasomy_peak(4)+hexasomy_peak(5))/2;   % 'aaabbb':'aabbbb'
hexasomy_cutoff(5)  = (hexasomy_peak(5)+hexasomy_peak(6))/2;   % 'aabbbb':'abbbbb'
hexasomy_cutoff(6)  = (hexasomy_peak(6)+hexasomy_peak(7))/2;   % 'abbbbb':'bbbbbb'end
