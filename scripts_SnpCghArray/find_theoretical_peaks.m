function [monosomy_peak,disomy_peak,trisomy_peak,tetrasomy_peak,pentasomy_peak,hexasomy_peak ] = ...
    find_theoretical_peaks(realHomozygous_peak)
% find_theoretical_peaks.m determines peak locations when given the homozygous peak
% location.

%% used for drawing locations of peaks in SNP histograms.
fprintf(['realHomozygous_peak = [' num2str(realHomozygous_peak) ']\n']);
monosomy_peak(1)  = realHomozygous_peak;         % 'a'
monosomy_peak(2)  = 1-monosomy_peak(1);          % 'b'
disomy_peak(1)    = realHomozygous_peak;         % 'aa'
disomy_peak(2)    = 0.5;                         % 'ab'
disomy_peak(3)    = 1-disomy_peak(1);            % 'bb'
trisomy_peak(1)   = realHomozygous_peak;         % 'aaa'
trisomy_peak(2)   = (1+realHomozygous_peak)/3;   % 'aab'
trisomy_peak(3)   = 1-trisomy_peak(2);           % 'abb'
trisomy_peak(4)   = 1-trisomy_peak(1);           % 'bbb'
tetrasomy_peak(1) = realHomozygous_peak;         % 'aaaa'
tetrasomy_peak(2) = (0.5+realHomozygous_peak)/2; % 'aaab'
tetrasomy_peak(3) = 0.5;                         % 'aabb'
tetrasomy_peak(4) = 1-tetrasomy_peak(2);         % 'abbb'
tetrasomy_peak(5) = 1-tetrasomy_peak(1);         % 'bbbb'
pentasomy_peak(1) = realHomozygous_peak;         % 'aaaaa'
pentasomy_peak(2) = 0.2+realHomozygous_peak*3/5; % 'aaaab'
pentasomy_peak(3) = 0.4+realHomozygous_peak/5;   % 'aaabb'
pentasomy_peak(4) = 1-pentasomy_peak(3);         % 'aabbb'
pentasomy_peak(5) = 1-pentasomy_peak(2);         % 'abbbb'
pentasomy_peak(6) = 1-pentasomy_peak(1);         % 'bbbbb'
hexasomy_peak(1) = realHomozygous_peak;          % 'aaaaaa'
hexasomy_peak(2) = 1/6+realHomozygous_peak*2/3;  % 'aaaaab'
hexasomy_peak(3) = 1/3+realHomozygous_peak/3;    % 'aaaabb'
hexasomy_peak(4) = 0.5;                          % 'aaabbb'
hexasomy_peak(5) = 1-hexasomy_peak(3);           % 'aabbbb'
hexasomy_peak(6) = 1-hexasomy_peak(2);           % 'abbbbb'
hexasomy_peak(7) = 1-hexasomy_peak(1);           % 'bbbbbb'
end
