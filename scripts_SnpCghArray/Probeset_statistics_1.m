%% ========================================================================
% Generate statistics based on probesets and calibration data.
% Generate tab-delimited-text file output of hapmap.
%==========================================================================
clear;

%% ========================================================================
% Define SNP probe pair polarity.
%==========================================================================

% ask the user which calibration file to use... or to use datasets.
v = 0;
while (v == 0)
    d   = dir;
    str = {d.name};
    j = 0;
    results = [];
    for i = length(str):-1:1
        if (isnan(strfind(str{i}, 'calibration')) == 0)
            j = j+1;
            results{j} = str{i};
        end;
    end;
    [s,v] = listdlg('PromptString','Choose calibration file:',...
        'SelectionMode','single',...
        'ListString',results);
end;
load(results{s});
load SNP_probeset_2.mat;
SNP_probeset_length = length(SNP_probeset);

% Assign SNP probes to homologs, based on calibration data collected earlier.
for i = 1:2:SNP_probeset_length
    if (SNP_probeset(i).assign_cyan > SNP_probeset(i).assign_magenta)
        SNP_probeset(i  ).probe_polarity = 1;
        SNP_probeset(i+1).probe_polarity = 1;
    elseif (SNP_probeset(i).assign_cyan < SNP_probeset(i).assign_magenta)
        SNP_probeset(i  ).probe_polarity = 2;
        SNP_probeset(i+1).probe_polarity = 2;
    elseif (SNP_probeset(i).assign_cyan == SNP_probeset(i).assign_magenta)
        SNP_probeset(i  ).probe_polarity = 0;
        SNP_probeset(i+1).probe_polarity = 0;
    end;
end;

% setup hapmap tab-delimited-text file.
resultsFileName = results{s};
for i = length(resultsFileName):-1:length(resultsFileName)-2
    resultsFileName(i) = [];
end;
for i = 11:-1:1
    resultsFileName(i) = [];
end;
resultsFileID = fopen(['hapmap' resultsFileName 'txt'], 'w');
fprintf(resultsFileID, 'SNP_chr\tSNP_bp\tprobe_A_allele\tprobe_B_allele\tmapped');

count_0a = 0;  % unassigned counter, no data.
count_0b = 0;  % unassigned counter, data.
count_1 = 0;  % a:magenta counter.
count_2 = 0;  % b:cyan counter.
count_4 = 0;  % probe identity counter.
for i = 1:2:SNP_probeset_length
    % determine offset between probes.
    string1 = probeset_2(i).probe_sequence;
    string2 = probeset_2(i+1).probe_sequence;
    counter = 0;
    for ii = 1:length(string2)
        counter = counter+1;
        string3 = '';
        for jj = 1:ceil(length(string1)/2);   string3 = [string3 '-'];   end;
        string3 = [string3 string2];
        string4 = '';
        for jj = 1:(counter);   string4 = [string4 '-'];   end;
        string4 = [string4 string1];
        count_ident(ii) = 0;
        if (length(string4)<length(string3))
            for jj = 1:length(string4)
                if (string4(jj) ~= '-')
                    if (string4(jj) == string3(jj))
                        count_ident(ii) = count_ident(ii)+1;
                    end;
                end;
            end;
        else
            for jj = 1:length(string3)
                if (string3(jj) ~= '-')
                    if (string4(jj) == string3(jj))
                        count_ident(ii) = count_ident(ii)+1;
                    end;
                end;
            end;
        end;
    end;
    Offset = find(count_ident==max(count_ident))-ceil(length(string1)/2);
    string1_ = '';
    string2_ = '';
    %i
    if (Offset >= 0)
        for jj = 1:length(string1)
            if (jj+Offset <= length(string2))
                if (string1(jj) == string2(jj+Offset))
                    string1_ = [string1_ string1(jj)];
                else
                    string1_ = [string1_ '[' string1(jj) ']'];
                end;
            else
                string1_ = [string1_ string1(jj)];
            end;
        end;
        for jj = 1:length(string2)
            if (jj+Offset <= length(string1)) && (jj-Offset > 0)
                if (string2(jj) == string1(jj-Offset))
                    string2_ = [string2_ string2(jj)];
                else
                    string2_ = [string2_ '[' string2(jj) ']'];
                end;
            else
                string2_ = [string2_ string2(jj)];
            end;
        end;
    else % (Offset < 0)
        for jj = 1:length(string1)
            if (jj+Offset > 0) && (jj+Offset <= length(string2))
                if (string1(jj) == string2(jj+Offset))
                    string1_ = [string1_ string1(jj)];
                else
                    string1_ = [string1_ '[' string1(jj) ']'];
                end;
            else
                string1_ = [string1_ string1(jj)];
            end;
        end;
        for jj = 1:length(string2)
            if (jj-Offset <= length(string1)) && (jj+Offset > 0)
                if (string2(jj) == string1(jj-Offset))
                    string2_ = [string2_ string2(jj)];
                else
                    string2_ = [string2_ '[' string2(jj) ']'];
                end;
            else
                string2_ = [string2_ string2(jj)];
            end;
        end;
    end;
    %i
    if (strcmp(string1_,string2_) == 1)
        %if (SNP_probeset(i).probe_chromosome == 8)
        %    fprintf(resultsFileID, '\n*R*');
        %else
        %    fprintf(resultsFileID, ['\n*' num2str(SNP_probeset(i).probe_chromosome) '*']);
        %end;
        %fprintf(resultsFileID, ['\t' num2str(SNP_probeset(i).probe_location) '\t']);
        %fprintf(resultsFileID, [string1_ '\t' string2_ '\t' num2str(Offset)]);
        % 
        %fprintf(resultsFileID, '\t4');
        %count_4 = count_4+1;
    else
        if (SNP_probeset(i).probe_polarity == 0)
            if (SNP_probeset(i).probe_chromosome == 8)
                fprintf(resultsFileID, '\nR');
            else
                fprintf(resultsFileID, ['\n' num2str(SNP_probeset(i).probe_chromosome)]);
            end;
            fprintf(resultsFileID, ['\t' num2str(SNP_probeset(i).probe_location) '\t']);
            fprintf(resultsFileID, [string1_ '\t' string2_ '\t0']);
            count_NaN = 0;
            for j = 1:length(SNP_probeset(i).probe_ch1)
                if (isnan(SNP_probeset(i).probe_ch1(j)) == 1)
                    count_NaN = count_NaN+1;
                end;
            end;
            if (count_NaN == 0)
                count_0b = count_0b+1;
            else
                count_0a = count_0a+1;
            end;
        elseif (SNP_probeset(i).probe_polarity == 1)
            if (SNP_probeset(i).probe_chromosome == 8)
                fprintf(resultsFileID, '\nR');
            else
                fprintf(resultsFileID, ['\n' num2str(SNP_probeset(i).probe_chromosome)]);
            end;
            fprintf(resultsFileID, ['\t' num2str(SNP_probeset(i).probe_location) '\t']);
            fprintf(resultsFileID, [string1_ '\t' string2_ '\t1']);
            count_1 = count_1+1;
        elseif (SNP_probeset(i).probe_polarity == 2)
            if (SNP_probeset(i).probe_chromosome == 8)
                fprintf(resultsFileID, '\nR');
            else
                fprintf(resultsFileID, ['\n' num2str(SNP_probeset(i).probe_chromosome)]);
            end;
            fprintf(resultsFileID, ['\t' num2str(SNP_probeset(i).probe_location) '\t']);
            fprintf(resultsFileID, [string2_ '\t' string1_ '\t2']);
            count_2 = count_2+1;
        end;
    end;
end;
fclose(resultsFileID);

fprintf(['\nUnnasigned, no data.          : ' num2str(count_0a) '\n']);
fprintf(['Unnasigned, never homozygous. : ' num2str(count_0b) '\n']);
fprintf(['Assigned. (1st=a, 2nd=b)      : ' num2str(count_1) '\n']);
fprintf(['Assigned. (1st=b, 2nd=a)      : ' num2str(count_2) '\n']);
fprintf(['probes identical.             : ' num2str(count_4) '\n']);
% cleanup.
%clear d i j s str v;
% delete any hanging waitbar dialogs or already saved figures.
set(0,'ShowHiddenHandles','on');
delete(get(0,'Children'));