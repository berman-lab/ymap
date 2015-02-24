function [CGH_data_all SNP_data_all] = Load_reference_data_1()

% Load currently used calibration data.
load('cal_data.mat');
Cal_data = SNP_probeset;
Cal_data_length = length(Cal_data);

% Choose directories.
dirs = uipickfiles;

% setup data stores.
SNP_data_all = [];
CGH_data_all = [];

for ii = 1:length(dirs)
    if ispc  % Windows
        file_dir = [char(dirs(ii)) '\'];
    else     % MacOS
        file_dir = [char(dirs(ii)) '/'];
    end;

    % Ensure data files exist.
    if (exist([file_dir 'SNP_data.mat'],'file') == 0) || (exist([file_dir 'CGH_data.mat'],'file') == 0)
        files_found         = BlueFuse_file_load_4(file_dir);
    end;
    
    % Load data files.
    load([file_dir 'SNP_data.mat']);
    load([file_dir 'CGH_data.mat']);

    if (ii == 1)
        for i = 1:2:Cal_data_length
            if (Cal_data(i).probe_polarity == 1) || (Cal_data(i).probe_polarity == 0)   % [no-swap]
                SNP_data_all(i).probe_ID           = Cal_data(i).probe_ID;
                SNP_data_all(i+1).probe_ID         = Cal_data(i).probe_ID;
                SNP_data_all(i).probe_chromosome   = Cal_data(i).probe_chromosome;
                SNP_data_all(i+1).probe_chromosome = Cal_data(i).probe_chromosome;
                SNP_data_all(i).probe_location     = Cal_data(i).probe_location;
                SNP_data_all(i+1).probe_location   = Cal_data(i).probe_location;
                SNP_data_all(i).probe_data         = Cal_data(i+1).probe_ch2;
                SNP_data_all(i+1).probe_data       = Cal_data(i).probe_ch2;
            elseif (Cal_data(i).probe_polarity == 2)   % [swap]
                SNP_data_all(i).probe_ID           = Cal_data(i).probe_ID;
                SNP_data_all(i+1).probe_ID         = Cal_data(i).probe_ID;
                SNP_data_all(i).probe_chromosome   = Cal_data(i).probe_chromosome;
                SNP_data_all(i+1).probe_chromosome = Cal_data(i).probe_chromosome;
                SNP_data_all(i).probe_location     = Cal_data(i).probe_location;
                SNP_data_all(i+1).probe_location   = Cal_data(i).probe_location;
                SNP_data_all(i).probe_data         = Cal_data(i).probe_ch2;
                SNP_data_all(i+1).probe_data       = Cal_data(i+1).probe_ch2;
            end;
        end;
        for i = 1:length(probeset2)
            CGH_data_all(i).probe_ID         = probeset2(i).probe_ID;
            CGH_data_all(i).probe_chromosome = probeset2(i).probe_chromosome;
            CGH_data_all(i).probe_location   = probeset2(i).probe_location;
            CGH_data_all(i).probe_data       = probeset2(i).probe_ch2;
        end;
    else
        for i = 1:2:Cal_data_length
            if (Cal_data(i).probe_polarity == 1) || (Cal_data(i).probe_polarity == 0)   % [no-swap]
                SNP_data_all(i).probe_data         = [SNP_data_all(i).probe_data   Cal_data(i+1).probe_ch2  ];
                SNP_data_all(i+1).probe_data       = [SNP_data_all(i+1).probe_data Cal_data(i).probe_ch2    ];
            elseif (Cal_data(i).probe_polarity == 2)   % [swap]
                SNP_data_all(i).probe_data         = [SNP_data_all(i).probe_data   Cal_data(i).probe_ch2    ];
                SNP_data_all(i+1).probe_data       = [SNP_data_all(i+1).probe_data Cal_data(i+1).probe_ch2  ];
            end;
        end;
        for i = 1:length(probeset2)
            CGH_data_all(i).probe_data = [CGH_data_all(i).probe_data probeset2(i).probe_ch2];
        end;
    end;
end;

for i = 1:Cal_data_length
    count = 1;
    for j = 1:length(SNP_data_all(i).probe_data)
        if (isnan(SNP_data_all(i).probe_data(j)) == 0)
            SNP_data_all(i).probe_data_trimmed(count) = SNP_data_all(i).probe_data(j);
            count = count+1;
        end;
    end;
end;
for i = 1:length(probeset2)
    count = 1;
    for j = 1:length(CGH_data_all(i).probe_data)
        if (isnan(CGH_data_all(i).probe_data(j)) == 0)
            CGH_data_all(i).probe_data_trimmed(count) = CGH_data_all(i).probe_data(j);
            count = count+1;
        end;
    end;
end;

for i = 1:Cal_data_length
    if (length(SNP_data_all(i).probe_data) > 0)
        SNP_data_all(i).probe_ave = sum(SNP_data_all(i).probe_data)/length(SNP_data_all(i).probe_data_trimmed);
    else
        SNP_data_all(i).probe_ave = NaN;
    end;
end;
for i = 1:length(probeset2)
    if (length(CGH_data_all(i).probe_data) > 0)
        CGH_data_all(i).probe_ave = sum(CGH_data_all(i).probe_data_trimmed)/length(CGH_data_all(i).probe_data_trimmed);
    else
        CGH_data_all(i).probe_ave = NaN;
    end;
end;

% normalize each SNP dataset to average of 1.
total(1:40) = 0;
counts(1:40) = 0;
for i = 1:length(SNP_data_all)
    for j = 1:length(SNP_data_all(i).probe_data)
        if (isnan(SNP_data_all(i).probe_data(j)) == 0)
            total(j) = total(j)+SNP_data_all(i).probe_data(j);
            counts(j) = counts(j)+1;
        end;
    end;
end;
for j = 1:length(total)
    total(j) = total(j)/counts(j);
end;
for i = 1:length(SNP_data_all)
    for j = 1:length(SNP_data_all(i).probe_data)
        SNP_data_all(i).probe_data(j) = SNP_data_all(i).probe_data(j)/total(j);
    end;
end;

% normalize each CGH dataset to average of 1.
total(1:40) = 0;
counts(1:40) = 0;
for i = 1:length(CGH_data_all)
    for j = 1:length(CGH_data_all(i).probe_data)
        if (isnan(CGH_data_all(i).probe_data(j)) == 0)
            total(j) = total(j)+CGH_data_all(i).probe_data(j);
            counts(j) = counts(j)+1;
        end;
    end;
end;
for j = 1:length(total)
    total(j) = total(j)/counts(j);
end;
for i = 1:length(CGH_data_all)
    for j = 1:length(CGH_data_all(i).probe_data)
        CGH_data_all(i).probe_data(j) = CGH_data_all(i).probe_data(j)/counts(j);
    end;
end;

end


%d = [];
%for i = 1:2:Cal_data_length
%    for j = 1:40;
%        if (isnan(SNP_data_all(i).probe_data(j)) == 0) && (isnan(SNP_data_all(i+1).probe_data(j)) == 0)
%            d = [d atan2(SNP_data_all(i+1).probe_data(j),SNP_data_all(i).probe_data(j))];
%        end;
%    end;
%end;
%plot(hist(d,100));