function [] = ploidyDetection(mainDir, user, genomeUser, project, genome, ploidy)
%    Detects ploidy for each segment from the ChARM script and generates
%    'chr_segments.txt' file
%    <Project Name>    : the name of the project.
%    <Project File>    : the name of the project file, including directory location.
%    <Ploidy Estimate> : a numerical value for ploidy from other work.
    addpath('../');
    projectDir = fullfile(mainDir, 'users/', user, '/projects/', project, '/');
    genomeDir  = fullfile(mainDir, 'users/', genomeUser, '/genomes/', genome , '/');
    
    [~, chrSizes, figureDetails, ~, ~] = Load_genome_information(genomeDir);
    [Aneuploidy] = Load_dataset_information(projectDir);
    
    % preallocating memory
    chrNames = cell(1, length(figureDetails));
    chrInUse = zeros(1, length(figureDetails));
    
    % getting used chrs
    for i = 1:length(figureDetails)
        if (figureDetails(i).chr == 0 && strcmp(figureDetails(i).label,'Key') == 1)
            continue;
        else
            chrNames(figureDetails(i).chr) = cellstr(figureDetails(i).name);
            chrInUse(figureDetails(i).chr) = str2double(figureDetails(i).useChr);
        end;
    end;
    numChrs = length(chrSizes);
    
    fprintf('\t|\tLoading "Common_CNV" data file, to be used in copy number estimation.\n');
    load([projectDir 'Common_CNV.mat']);   % 'CNVplot2', 'genome_CNV'
    % Ignoring ploidyAdjust since it's being used only for histogram
    % drawing
    [chrBreaks, chrCopyNum, ~] = FindChrSizes_4(Aneuploidy, CNVplot2, ploidy, numChrs, chrInUse);
    
    for chr = 1:numChrs
        % more than one segment, so lets examine if adjacent segments have different copyNums.
        if (chrInUse(chr) == 1 && length(chrCopyNum{chr}) > 1)
            %% Merge any adjacent segments with the same copy number.
            % add break representing left end of chromosome.
            breakCountNew = 1;
            chrBreaksNew{chr} = [];
            chrCopyNumNew{chr} = [];
            chrBreaksNew{chr}(1) = 0.0;
            chrCopyNumNew{chr}(1) = chrCopyNum{chr}(1);
            for segment = 1:(length(chrCopyNum{chr})-1)
                if (round(chrCopyNum{chr}(segment)) == round(chrCopyNum{chr}(segment+1)))
                    % two adjacent segments have identical copyNum and should be fused into one; don't add boundry to new list.
                else
                    % two adjacent segments have different copyNum; add boundry to new list.
                    breakCountNew = breakCountNew + 1;
                    chrBreaksNew{chr}(breakCountNew) = chrBreaks{chr}(segment+1);
                    chrCopyNumNew{chr}(breakCountNew) = chrCopyNum{chr}(segment+1);
                end;
            end;
            % add break representing right end of chromosome.
            breakCountNew = breakCountNew+1;
            chrBreaksNew{chr}(breakCountNew) = 1.0;		
            % copy new lists to old.
            chrBreaks{chr} = chrBreaksNew{chr};
            chrCopyNum{chr} = [];
            chrCopyNum{chr} = chrCopyNumNew{chr};  
        end;
    end;
    
    % creating chr_segments file
    segementFile = fopen([projectDir 'chr_segments.txt'], 'wt' );
    % header
    fprintf(segementFile, '#chr\tstartBP\tendBP\tploidy\n');
    for chr = 1:numChrs
        % avoid running over chromosomes with empty copy number
        if (chrInUse(chr) == 1 && ~isempty(chrCopyNum{chr}))
            for i = 1:length(chrCopyNum{chr})
                startBP = round(chrBreaks{chr}(i) * chrSizes(chr).size);
                if (startBP == 0)
                    startBP = 1;
                end
                endBP = round(chrBreaks{chr}(i+1) * chrSizes(chr).size);
                name = char(chrNames(chr));
                % write line
                fprintf(segementFile, '%s\t%d\t%d\t%d\n', name, startBP, endBP, round(chrCopyNum{chr}(i)));
            end
        end
    end
    fclose(segementFile);
    
    fprintf('\n\n#===========================#\n');
    fprintf(    '|END OF "ploidyDetection.m" script.|\n');
    fprintf(    '#===========================#\n');
end