function [] = ploidyDetection(main_dir, user, genomeUser, project, genome, ploidyEstimateString)
%    Detects ploidy for each segment from the ChARM script and generates
%    'chr_segments.txt' file
%    <Project Name>    : the name of the project.
%    <Project File>    : the name of the project file, including directory location.
%    <Ploidy Estimate> : a numerical value for ploidy from other work.
    addpath('../');
    projectDir = [main_dir 'users/' user '/projects/' project '/'];
    genomeDir  = [main_dir 'users/' genomeUser '/genomes/' genome '/'];
    
    ploidy = str2num(ploidyEstimateString);
    
    [centromeres, chr_sizes, figure_details, annotations, ploidy_default] = Load_genome_information(genomeDir);
    [Aneuploidy]                                                          = Load_dataset_information(projectDir);
    
    % getting used chars
    for i = 1:length(figure_details)
        if (figure_details(i).chr == 0)
            if (strcmp(figure_details(i).label,'Key') == 1)
                continue
            end;
        else
            chr_names(figure_details(i).chr) = cellstr(figure_details(i).name);
            chr_in_use(figure_details(i).chr) = str2num(figure_details(i).useChr);
        end;
    end;
    num_chrs                          = length(chr_sizes);
    
    fprintf('\t|\tLoading "Common_CNV" data file, to be used in copy number estimation.\n');
    load([projectDir 'Common_CNV.mat']);   % 'CNVplot2', 'genome_CNV'
    % Ignoring ploidyAdjust since it's being used only for histogram
    % drawing
    [chr_breaks, chrCopyNum, ploidyAdjust] = FindChrSizes_4(Aneuploidy,CNVplot2,ploidy,num_chrs,chr_in_use);
    
    for chr = 1:num_chrs
        if (chr_in_use(chr) == 1)
            if (length(chrCopyNum{chr}) > 1)  % more than one segment, so lets examine if adjacent segments have different copyNums.
                %% Merge any adjacent segments with the same copy number.
                % add break representing left end of chromosome.
                breakCount_new         = 1;
                chr_breaks_new{chr}    = [];
                chrCopyNum_new{chr}    = [];
                chr_breaks_new{chr}(1) = 0.0;
                if (length(chrCopyNum{chr}) > 0)
                    chrCopyNum_new{chr}(1) = chrCopyNum{chr}(1);
                    for segment = 1:(length(chrCopyNum{chr})-1)
                        if (round(chrCopyNum{chr}(segment)) == round(chrCopyNum{chr}(segment+1)))
                            % two adjacent segments have identical copyNum and should be fused into one; don't add boundry to new list.
                        else
                            % two adjacent segments have different copyNum; add boundry to new list.
                            breakCount_new                      = breakCount_new + 1;
                            chr_breaks_new{chr}(breakCount_new) = chr_breaks{chr}(segment+1);
                            chrCopyNum_new{chr}(breakCount_new) = chrCopyNum{chr}(segment+1);
                        end;
                    end;
                end;
                % add break representing right end of chromosome.
                breakCount_new = breakCount_new+1;
                chr_breaks_new{chr}(breakCount_new) = 1.0;		
                % copy new lists to old.
                chr_breaks{chr} = chr_breaks_new{chr};
                chrCopyNum{chr} = [];
                chrCopyNum{chr} = chrCopyNum_new{chr};
            end;
        end;
    end;
    
    % creating chr_segments file
    segement_file = fopen([projectDir 'chr_segments.txt'], 'wt' );
    % header
    fprintf(segement_file, '#chr\tstartBP\tendBP\tploidy\n');
    for chr = 1:num_chrs
        % avoid running over chromosomes with empty copy number
        if (chr_in_use(chr) == 1 && ~isempty(chrCopyNum{chr}))
            for i = 1:length(chrCopyNum{chr})
                startBP = round(chr_breaks{chr}(i) * chr_sizes(chr).size);
                if (startBP == 0)
                    startBP = 1;
                end
                endBP = round(chr_breaks{chr}(i+1) * chr_sizes(chr).size);
                name = char(chr_names(chr));
                % write line
                fprintf(segement_file, '%s\t%d\t%d\t%d\n', name, startBP, endBP, round(chrCopyNum{chr}(i)));
            end
        end
    end
    fclose(segement_file);
    
    fprintf('\n\n#===========================#\n');
    fprintf(    '|END OF "ploidyDetection.m" script.|\n');
    fprintf(    '#===========================#\n');
end