function analyze_CNVs_RADseq_2(projectName, referenceName, ploidyEstimateString,ploidyBaseString)
% analyze_CNVS_RADseq_1(<Project Name>,<Project File>,<Reference Name>,<Reference File>,<Ploidy Estimate)
% A componant of the sequence analysis pipeline, analyzing CNVs only.
%    <Project Name>    : the name of the project.
%    <Project File>    : the name of the project file, including directory location.
%    <Reference Name>  : the name of the project used as a normalization reference.
%    <Reference File>  : the name of the reference project file, including directory location.
%    <Ploidy Estimate> : a numerical value for ploidy from other work.

% log file start, for in-process analysis.
fprintf(['ProjectName : [[[' projectName '[[[\n']);

workingDir             = '/home/bermanj/shared/links/';
figureDir              = '~/';
CNV_verString          = 'v1';
INDEL_verString        = 'v1';
SNP_verString          = 'v4';
rDNA_verString         = 'v1';
displayBREAKS          = true;
referenceCHR           = 1;

if (strcmp(projectName,strrep(projectName,' ','')) == 0) || (strcmp(projectName,strrep(projectName,'/','')) == 0) || (strcmp(projectName,strrep(projectName,'\','')) == 0)
    fprintf(['ERROR:[analyze_CNVs_RADseq_1.m]: <Project Name> ''' projectName ''' cannot include spaces or slash characters.\n']);
    exit;
end;
if (strcmp(referenceName,strrep(referenceName,' ','')) == 0) || (strcmp(referenceName,strrep(referenceName,'/','')) == 0) || (strcmp(referenceName,strrep(referenceName,'\','')) == 0)
    fprintf(['ERROR:[analyze_CNVs_RADseq_1.m]: <Reference Name> ''' projectName ''' cannot include spaces or slash characters.\n']);
    exit;
end;

%% Grab the first column of the first line of the putative_CNV pileup file.
datafile            = [workingDir 'pileup_dir/' projectName '_putative_CNVs_' CNV_verString '.txt'];
refdatafile         = [workingDir 'pileup_dir/' referenceName '_putative_CNVs_' CNV_verString '.txt'];
data                = fopen(datafile);
ref                 = fopen(refdatafile);
lineData            = fgetl(data);           % grab the first line of the putative_CNV datafile.
lineRef             = fgetl(ref);
exampleChrName_data = sscanf(lineData, '%s',1);  % grab the first column, ex : 'ChrA_C_glabrata_CBS138';
exampleChrName_ref  = sscanf(lineRef, '%s',1);
fclose(data);
fclose(ref);

%% Compare the chromosome name found above to the registered genomes in order to determine which genome is being used in this project.
genomeList = dir('genomeSpecific/');
for i = length(genomeList):-1:1
    if (genomeList(i).isdir == false) || (strcmp(genomeList(i).name,'.') == true) || (strcmp(genomeList(i).name,'..') == true)
        genomeList(i) = [];
    end;
end;
for i = 1:length(genomeList)
    genomeTest{i} = genomeList(i).name;
end;
genome_data        = 'unknown';
genome_ref         = 'unknown';
for i = 1:length(genomeTest)
    [~,~,figure_details,~,~] = Load_genome_information_quiet_1(workingDir,figureDir,genomeTest{i});
    for j = 1:length(figure_details)
        if (strcmp(exampleChrName_data,figure_details(j).name) == 1)
            genome_data = genomeTest{i};
        end;
        if (strcmp(exampleChrName_ref,figure_details(j).name) == 1)
            genome_ref = genomeTest{i};
        end;
    end;
end;

if (strcmp(genome_data,genome_ref) == 0) 
    fprintf(['ERROR:[analyze_CNVs_RADseq_1.m]: The project ''' projectName ''' and reference ''' referenceName ''' appear to be for different genomes.\n']);
    exit;
end;

%% Generate CNV map from read count across genome.
% Process dataset by genome restriction fragments.

if (displayBREAKS == true)
    displayBREAKS_ = 'true';
else
    displayBREAKS_ = 'false';
end;
fprintf('Pre-processing data by restriction digested genome.\n');

%% Generate fragment_length_bias file for experimental dataset.
arguments  = [projectName ' ' genome_data ' "' workingDir '" "' figureDir '"'];
outputFile = [workingDir 'pileup_dir/' projectName '_RADseq_digest_analysis_CNV.txt'];
systemCall = ['python genome_process_for_RADseq.fragment_length_bias_1.py ' arguments ' > ' outputFile];
if (exist(outputFile,'file') == 0)
    fprintf('Processing project against resriction digested genome.\n');
    fprintf(['\tSystem call to Python from within Matlab: \n\t' systemCall '\n']);
    fprintf('Processing pileup file against restriction digested genome.');
    system(systemCall);
else
    fprintf('Pileup file for project already processed against restriction digested genome.\n');
end;

%% Generate fragment_length_bias file for reference dataset.
arguments  = [projectName ' ' genome_data ' "' workingDir '" "' figureDir '"'];
outputFile = [workingDir 'pileup_dir/' referenceName '_RADseq_digest_analysis_CNV.txt'];
systemCall = ['python genome_process_for_RADseq.fragment_length_bias_1.py ' arguments ' > ' outputFile];
if (exist(outputFile,'file') == 0)   
    fprintf('Processing reference against resriction digested genome.\n');
    fprintf(['\tSystem call to Python from within Matlab: \n\t' systemCall '\n']);
    fprintf('Processing pileup file against restriction digested genome.');
    system(systemCall);
else
    fprintf('Pileup file for reference already processed against restriction digested genome.\n');
end;

%% This method normalizes the experimental read data by the read data from SC5314.  [in 'old_scripts' directory.]
% CNV_v6_normalized(projectName,genome_data,referenceName,genome_ref,ploidyEstimateString, ...
%                  CNV_verString,rDNA_verString,workingDir,figureDir,displayBREAKS, referenceCHR);

%% This method corrects the read count bias seen in relation to the restriction fragment lengths.   [in 'old_scripts' directory.]
% CNV_v6_fragmentLengthCorrected_1(projectName,genome_data,referenceName,genome_ref,ploidyEstimateString, ...
%                                  CNV_verString,rDNA_verString,workingDir,figureDir,displayBREAKS, referenceCHR);

%% This method corrects the read count bias and normalizes against SC5314.
% CNV_v6_fragmentLengthCorrected_2(projectName,genome_data,referenceName,genome_ref,ploidyEstimateString,ploidyBaseString, ...
%                                  CNV_verString,rDNA_verString,workingDir,figureDir,displayBREAKS, referenceCHR);

%% This method corrects the read count bias and normalizes against SC5314.
% CNV_v6_fragmentLengthCorrected_3(projectName,genome_data,referenceName,genome_ref,ploidyEstimateString,ploidyBaseString, ...
%                                  CNV_verString,rDNA_verString,workingDir,figureDir,displayBREAKS, referenceCHR);

%% Same as above, with better annotation display.
CNV_v6_fragmentLengthCorrected_4(projectName,genome_data,referenceName,genome_ref,ploidyEstimateString,ploidyBaseString, ...
                                 CNV_verString,rDNA_verString,workingDir,figureDir,displayBREAKS, referenceCHR);

fprintf('*--- End of ''analyze_CNVs_RADseq_1.m'' was reached ---*\n');
% log file end, for in-process analysis.
fprintf(['ProjectName : ]]]' projectName ']]]\n']);
end
