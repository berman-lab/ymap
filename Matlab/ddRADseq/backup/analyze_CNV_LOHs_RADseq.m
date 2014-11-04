function analyze_CNV_LOHs_RADseq(projectName1_parent,projectName2_child,file1_parent_childLOH,file2_childCNV, ploidyString)
%    proj_name1=$1    parent
%    proj_name2=$2    child
%    proj1=$file1_fix : LOH from parent to child.
%    proj2=$file2_fix : child CNV.
%    ploidy=$ploidy

% log file start, for in-process analysis.
if (strcmp(projectName1_parent,projectName1_parent) == 1)
	fprintf(['ProjectName : [[[' projectName1_parent '[[[\n']);
else
	fprintf(['ProjectName : [[[' projectName1_parent '->' projectName2_child '[[[\n']);
end;

workingDir      = '/home/bermanj/shared/links/';
figureDir       = '~/';
CNV_verString   = 'v1';
INDEL_verString = 'v1';
SNP_verString   = 'v4';
LOH_verString   = 'v2';
displayBREAKS   = true;

fprintf('\n');
fprintf(['[analyze_CNV_LOHs.m data file inputs:\n']);
fprintf(['    ProjectName1_parent           : ' projectName1_parent       '\n']);
fprintf(['    ProjectName2_child            : ' projectName2_child        '\n']);
fprintf(['    LOH from parent to child file : ' file1_parent_childLOH     '\n']);
fprintf(['    child CNV file                : ' file2_childCNV            '\n']);
fprintf(['    Ploidy                        : ' ploidyString              '\n']);
fprintf('\n');

if (strcmp(projectName1_parent,strrep(projectName1_parent,' ','')) == 0) || ...
	(strcmp(projectName1_parent,strrep(projectName1_parent,'/','')) == 0) || ...
	(strcmp(projectName1_parent,strrep(projectName1_parent,'\','')) == 0)
    fprintf(['ERROR:[analyze_CNV_LOHs_RADseq_1.m]: <Reference Name> ''' projectName1_parent ''' cannot include spaces or slash characters.\n']);
    exit 1;
end;
if (strcmp(projectName2_child,strrep(projectName2_child,' ','')) == 0) || ... 
	(strcmp(projectName2_child,strrep(projectName2_child,'/','')) == 0) || ...
	(strcmp(projectName2_child,strrep(projectName2_child,'\','')) == 0)
    fprintf(['ERROR:[analyze_CNV_LOHs_RADseq_1.m]: <Project Name> ''' projectName2_child ''' cannot include spaces or slash characters.\n']);
    exit 1;
end;

if (strcmp(projectName1_parent,projectName2_child) == 0)
    %% parent and child are different strains.
    if (exist(file1_parent_childLOH,'file') == 0)
	fprintf(['ERROR:[analyze_CNV_LOHs_RADseq_1.m]: LOH data comparing parent project ''' projectName1_parent ''' to child project ''' projectName2_child ''' is not found.\n']);
	exit 1;
    end;
    if (exist(file2_childCNV,'file') == 0)
	fprintf(['ERROR:[analyze_CNV_LOHs_RADseq_1.m]: (child CNV) ''' projectName2_child ''' is not found.\n']);
	exit 1;
    end;
else
    %% parent and child are the same strain.
    if (exist(file1_parent_childLOH,'file') == 0)
	fprintf(['ERROR:[analyze_CNV_LOHs.m]: Project ''' projectName1 ''' SNP data is not found.\n']);
	exit 1;
    end;
    if (exist(file2_childCNV,'file') == 0)
	fprintf(['ERROR:[analyze_CNV_LOHs.m]: Project ''' projectName1 ''' CNV data is not found.\n']);
	exit 1;
    end;
end;

%% Grab the first column of the first line of the putative_CNV pileup file.
datafile            = [workingDir 'pileup_dir/' projectName2_child '_putative_CNVs_' CNV_verString '.txt'];
refdatafile         = [workingDir 'pileup_dir/' projectName1_parent '_putative_CNVs_' CNV_verString '.txt'];
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
    fprintf(['ERROR:[analyze_CNVs_RADseq_1.m]: The project ''' projectName2_child ''' and reference ''' projectName1_parent ''' appear to be for different genomes.\n']);
    exit;
end;

%% Generate CNV map from read count across genome.
ave_copy_num = 34; % filler number.
CNV_SNP_normalized_v1(projectName1_parent,projectName2_child,genome_data,ploidyString, ...
	ave_copy_num,CNV_verString,SNP_verString,LOH_verString,workingDir,figureDir,displayBREAKS);

fprintf('*--- End of ''analyze_CNV_LOHs_RADseq_1.m'' was reached ---*\n');
% log file end, for in-process analysis.
if (strcmp(projectName1_parent,projectName1_parent) == 1)
	fprintf(['ProjectName : ]]]' projectName1_parent ']]]\n']);
else
	fprintf(['ProjectName : ]]]' projectName1_parent '->' projectName2_child ']]]\n']);
end;

end
