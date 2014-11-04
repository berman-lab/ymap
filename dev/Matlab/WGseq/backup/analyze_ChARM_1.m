function  analyze_ChARM_1(projectName, projectFile)
% analyze_CNVS(<Project Name>,<Project File>)
% A componant of the sequence analysis pipeline, analyzing CNVs only.
%    <Project Name>    : the name of the project.
%    <Project File>    : the name of the project file, including directory location.

% log file start, for in-process analysis.
fprintf(['ProjectName : [[[' projectName '[[[\n']);

workingDir             = '/home/bermanj/shared/links/';
figureDir              = '~/';
CNV_verString          = 'v1';
displayBREAKS          = true;

%% ==========================================================================================================================

%% Generate ChARM map from previously analyzed CNV data mapped across the genome.
ChARM_v4(projectName, workingDir,figureDir);

fprintf('*--- End of ''analyze_ChARM.m'' was reached ---*\n');
%log file end, for in-process analysis
fprintf(['ProjectName : ]]]' projectName ']]]\n']);
end
