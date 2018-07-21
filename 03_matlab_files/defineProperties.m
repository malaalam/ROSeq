function [ myhandles ] = defineProperties( source, eventdata, varargin )
% *************************************************************************
% This function 'defineProperties' defines the default properties of the
% ROSeq software
%
% *************************************************************************
% Define a structure that contains handles to all the objects in
% figureHandle
% *************************************************************************
  
  myhandles = guidata(source);
  myhandles.notesString= ...
      [newline '<Space for making some notes about the current project>'];
  myhandles.statusString=''; 
  myhandles.projectFile='';
  myhandles.projectPath='';
  myhandles.results='';

  
  myhandles.defaultFileLocation='';
  myhandles.rawFileLocation='';
  myhandles.rawData='';
  myhandles.filterThreshold=6;
  myhandles.normalizedFileLocation='';
  myhandles.normalizedData='';
  myhandles.numberGenes=''; % Number of Genes specified by the user
  myhandles.geneCount='';   % Number of Genes left after the filtering process
  myhandles.groupOneStart='';
  myhandles.groupOneEnd='';
  myhandles.groupTwoStart='';
  myhandles.groupTwoEnd='';
% *************************************************************************
% Save the structure
% *************************************************************************
guidata(source, myhandles)