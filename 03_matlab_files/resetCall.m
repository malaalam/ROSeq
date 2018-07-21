function []=resetCall(source, eventdata, varargin)
% *************************************************************************
% This function "resetCall" define the events that ensue when the user
% decides to reset the Pre-Processing Panel.
%
% Last update: 27 July, 2017
% *************************************************************************
% Define a structure that contains handles to all the objects in
% figureHandle
% *************************************************************************
% *************************************************************************
myhandles=guidata(varargin{1});

% Default Properties
  myhandles.projectFile='';
  myhandles.projectPath='';
  myhandles.results='';
  myhandles.defaultFileLocation='';
  myhandles.rawFileLocation='';
  myhandles.rawData='';
  myhandles.filterThreshold='';
  myhandles.normalizedFileLocation='';
  myhandles.normalizedData='';
  myhandles.numberGenes='';
  myhandles.groupOneStart='';
  myhandles.groupOneEnd='';
  myhandles.groupTwoStart='';
  myhandles.groupTwoEnd='';
% Set the fields equal to default
set(myhandles.preProcessing(1), 'Value', 0);
set(myhandles.preProcessing(2), 'Value', 0);
set(myhandles.preProcessing(3), 'Visible', 'off');
set(myhandles.preProcessing(4), 'Visible', 'off');
set(myhandles.preProcessing(5), 'Visible', 'off');
set(myhandles.preProcessing(6), 'Visible', 'off');
set(myhandles.preProcessing(7), 'Visible', 'off');
set(myhandles.preProcessing(8), 'Visible', 'off');
set(myhandles.preProcessing(8), 'String', '6');  
set(myhandles.preProcessing(10), 'String', 'Start (Gr 1)');
set(myhandles.preProcessing(11), 'String', 'End (Gr 1)');
set(myhandles.preProcessing(12), 'String', 'Start (Gr 2)');
set(myhandles.preProcessing(13), 'String', 'End (Gr 2)');
set(myhandles.preProcessing(14), 'Enable', 'off');
set(myhandles.preProcessing(15), 'Enable', 'off');


myhandles.statusString=[myhandles.statusString newline ...
        'Pre-processing Panel reset to default successfully.' ...
        newline];
    set(myhandles.status(1), 'String', myhandles.statusString, ...
        'ForegroundColor', [0 0.6 0]);
    jhEdit=findjobj(myhandles.status(1));
    jEdit=jhEdit.getComponent(0).getComponent(0);
    jEdit.setCaretPosition(jEdit.getDocument.getLength);
    % *************************************************************************
    % Save the structure/data 'myhandles' in the 'figureHandle - varargin{1}'
    % *************************************************************************
    guidata(varargin{1}, myhandles);
 
guidata(source, myhandles);