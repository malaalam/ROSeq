function []=selectRadio(source, eventdata, varargin)
% *************************************************************************
% This function "selectRadio" define the events that ensue when the user
% chooses one of the two radio buttons for uploading raw or normalized
% read count data file.
%
% Last update: 27 July, 2017
% *************************************************************************
% Define a structure that contains handles to all the objects in
% figureHandle
% *************************************************************************
% *************************************************************************
myhandles=guidata(varargin{1});
if source==myhandles.preProcessing(1) && get(myhandles.preProcessing(1), 'Value')==1
    set(myhandles.preProcessing(2), 'Value', 0);
    set(myhandles.preProcessing(3), 'Visible', 'on');
    set(myhandles.preProcessing(4), 'Visible', 'off');
    myhandles.normalizedFileLocation='';
    myhandles.statusString=[myhandles.statusString newline ...
        'Raw Read Counts File radio button selected successfully.' ...
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
elseif source==myhandles.preProcessing(1) && get(myhandles.preProcessing(1), 'Value')==0
    set(myhandles.preProcessing(2), 'Value', 1);
    set(myhandles.preProcessing(3), 'Visible', 'on');
    set(myhandles.preProcessing(4), 'Visible', 'off');
    myhandles.rawFileLocation='';
    myhandles.statusString=[myhandles.statusString newline ...
        'Normalized Read Counts File radio button selected successfully.' ...
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
elseif source==myhandles.preProcessing(2) && get(myhandles.preProcessing(2), 'Value')==1
    set(myhandles.preProcessing(1), 'Value', 0);
    set(myhandles.preProcessing(3), 'Visible', 'on');
    set(myhandles.preProcessing(4), 'Visible', 'off');
    myhandles.rawFileLocation='';
    myhandles.statusString=[myhandles.statusString newline ...
        'Normalized Read Counts File radio button selected successfully.' ...
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
elseif source==myhandles.preProcessing(2) && get(myhandles.preProcessing(2), 'Value')==0
    set(myhandles.preProcessing(1), 'Value', 1);
    set(myhandles.preProcessing(3), 'Visible', 'on');
    set(myhandles.preProcessing(4), 'Visible', 'off');
    myhandles.normalizedFileLocation='';
    myhandles.statusString=[myhandles.statusString newline ...
        'Raw Read Counts File radio button selected successfully.' ...
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
end
set(myhandles.preProcessing(5), 'Visible', 'off');
set(myhandles.preProcessing(6), 'Visible', 'off');
set(myhandles.preProcessing(7), 'Visible', 'off');
set(myhandles.preProcessing(8), 'Visible', 'off');

guidata(source, myhandles);