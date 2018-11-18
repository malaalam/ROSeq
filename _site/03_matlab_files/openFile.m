function []=openFile(source, eventdata, varargin)
% *************************************************************************
% This function "openFile" define the events that ensue when the user
% decides to open a normalized read count data file.
%
% Last update: 27 July, 2017
% *************************************************************************
% Define a structure that contains handles to all the objects in
% figureHandle
% *************************************************************************
% *************************************************************************
myhandles=guidata(varargin{1});
% *************************************************************************
% Break the complete name of the file into two parts -
% the long Path name and the immediate fileName
% *************************************************************************
if source ==myhandles.preProcessing(3)
[fileName, pathName] = ...
    uigetfile({'*.csv'}, ...
    'Select the Read Count csv File', myhandles.defaultFileLocation);

% *****************************************************************
% Add an error message in the Status Panel. The error message
% concerns the failure to load the file.

% Use the function "findjobj" to find the underlying java handle to
% the "Status Panel". Use the function "setCaretPosition" to scroll
% the "Status Panel" to the bottom
% *****************************************************************
if isempty(fileName) || ~ischar(fileName)
    myhandles.statusString=[myhandles.statusString newline ...
        'Read Count File upload not successful. Try Again?' ...
        newline];
    set(myhandles.status(1), 'String', myhandles.statusString, ...
        'ForegroundColor', [1 0 0]);
    jhEdit=findjobj(myhandles.status(1));
    jEdit=jhEdit.getComponent(0).getComponent(0);
    jEdit.setCaretPosition(jEdit.getDocument.getLength);
    % *************************************************************************
    % Save the structure/data 'myhandles' in the 'figureHandle - varargin{1}'
    % *************************************************************************
    guidata(varargin{1}, myhandles);
else
    myhandles.defaultFileLocation=[pathName fileName];
    if get(myhandles.preProcessing(2), 'Value') ==1
        myhandles.normalizedFileLocation=[pathName fileName];
        myhandles.normalizedData=uiimport([pathName fileName]);
        set(myhandles.preProcessing(5), 'Visible', 'on');
        set(myhandles.preProcessing(6), 'Visible', 'on');
        set(myhandles.preProcessing(6), 'String', num2str(size(myhandles.normalizedData.data,1)));
        myhandles.numberGenes=size(myhandles.normalizedData.data,1);
        myhandles.rawData='';
        set(myhandles.preProcessing(7), 'Visible', 'off');
        set(myhandles.preProcessing(8), 'Visible', 'off');
        
    else
        myhandles.rawFileLocation=[pathName fileName];
        myhandles.rawData=uiimport([pathName fileName]);
        set(myhandles.preProcessing(5), 'Visible', 'on');
        set(myhandles.preProcessing(6), 'Visible', 'on');
        set(myhandles.preProcessing(6), 'String', num2str(size(myhandles.rawData.data,1)));
        myhandles.numberGenes=size(myhandles.rawData.data,1);
        myhandles.normalizedData='';
        set(myhandles.preProcessing(7), 'Visible', 'on');
        set(myhandles.preProcessing(8), 'Visible', 'on');
     end
    set(myhandles.preProcessing(3), 'Visible', 'off');
    set(myhandles.preProcessing(4), ...
        'String', ['File loaded from - ' pathName fileName], ...
        'Visible', 'on');
   

    % *****************************************************************
    % Add a message indicating success in the Status Panel.
    % The success message concerns the fact that the file name was
    % successfully processed.
    
    % Use the function "findjobj" to find the underlying java handle to
    % the "Status Panel". Use the function "setCaretPosition" to scroll
    % the "Status Panel" to the bottom
    % *****************************************************************
    myhandles.statusString=[myhandles.statusString newline ...
        'Read Count File successfully loaded.' newline];
    set(myhandles.status(1), 'String', myhandles.statusString, ...
        'ForegroundColor', [0 0.6 0]);
    jhEdit=findjobj(myhandles.status(1));
    jEdit=jhEdit.getComponent(0).getComponent(0);
    jEdit.setCaretPosition(jEdit.getDocument.getLength);
    % *****************************************************************
    % Display the File name in the ``Build Element Panel"
    % *****************************************************************
    guidata(varargin{1}, myhandles);
end
elseif source == myhandles.postProcessing2(2)
    [fileName, pathName] = ...
    uigetfile({'*.txt'}, ...
    'Select the Read Count csv File', myhandles.defaultFileLocation);
    myhandles.TP_path=pathName;
    myhandles.TP_file=fileName;
    set(myhandles.postProcessing2(2), 'Visible', 'off');
    set(myhandles.postProcessing2(5), 'Visible', 'on', 'string', fileName);

elseif source == myhandles.postProcessing2(4)
    [fileName, pathName] = ...
    uigetfile({'*.txt'}, ...
    'Select the Read Count csv File', myhandles.defaultFileLocation);
    myhandles.TN_path=pathName;
    myhandles.TN_file=fileName;
    set(myhandles.postProcessing2(4), 'Visible', 'off');
    set(myhandles.postProcessing2(6), 'Visible', 'on', 'string', fileName);
end
% *************************************************************************
% Save the structure/data 'myhandles' in the 'figureHandle - varargin{1}'
% *************************************************************************
guidata(source, myhandles);
end


