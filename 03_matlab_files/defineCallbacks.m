function [myhandles] = defineCallbacks( source, eventdata, varargin )
% *************************************************************************
% This function 'defineCallbacks' specifies the callbacks i.e. what happens
% when the user clicks a button, selects an option in a checkbox et cetera
%
% *************************************************************************
% Define a structure that contains handles to all the objects in
% the figure window - source
% *************************************************************************
myhandles = guidata(source);
% *************************************************************************
% Define the first instruction to the user in the Status Panel. The instruction  
% concerns selecting one of the three acoustic elements (Source, Transfer 
% Matrix and Open End). Use the "findjobj" function to find the underlying 
% java handle to the panel. Use the "setCaretPosition" function to always 
% scroll the Status Panel to the bottom.
% *************************************************************************
myhandles.statusString=[myhandles.statusString newline ...
    'New Project loaded.' newline...
    'Choose a *.csv file containing raw read counts or normalized read counts' ...
    ' for uploading.' newline];
jhEdit=findjobj(myhandles.status(1));
jEdit=jhEdit.getComponent(0).getComponent(0);
jEdit.setCaretPosition(jEdit.getDocument.getLength);
% *************************************************************************
% Save the structure/data 'myhandles' in the 'figureHandle - source'
% *************************************************************************
guidata(source, myhandles)
% *************************************************************************
% Define the callbacks to the various uicontrol features. In order to
% look up these features (myhandles.status(1), myhandles.About et cetera) 
% check the function named "definePanels"
% *************************************************************************
set(myhandles.status(1), 'String', myhandles.statusString);
set(myhandles.SaveAs_NewProject, 'Callback', {@saveAsCall, source});
set(myhandles.Save_NewProject, 'Callback', {@saveCall, source});
set(myhandles.File_NewProject, 'Callback', {@initializeGUI, source});
set(myhandles.Open_Project, 'Callback', {@openCall, source});
set(myhandles.HelpContents, 'Callback', {@helpCall, source});
set(myhandles.About, 'Callback', {@aboutCall, source});
set(myhandles.preProcessing(1), 'Callback', {@selectRadio, source});
set(myhandles.preProcessing(2), 'Callback', {@selectRadio, source});
set(myhandles.preProcessing(3), 'Callback', {@openFile, source});
set(myhandles.preProcessing(6), 'Callback', {@inputNumeric, source});
set(myhandles.preProcessing(8), 'Callback', {@inputNumeric, source});
set(myhandles.preProcessing(10), 'Callback', {@inputNumeric, source});
set(myhandles.preProcessing(11), 'Callback', {@inputNumeric, source});
set(myhandles.preProcessing(12), 'Callback', {@inputNumeric, source});
set(myhandles.preProcessing(13), 'Callback', {@inputNumeric, source});
set(myhandles.preProcessing(14), 'Callback', {@callMainFile, source});
set(myhandles.preProcessing(15), 'Callback', {@saveResultsCall, source});
set(myhandles.preProcessing(16), 'Callback', {@resetCall, source});
% *************************************************************************
% Save the structure/data 'myhandles' in the 'figureHandle - source'
% *************************************************************************
guidata(source, myhandles)
end

