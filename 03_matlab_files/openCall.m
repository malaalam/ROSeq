function []=openCall(source, eventdata, varargin)
% *************************************************************************
% This function "openCall" define the events that ensue when the user
% decides to open an already existing session of ROSeq.
%
% *************************************************************************
% Define a structure that contains handles to all the objects in
% figureHandle
% *************************************************************************
myhandles=guidata(varargin{1});
try
[fileName,pathName] = uigetfile('*.mat','Select the Project File to Open (*.mat)');
if isempty(fileName) || ~ischar(fileName)
myhandles.statusString = [myhandles.statusString newline ...
'Project File not found. Try again?' newline];
set(myhandles.status(1), 'String', myhandles.statusString, 'ForegroundColor', [1 0 0]);
jhEdit=findjobj(myhandles.status(1));
jEdit=jhEdit.getComponent(0).getComponent(0);
jEdit.setCaretPosition(jEdit.getDocument.getLength);
guidata(source, myhandles);
else
close(gcf);
structure=load([pathName fileName]);
myhandles=getfield(structure, 'myhandles');
myhandles.statusString=[myhandles.statusString newline 'Project File "' fileName '" successfully loaded.' newline];
set(myhandles.status(1), 'String', myhandles.statusString, 'ForegroundColor', [0 0.6 0]);
jhEdit=findjobj(myhandles.status(1));
jEdit=jhEdit.getComponent(0).getComponent(0);
jEdit.setCaretPosition(jEdit.getDocument.getLength);
guidata(gcf, myhandles);
 end
catch e
    myhandles.statusString = [myhandles.statusString newline ...
'The identifier was: ' e.identifier newline...
'There was an error! The message was: ' e.message newline ...
'Check "openCall.m" for further details.' newline];
set(myhandles.status(1), 'String', myhandles.statusString);
jhEdit=findjobj(myhandles.status(1));
jEdit=jhEdit.getComponent(0).getComponent(0);
jEdit.setCaretPosition(jEdit.getDocument.getLength);

% *************************************************************************
% Save the structure/data 'myhandles' in the 'figureHandle - varargin{1}'
% *************************************************************************
guidata(source, myhandles);
end


end