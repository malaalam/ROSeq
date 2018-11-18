function []=saveCall(source, eventdata, varargin)
% *************************************************************************
% This function "saveCall" define the events that ensue when the user
% decides to save the current project.

% *************************************************************************
% Define a structure that contains handles to all the objects in
% figureHandle
% *************************************************************************

myhandles=guidata(varargin{1});
% save('newProject.mat', 'myhandles') ;
% *************************************************************************
% If the project hasn't been saved before, it goes in the Save As mode
% *************************************************************************
if(isempty(myhandles.projectFile))
    [file,path] = uiputfile('newProject.mat','Save As Project');
    myhandles.projectPath=path;
    myhandles.projectFile=file;
    if (ischar(myhandles.projectFile) && ~isempty(myhandles.projectFile) )
        myhandles.statusString=[myhandles.statusString newline ...
            'Project File "' myhandles.projectFile '" successfully saved at :-' ...
            myhandles.projectPath newline];
        set(myhandles.status(1), 'String', myhandles.statusString, ...
            'ForegroundColor', [0 0.6 0]);
        jhEdit=findjobj(myhandles.status(1));
        jEdit=jhEdit.getComponent(0).getComponent(0);
        jEdit.setCaretPosition(jEdit.getDocument.getLength);
        set(varargin{1}, 'name', ['ROSeq: ' myhandles.projectFile(1:end-4)]);
        save([myhandles.projectPath myhandles.projectFile], 'myhandles');
        
    else
        myhandles.statusString=[myhandles.statusString newline 'Project File' ...
            ' was not saved. Try again?'  ...
            newline];
        set(myhandles.status(1), 'String', myhandles.statusString, ...
            'ForegroundColor', [1 0 0]);
        jhEdit=findjobj(myhandles.status(1));
        jEdit=jhEdit.getComponent(0).getComponent(0);
        jEdit.setCaretPosition(jEdit.getDocument.getLength);
    end
else
    if (ischar(myhandles.projectFile) && ~isempty(myhandles.projectFile) )
        myhandles.statusString=[myhandles.statusString newline 'Project File "' ...
            myhandles.projectFile '" successfully saved at :-' myhandles.projectPath ...
            newline];
        set(myhandles.status(1), 'String', myhandles.statusString, ...
            'ForegroundColor', [0 0.6 0]);
        jhEdit=findjobj(myhandles.status(1));
        jEdit=jhEdit.getComponent(0).getComponent(0);
        jEdit.setCaretPosition(jEdit.getDocument.getLength);
        set(varargin{1}, 'name', ['ROSeq: ' myhandles.projectFile(1:end-4)]);
        save([myhandles.projectPath myhandles.projectFile], 'myhandles');
       
    else
        myhandles.statusString=[myhandles.statusString newline 'Project File' ...
             ' was not saved. Try again?'  ...
            newline];
        set(myhandles.status(1), 'String', myhandles.statusString, ...
            'ForegroundColor', [1 0 0]);
        set(myhandles.confirmNetwork(1), 'Data', myhandles.stageData);
        jhEdit=findjobj(myhandles.status(1));
        jEdit=jhEdit.getComponent(0).getComponent(0);
        jEdit.setCaretPosition(jEdit.getDocument.getLength);
    end
end
% *************************************************************************
% Save the structure/data 'myhandles' in the 'figureHandle - varargin{1}'
% *************************************************************************
guidata(varargin{1}, myhandles);