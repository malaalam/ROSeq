function []=saveResultsCall(source, eventdata, varargin)
% *************************************************************************
% This function "saveResults" define the events that ensue when the user
% decides to save the results from the differential expresssion analysis,
% as a separate csv file. The three variables (pVal, pAdjusted, log2FC )
% are exported as three columns in a csv file.
%
%
% Last update: 27 July, 2017
% *************************************************************************
% Define a structure that contains handles to all the objects in
% figureHandle
% *************************************************************************

myhandles=guidata(varargin{1});
try
    [file,path] = uiputfile('*.csv','Results');
    if isempty(file) || ~ischar(file)
        myhandles.statusString = [myhandles.statusString newline ...
            'No file location selected. Try again?' newline];
        set(myhandles.status(1), 'String', myhandles.statusString, 'ForegroundColor', [1 0 0]);
        jhEdit=findjobj(myhandles.status(1));
        jEdit=jhEdit.getComponent(0).getComponent(0);
        jEdit.setCaretPosition(jEdit.getDocument.getLength);
        guidata(source, myhandles);
    else
        geneNames=myhandles.results.geneNames;
        pVal=myhandles.results.pVal;
        pAdj=myhandles.results.padj;
        log2FC=myhandles.results.log2FC;
        T = table(geneNames,pVal,pAdj, log2FC);
        filename = [path file];
        writetable(T,filename)
        myhandles.statusString=[myhandles.statusString newline ...
            'Results saved in a *.csv file format at the desired location.' ...
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
catch e
    myhandles.statusString = [myhandles.statusString newline ...
        'The identifier was: ' e.identifier newline...
        'There was an error! The message was: ' e.message newline ...
        'Check "openCall.m" for further details.' newline];
    set(myhandles.status(1), 'String', myhandles.statusString);
    jhEdit=findjobj(myhandles.status(1));
    jEdit=jhEdit.getComponent(0).getComponent(0);
    jEdit.setCaretPosition(jEdit.getDocument.getLength);
    guidata(source, myhandles);
end

