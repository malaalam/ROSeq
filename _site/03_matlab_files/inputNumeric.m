function []=inputNumeric(source, eventdata, varargin)
% *************************************************************************
% This function "inputNumeric" define the events that ensue when the user
% decides to enter a numeric input in a edit panel.
%
% Last update: 27 July, 2017
% *************************************************************************
% Define a structure that contains handles to all the objects in
% figureHandle
% *************************************************************************
% *************************************************************************
myhandles=guidata(varargin{1});
if source==myhandles.preProcessing(10)
    myhandles.groupOneStart=str2double(get(myhandles.preProcessing(10), 'string'));
    if(isnan(str2double(get(myhandles.preProcessing(10), 'string'))))
        set(myhandles.preProcessing(10), 'string', 'Start (Gr 1)');
        myhandles.groupOneStart='';
        myhandles.statusString=[myhandles.statusString newline ...
        'Please enter a numeric input.' newline];
        set(myhandles.status(1), 'String', myhandles.statusString, ...
        'ForegroundColor', [1 0 0]);
        jhEdit=findjobj(myhandles.status(1));
        jEdit=jhEdit.getComponent(0).getComponent(0);
        jEdit.setCaretPosition(jEdit.getDocument.getLength);
        guidata(varargin{1}, myhandles);
    elseif myhandles.groupOneStart>myhandles.groupOneEnd
        set(myhandles.preProcessing(10), 'string', 'Start (Gr 1)');
        myhandles.groupOneStart='';
        myhandles.statusString=[myhandles.statusString newline ...
        'The Starting Column Index should be less than the Ending Column Index.' newline];
        set(myhandles.status(1), 'String', myhandles.statusString, ...
        'ForegroundColor', [1 0 0]);
        jhEdit=findjobj(myhandles.status(1));
        jEdit=jhEdit.getComponent(0).getComponent(0);
        jEdit.setCaretPosition(jEdit.getDocument.getLength);
        guidata(varargin{1}, myhandles);
    else
        myhandles.statusString=[myhandles.statusString newline ...
        'Index of Starting Column for Group One set successfully. ' newline];
       set(myhandles.status(1), 'String', myhandles.statusString, ...
        'ForegroundColor', [0 0.8 0]);
        jhEdit=findjobj(myhandles.status(1));
        jEdit=jhEdit.getComponent(0).getComponent(0);
        jEdit.setCaretPosition(jEdit.getDocument.getLength);
        guidata(varargin{1}, myhandles);
    end
elseif source==myhandles.preProcessing(11)
    myhandles.groupOneEnd=str2double(get(myhandles.preProcessing(11), 'string'));
    if(isnan(str2double(get(myhandles.preProcessing(11), 'string'))))
        set(myhandles.preProcessing(11), 'string', 'End (Gr 1)');
        myhandles.groupOneEnd='';
        
        myhandles.statusString=[myhandles.statusString newline ...
        'Please enter a numeric input.' newline];
        set(myhandles.status(1), 'String', myhandles.statusString, ...
        'ForegroundColor', [1 0 0]);
        jhEdit=findjobj(myhandles.status(1));
        jEdit=jhEdit.getComponent(0).getComponent(0);
        jEdit.setCaretPosition(jEdit.getDocument.getLength);
        guidata(varargin{1}, myhandles);
    elseif myhandles.groupOneStart>myhandles.groupOneEnd
        set(myhandles.preProcessing(11), 'string', 'End (Gr 1)');
        myhandles.groupOneEnd='';
        myhandles.statusString=[myhandles.statusString newline ...
        'The Ending Column Index should be greater than the Starting Column Index.' newline];
        set(myhandles.status(1), 'String', myhandles.statusString, ...
        'ForegroundColor', [1 0 0]);
        jhEdit=findjobj(myhandles.status(1));
        jEdit=jhEdit.getComponent(0).getComponent(0);
        jEdit.setCaretPosition(jEdit.getDocument.getLength);
        guidata(varargin{1}, myhandles);
    else
        
        myhandles.statusString=[myhandles.statusString newline ...
        'Index of Ending Column for Group One set successfully. ' newline];
       set(myhandles.status(1), 'String', myhandles.statusString, ...
        'ForegroundColor', [0 0.8 0]);
        jhEdit=findjobj(myhandles.status(1));
        jEdit=jhEdit.getComponent(0).getComponent(0);
        jEdit.setCaretPosition(jEdit.getDocument.getLength);
        guidata(varargin{1}, myhandles);
    end
elseif source==myhandles.preProcessing(12)
    myhandles.groupTwoStart=str2double(get(myhandles.preProcessing(12), 'string'));   
    if(isnan(str2double(get(myhandles.preProcessing(12), 'string'))))
        set(myhandles.preProcessing(12), 'string', 'Start (Gr 2)');
        myhandles.groupTwoStart='';
        
        myhandles.statusString=[myhandles.statusString newline ...
        'Please enter a numeric input.' newline];
        set(myhandles.status(1), 'String', myhandles.statusString, ...
        'ForegroundColor', [1 0 0]);
        jhEdit=findjobj(myhandles.status(1));
        jEdit=jhEdit.getComponent(0).getComponent(0);
        jEdit.setCaretPosition(jEdit.getDocument.getLength);
        guidata(varargin{1}, myhandles);
    elseif myhandles.groupTwoStart>myhandles.groupTwoEnd
        set(myhandles.preProcessing(12), 'string', 'Start (Gr 2)');
        myhandles.groupTwoStart='';
        myhandles.statusString=[myhandles.statusString newline ...
        'The Starting Column Index should be less than the Ending Column Index.' newline];
        set(myhandles.status(1), 'String', myhandles.statusString, ...
        'ForegroundColor', [1 0 0]);
        jhEdit=findjobj(myhandles.status(1));
        jEdit=jhEdit.getComponent(0).getComponent(0);
        jEdit.setCaretPosition(jEdit.getDocument.getLength);
        guidata(varargin{1}, myhandles);
    else
        myhandles.statusString=[myhandles.statusString newline ...
        'Index of Starting Column for Group Two set successfully. ' newline];
       set(myhandles.status(1), 'String', myhandles.statusString, ...
        'ForegroundColor', [0 0.8 0]);
        jhEdit=findjobj(myhandles.status(1));
        jEdit=jhEdit.getComponent(0).getComponent(0);
        jEdit.setCaretPosition(jEdit.getDocument.getLength);
        guidata(varargin{1}, myhandles);
    end
elseif source==myhandles.preProcessing(13)
   myhandles.groupTwoEnd=str2double(get(myhandles.preProcessing(13), 'string'));
   if(isnan(str2double(get(myhandles.preProcessing(13), 'string'))))
        set(myhandles.preProcessing(13), 'string', 'End (Gr 2)');
        myhandles.groupTwoEnd='';
        
        myhandles.statusString=[myhandles.statusString newline ...
        'Please enter a numeric input.' newline];
        set(myhandles.status(1), 'String', myhandles.statusString, ...
        'ForegroundColor', [1 0 0]);
        jhEdit=findjobj(myhandles.status(1));
        jEdit=jhEdit.getComponent(0).getComponent(0);
        jEdit.setCaretPosition(jEdit.getDocument.getLength);
        guidata(varargin{1}, myhandles);
   elseif myhandles.groupTwoEnd<myhandles.groupTwoStart
        set(myhandles.preProcessing(13), 'string', 'End (Gr 2)');
        myhandles.groupTwoEnd='';
        myhandles.statusString=[myhandles.statusString newline ...
        'The Ending Column Index should be greater than the Starting Column Index.' newline];
        set(myhandles.status(1), 'String', myhandles.statusString, ...
        'ForegroundColor', [1 0 0]);
        jhEdit=findjobj(myhandles.status(1));
        jEdit=jhEdit.getComponent(0).getComponent(0);
        jEdit.setCaretPosition(jEdit.getDocument.getLength);
        guidata(varargin{1}, myhandles);
      
   else
       
        myhandles.statusString=[myhandles.statusString newline ...
        'Index of Ending Column for Group Two set successfully. ' newline];
        set(myhandles.status(1), 'String', myhandles.statusString, ...
        'ForegroundColor', [0 0.8 0]);
        jhEdit=findjobj(myhandles.status(1));
        jEdit=jhEdit.getComponent(0).getComponent(0);
        jEdit.setCaretPosition(jEdit.getDocument.getLength);
        guidata(varargin{1}, myhandles);
   end
elseif source==myhandles.preProcessing(8)
    myhandles.filterThreshold=str2double(get(myhandles.preProcessing(8), 'string'));
   if(isnan(str2double(get(myhandles.preProcessing(8), 'string'))))
        set(myhandles.preProcessing(8), 'string', '');
        myhandles.filterThreshold='';
        
        myhandles.statusString=[myhandles.statusString newline ...
        'Please enter a numeric input.' newline];
        set(myhandles.status(1), 'String', myhandles.statusString, ...
        'ForegroundColor', [1 0 0]);
        jhEdit=findjobj(myhandles.status(1));
        jEdit=jhEdit.getComponent(0).getComponent(0);
        jEdit.setCaretPosition(jEdit.getDocument.getLength);
        guidata(varargin{1}, myhandles);
   else
        myhandles.statusString=[myhandles.statusString newline ...
        'Filtering Threshold set sucessfully. ' newline ...
        'Genes with number of non-zero read counts less than the threshold shall be removed.' newline];
        set(myhandles.status(1), 'String', myhandles.statusString, ...
        'ForegroundColor', [0 0.8 0]);
        jhEdit=findjobj(myhandles.status(1));
        jEdit=jhEdit.getComponent(0).getComponent(0);
        jEdit.setCaretPosition(jEdit.getDocument.getLength);
        guidata(varargin{1}, myhandles);
   end
elseif source==myhandles.preProcessing(6)
    myhandles.numberGenes=str2double(get(myhandles.preProcessing(6), 'string'));
   if(isnan(str2double(get(myhandles.preProcessing(6), 'string'))))
        set(myhandles.preProcessing(6), 'string', '');
        myhandles.numberGenes='';
        
        myhandles.statusString=[myhandles.statusString newline ...
        'Please enter a numeric input.' newline];
        set(myhandles.status(1), 'String', myhandles.statusString, ...
        'ForegroundColor', [1 0 0]);
        jhEdit=findjobj(myhandles.status(1));
        jEdit=jhEdit.getComponent(0).getComponent(0);
        jEdit.setCaretPosition(jEdit.getDocument.getLength);
        guidata(varargin{1}, myhandles);
   else
        myhandles.statusString=[myhandles.statusString newline ...
        'Number of genes to be selected from the original file, set sucessfully. ' newline];
        set(myhandles.status(1), 'String', myhandles.statusString, ...
        'ForegroundColor', [0 0.8 0]);
        jhEdit=findjobj(myhandles.status(1));
        jEdit=jhEdit.getComponent(0).getComponent(0);
        jEdit.setCaretPosition(jEdit.getDocument.getLength);
        guidata(varargin{1}, myhandles);
   end   
end
% *************************************************************************
% Save the source in the myhandles structure
% *************************************************************************
if(~isempty(myhandles.groupOneStart) && ~isempty(myhandles.groupOneEnd) ...
      && ~isempty(myhandles.groupTwoStart) && ~isempty(myhandles.groupTwoEnd))
    set(myhandles.preProcessing(14), 'enable', 'on');
    myhandles.statusString=[myhandles.statusString newline ...
        '"Solve for DE" Push Button is enabled. Click for evaluating Differential ' ...
        'Expression.' newline];
        set(myhandles.status(1), 'String', myhandles.statusString, ...
        'ForegroundColor', [0 0.8 0]);
        jhEdit=findjobj(myhandles.status(1));
        jEdit=jhEdit.getComponent(0).getComponent(0);
        jEdit.setCaretPosition(jEdit.getDocument.getLength);
        guidata(varargin{1}, myhandles);
end


guidata(source, myhandles);