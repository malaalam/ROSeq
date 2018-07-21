function []=helpCall(source, eventdata, varargin)
% *************************************************************************
% This function "helpCall" defines the events that ensue when the user
% decides to open the Help Manual.
%
% *************************************************************************
% Define a structure that contains handles to all the objects in
% figureHandle
% *************************************************************************
myhandles=guidata(varargin{1});
if ismac
    open('ROSeq.pdf');
elseif ispc
    winopen('ROSeq.pdf');
end
% *************************************************************************
% Save the structure/data 'myhandles' in the 'figureHandle - varargin{1}'
% *************************************************************************
guidata(varargin{1}, myhandles);
end