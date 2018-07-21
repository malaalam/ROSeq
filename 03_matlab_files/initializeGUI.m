function[]=initializeGUI(source, eventdata, varargin)

% *************************************************************************
% This function 'initializeGUI()' creates a new instance of the Graphical
% User Interface - ROSeq: Rank-Ordered Distribution for scRNA-Seq Read 
% Counts.
%
% *************************************************************************
% Take away any previous fields or previous instances of ESA
% *************************************************************************
clc
clear all
close all
% *************************************************************************
% Set the units for position and sizing to be pixels subsequently
% *************************************************************************
set(0,'units','pixels');
% *************************************************************************
% Set the Width and Height of ESA to be 80% of the Screen Size
% *************************************************************************
pixelsScreenSize = get(0,'screensize');
edgeLength=min(pixelsScreenSize(3), pixelsScreenSize(4));
figureHandle = figure('units','pixels', ...
    'position', ...
    [0.1*edgeLength 0.1*edgeLength ...
    0.8*edgeLength 0.8*edgeLength], ...
    'menubar','none', ...
    'name','ROSeq', ...
    'numbertitle','off', ...
    'resize','on');

% *************************************************************************
% Define a structure that contains handles to all the objects in
% figureHandle
% *************************************************************************
myhandles = guihandles(figureHandle);
% *************************************************************************
% Define the various menus and submenus. Provide an option to open these
% menus using the keyboard as well.
% *************************************************************************
myhandles.File = ...
    uimenu(figureHandle, 'Label', 'File', ...
    'Accelerator', 'F', 'Enable', 'on');
myhandles.File_NewProject= uimenu(myhandles.File, 'Label', 'New Project', ...
    'Accelerator', 'N', 'Enable','on');
myhandles.Open_Project= uimenu(myhandles.File, 'Label', 'Open Project', ...
    'Accelerator', 'O', 'Enable','on');
myhandles.Save_NewProject= uimenu(myhandles.File, 'Label', 'Save Project', ...
    'Accelerator', 'S', 'Enable', 'on');
myhandles.SaveAs_NewProject= uimenu(myhandles.File, 'Label', 'Save As', ...
    'Accelerator', 'R', 'Enable', 'on');
myhandles.Help = uimenu(figureHandle, 'Label', 'Help', ...
    'Accelerator', 'H', 'Enable', 'on');
myhandles.HelpContents = uimenu(myhandles.Help, 'Label', 'Help Contents', ...
    'Accelerator', 'C', 'Enable', 'on');
myhandles.About = uimenu(myhandles.Help, 'Label', 'About', ...
    'Accelerator', 'E', 'Enable', 'on');
% *************************************************************************

% Save the structure/data 'myhandles' in the 'figureHandle'
% *************************************************************************
guidata(figureHandle, myhandles);
% *************************************************************************
% Define some initial fields for use later in ESA
% *************************************************************************
defineProperties(figureHandle);
% *************************************************************************
% Define the various uicontrol panels, buttons et cetera
% *************************************************************************
definePanels(figureHandle);
% *************************************************************************
% Define the callbacks i.e. how the various panels interact with each other
% *************************************************************************
defineCallbacks(figureHandle);
end





