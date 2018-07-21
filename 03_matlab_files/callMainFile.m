function []=callMainFile(source, eventdata, varargin)
% *************************************************************************
% This function "callMainFile" define the events that ensue when the user
% presses the "Solve for DE" push button.
try
% *************************************************************************
% Define a structure that contains handles to all the objects in
% figureHandle
% *************************************************************************
myhandles=guidata(varargin{1});

% *************************************************************************
% Input Data for Data One and Two - Make sure the data is normalized
% *************************************************************************

if isempty(myhandles.rawFileLocation)
    filterThreshold=0;
    subPopulation=myhandles.normalizedData.data;
    geneNames=myhandles.normalizedData.textdata(2:end,1);
elseif isempty(myhandles.normalizedFileLocation)
    filterThreshold=myhandles.filterThreshold;
    subPopulation=myhandles.rawData.data;
    geneNames=myhandles.rawData.textdata(2:end,1);
end
   
% *************************************************************************
% FILTER
% *************************************************************************

rowIndex=sum(subPopulation > 0, 2) >= filterThreshold;
subPopulation= subPopulation(rowIndex, :);
geneCount=size(subPopulation, 1);
%set(myhandles.preProcessing(6), 'string', num2str(min(geneCount, myhandles.numberGenes)));
myhandles.statusString=[myhandles.statusString newline ...
        'Total number of genes ' ...
        'considered for Differential Expression equals: ' ...
        num2str(min(geneCount, myhandles.numberGenes)) newline];
       set(myhandles.status(1), 'String', myhandles.statusString, ...
        'ForegroundColor', [0 0.8 0]);
        jhEdit=findjobj(myhandles.status(1));
        jEdit=jhEdit.getComponent(0).getComponent(0);
        jEdit.setCaretPosition(jEdit.getDocument.getLength);
        guidata(varargin{1}, myhandles);
guidata(varargin{1}, myhandles);
geneNames=geneNames(rowIndex, :);

% *************************************************************************
% MEDIAN NORMALIZATION
% *************************************************************************
if isempty(myhandles.normalizedFileLocation)
    libSize=sum(subPopulation,1);
    subPopulation=(subPopulation./libSize*median(libSize));
end

geneCount=min(myhandles.numberGenes, geneCount);
startRow1=1;
startColumn1=myhandles.groupOneStart;
endColumn1=myhandles.groupOneEnd;
numberRow1=geneCount;
endRow1=startRow1+numberRow1-1;
startRow2=1;
startColumn2=myhandles.groupTwoStart;
endColumn2=myhandles.groupTwoEnd;
numberRow2=geneCount;
endRow2=startRow2+numberRow2-1;
% *************************************************************************
% Save the structure/data 'myhandles' in the 'figureHandle - source'
% *************************************************************************
myhandles.geneCount=geneCount;
guidata(source, myhandles)
subPopulation_One=subPopulation(startRow1:endRow1, startColumn1:endColumn1);
subPopulation_Two=subPopulation(startRow2:endRow2, startColumn2:endColumn2);
geneNames=geneNames(startRow1:endRow1);

% *************************************************************************
% Find statistics for the available data. These shall be used later while
% setting the bins
% *************************************************************************

for gene=1:geneCount
   maxds(gene,1)=max(max(subPopulation_One(gene,:)), ...
       max(subPopulation_Two(gene,:))); 
   minds(gene,1)=min(min(subPopulation_One(gene,:)), ...
       min(subPopulation_Two(gene,:)));
   meands(gene,1)=mean([subPopulation_One(gene,:) ...
       subPopulation_Two(gene,:)]);
   stdds(gene,1)=std([subPopulation_One(gene, :) ...
       subPopulation_Two(gene, :)]);
   ceilds(gene,1)=ceil((maxds(gene,1)-meands(gene,1))/stdds(gene,1));
   floords(gene,1)=floor((minds(gene,1)-meands(gene,1))/stdds(gene,1));
   log2FC(gene,1)=abs(log2(mean(subPopulation_One(gene, :))/...
   mean(subPopulation_Two(gene, :))));
end

save('stats', 'meands','stdds', 'ceilds', 'floords', 'minds', 'maxds')

% *************************************************************************
% Find the Best Fitingt Parameters individually for Group One and Two
% *************************************************************************

results_groupOne=analyse_group_Zipf(subPopulation_One);
myhandles.results_groupOne=results_groupOne;
results_groupTwo=analyse_group_Zipf(subPopulation_Two);
myhandles.results_groupTwo=results_groupTwo;


% *************************************************************************
% Determine Differential Expression between Groups One and Two
% *************************************************************************
results_DE=get_DE(results_groupOne, results_groupTwo);
myhandles.results_DE=results_DE;
% *************************************************************************
% Store the T-statistic as a separate column
% *************************************************************************

Tset=[];
for i=1:length(results_DE)
    Tset=[Tset; results_DE{i}.T];
end

% *************************************************************************
% Convert the T-Statistic into a pValue. Find p Adjusted Value
% *************************************************************************

pVal = chi2cdf(Tset, 2, 'upper');
[padj]=mafdr(pVal,'BHFDR', true);



% *************************************************************************
% Save pValue, pAdjusted Value and log2FC in a separate structure
% *************************************************************************

results.T=Tset;
results.pVal=pVal;
results.padj=padj;
results.geneNames=geneNames;
results.log2FC=log2FC;
myhandles.results=results;

% *************************************************************************
% Save Meta Data Results - Optional (Helps in Troubleshooting!)
% *************************************************************************
saveMetaData(myhandles);
saveMetaImages(myhandles);
saveSingleMetaImage(myhandles);
saveFitImage(myhandles);

% *************************************************************************
% Add a message to indicate that DE analysis is complete
% *************************************************************************
set(myhandles.preProcessing(15), 'enable', 'on');
 myhandles.statusString=[myhandles.statusString newline ...
        'DE analysis completed successfully. ' newline ...
        'Click to save results as a *.csv file format.' newline];
       set(myhandles.status(1), 'String', myhandles.statusString, ...
        'ForegroundColor', [0 0.8 0]);
        jhEdit=findjobj(myhandles.status(1));
        jEdit=jhEdit.getComponent(0).getComponent(0);
        jEdit.setCaretPosition(jEdit.getDocument.getLength);
        guidata(varargin{1}, myhandles);
        
% *************************************************************************
% Save the structure/data 'myhandles' in the 'figureHandle - source'
% *************************************************************************
guidata(source, myhandles)
catch e
    myhandles.statusString = [myhandles.statusString newline ...
    'The identifier was: ' e.identifier newline...
    'There was an error! The message was: ' e.message newline ...
    'Check "callMainFile.m" for further details.' newline];
       set(myhandles.status(1), 'String', myhandles.statusString, ...
        'ForegroundColor', [240/255,128/255,128/255]);
        jhEdit=findjobj(myhandles.status(1));
        jEdit=jhEdit.getComponent(0).getComponent(0);
        jEdit.setCaretPosition(jEdit.getDocument.getLength);
        guidata(varargin{1}, myhandles);
end