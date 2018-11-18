function []=saveMetaData(myhandles)
% *************************************************************************
% This function "saveMetaData" saves some intermediate results as a csv
% file by default. This is, in case, one needs to perform some further
% post-processing.
%
%
% Last update: 27 July, 2017

%[file,path] = uiputfile('*.csv','Results');
a1=zeros(length(myhandles.results_groupOne),1);
b1=zeros(length(myhandles.results_groupOne),1);
A1=zeros(length(myhandles.results_groupOne),1);
n1=zeros(length(myhandles.results_groupOne),1);
R21=zeros(length(myhandles.results_groupOne),1);
a2=zeros(length(myhandles.results_groupOne),1);
b2=zeros(length(myhandles.results_groupOne),1);
A2=zeros(length(myhandles.results_groupOne),1);
n2=zeros(length(myhandles.results_groupOne),1);
R22=zeros(length(myhandles.results_groupOne),1);

for gene=1:length(myhandles.results_groupOne)
    a1(gene,1)=myhandles.results_groupOne{gene,1}.a;
    b1(gene,1)=myhandles.results_groupOne{gene,1}.b;
    A1(gene,1)=myhandles.results_groupOne{gene,1}.A;
    n1(gene,1)=myhandles.results_groupOne{gene,1}.number_of_bins;
    R21(gene,1) = myhandles.results_groupOne{gene,1}.R2;
    
    a2(gene,1)=myhandles.results_groupTwo{gene,1}.a;
    b2(gene,1)=myhandles.results_groupTwo{gene,1}.b;
    A2(gene,1)=myhandles.results_groupTwo{gene,1}.A;
    n2(gene,1)=myhandles.results_groupTwo{gene,1}.number_of_bins;
    R22(gene,1) = myhandles.results_groupTwo{gene,1}.R2;
end




T = table(a1,b1,A1, n1, R21, a2, b2, A2, n2, R22);

writetable(T, 'metadata.csv')

