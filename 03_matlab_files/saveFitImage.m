function []=saveFitImage(myhandles)
% *************************************************************************
% This function "saveFitImage" creates and saves an image that presents the
% R2fit as a histogram.
%
% Last update: 27 July, 2017
% *************************************************************************
% Define a structure that contains handles to all the objects in
% figureHandle
% *************************************************************************
h=figure; % Group One
set(h, 'Visible', 'off');

R2_set_groupOne=zeros(myhandles.geneCount,1);
R2_set_groupTwo=zeros(myhandles.geneCount,1);
for gene=1:myhandles.geneCount
   R2_set_groupOne(gene,1)=myhandles.results_groupOne{gene,1}.R2;
   R2_set_groupTwo(gene,1)=myhandles.results_groupTwo{gene,1}.R2; 
end

edges=0:0.05:1;
histogram(R2_set_groupOne, edges, 'facecolor', 'magenta');
xlabel('R^{2} - Coefficient of Determination', 'FontSize', 20,'FontWeight', 'bold' )
ylabel('Count of Genes', 'FontSize', 20, 'FontWeight', 'bold')
ylim([0 9000])
savefig(h, 'fit_groupOne_image.fig')
close(h)

i=figure; % Group Two
set(i, 'Visible', 'off');
histogram(R2_set_groupTwo, edges, 'facecolor', 'magenta');
xlabel('R^{2} - Coefficient of Determination', 'FontSize', 20, 'FontWeight', 'bold')
ylabel('Count of Genes', 'FontSize', 20, 'FontWeight', 'bold')
ylim([0 9000])
savefig(i, 'fit_groupTwo_image.fig')

close(i)