function []=saveMetaImages(myhandles)
% *************************************************************************
% This function "saveMetaImages" creates some meta images by randomly
% picking genes and showing their quality of fit.
%
%
% Last update: 27 July, 2017
% *************************************************************************
% Define a structure that contains handles to all the objects in
% figureHandle
% *************************************************************************
h=figure;
set(h, 'Visible', 'off');
for image_iter=1:5
    subplot(2, 5, image_iter)
    gene=ceil(rand(1)*myhandles.geneCount);
    a1=myhandles.results_groupOne{gene,1}.a;
    b1=myhandles.results_groupOne{gene,1}.b;
    A1=myhandles.results_groupOne{gene,1}.A;
    n1=myhandles.results_groupOne{gene,1}.number_of_bins;
    model1=myhandles.results_groupOne{gene,1}.model;
    actual1=myhandles.results_groupOne{gene,1}.actual;
    R21 = myhandles.results_groupOne{gene,1}.R2;
    rank1=myhandles.results_groupOne{gene,1}.rank;
    if(isnan(R21) || isinf(R21))
    else
    plot(rank1, model1, 'k');
    hold on
    plot(rank1, actual1, 'o');
    title([' SP1. Gene = ', num2str(gene), '. R2 = ', num2str(R21), '. a = ', num2str(a1), '. b = ', num2str(b1)])
    xlabel('Rank of Bin')
    end
   
    
    subplot(2, 5, image_iter+5)
    a2=myhandles.results_groupTwo{gene,1}.a;
    b2=myhandles.results_groupTwo{gene,1}.b;
    A2=myhandles.results_groupTwo{gene,1}.A;
    n2=myhandles.results_groupTwo{gene,1}.number_of_bins;
    model2=myhandles.results_groupTwo{gene,1}.model;
    actual2=myhandles.results_groupTwo{gene,1}.actual;
    R22 = myhandles.results_groupTwo{gene,1}.R2;
    rank2=myhandles.results_groupTwo{gene,1}.rank;
    if(isnan(R22) || isinf(R22))
    else
    plot(rank2, model2, 'k');
    hold on
    plot(rank2, actual2, 'o');
    title(['SP2. Gene = ', num2str(gene), '. R2 = ', num2str(R22), '. a = ', num2str(a2), '. b = ', num2str(b2)])
    xlabel('Rank of Bin')
    end

end
savefig(h, 'metaimage.fig')
close(h)