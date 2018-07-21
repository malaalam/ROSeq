function []=saveSingleMetaImage(myhandles)
% *************************************************************************
% This function "saveSingleMetaImage" creates just one meta image for
% export to the Bioinformatics Paper.
%
%
% Last update: 27 July, 2017
% *************************************************************************
% Define a structure that contains handles to all the objects in
% figureHandle
% *************************************************************************
h=figure;
set(h, 'Visible', 'off');
    set(h,'defaulttextinterpreter','latex');
    gene=340;
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
    plot(rank1, model1, '-g',...
    'LineWidth',4);
    hold on
    plot(rank1, actual1, '-r', 'LineWidth', 4);
    %text_Display = ['Gene Name :', myhandles.results.geneNames(gene, 1) ,...
    %   '$$R^{2}$$ = ', R21, '$$\widehat{a}$$ = ', a1, ...
    %    '$$\widehat{b}$$ = ', b1];
    text_Display = ['Gene Name : ', 'SRSF4 ', ...
        '$$R^{2}$$ = ', round(R21,3), '$$\widehat{a}$$ = ', round(a1,3), ...
        '$$\widehat{b}$$ = ', round(b1,3)];
    text(80,0.02,text_Display,'Interpreter','latex', 'FontSize', 20)
    xlabel('Rank of Bin', 'FontSize', 20, 'FontWeight', 'bold')
    ylabel('Normalized Frequency', 'FontSize', 20, 'FontWeight', 'bold')
    legend({'Model', 'Actual'}, 'FontSize', 20, 'FontWeight', 'bold')
    end
    savefig(h, 'singlemetaimage.fig')
close(h)