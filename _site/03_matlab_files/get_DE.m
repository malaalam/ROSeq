function [ results_DE ] = get_DE( results_1, results_2)

% *************************************************************************
% This function 'get_DE' finds the differential expression between results
% arising from the two sub-populations
% *************************************************************************

I_1=get_I(results_1);
I_2=get_I(results_2);

% *************************************************************************
% The results from differential expression are saved in a new structure
% *************************************************************************
for gene=1:length(I_1)
    
    I1=[I_1{gene}.I11 I_1{gene}.I12; I_1{gene}.I21 I_1{gene}.I22];
    I2=[I_2{gene}.I11 I_2{gene}.I12; I_2{gene}.I21 I_2{gene}.I22];
    results_DE{gene,1}.I1=I1;
    results_DE{gene,1}.I2=I2;
    results_DE{gene, 1}.V1=inv(I1);
    results_DE{gene, 1}.V2=inv(I2);
    results_DE{gene, 1}.m=results_1{gene}.number_of_bins;
    results_DE{gene, 1}.n=results_2{gene}.number_of_bins;
    results_DE{gene, 1}.a1=results_1{gene}.a;
    results_DE{gene, 1}.b1=results_1{gene}.b;
    results_DE{gene, 1}.a2=results_2{gene}.a;
    results_DE{gene, 1}.b2=results_2{gene}.b;
    results_DE{gene, 1}.w=results_DE{gene, 1}.n/ (results_DE{gene, 1}.m + results_DE{gene, 1}.n);
    results_DE{gene, 1}.T=...
        (results_DE{gene}.m * results_DE{gene}.n/(results_DE{gene}.m + results_DE{gene}.n))* ...
        [results_DE{gene}.a1 - results_DE{gene}.a2; results_DE{gene}.b1 - results_DE{gene}.b2].' * ...
        inv(results_DE{gene}.w * results_DE{gene}.V1 + (1-results_DE{gene}.w)* results_DE{gene}.V2) * ...
        [results_DE{gene}.a1 - results_DE{gene}.a2; results_DE{gene}.b1 - results_DE{gene}.b2];
    
end

