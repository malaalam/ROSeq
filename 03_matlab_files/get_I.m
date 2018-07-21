function [ I ] = get_I(results)

% *************************************************************************
% This function 'get_I' finds the Fisher Information matrix.
% *************************************************************************

number_of_genes=size(results, 1);

for gene=1:number_of_genes
   
    %number_of_bins=results{gene}.number_of_bins;
    %rank=[1:1:number_of_bins];
    rank=results{gene}.rank;
    coefficients=[results{gene}.a, results{gene}.b];
    
    
    u1=get_u1(coefficients, rank);
    v1=get_v1(coefficients, rank);
    u2=get_u2(coefficients, rank);
    v2=get_v2(coefficients, rank);
    
    
    
    du1_da=get_du1_da(coefficients, rank);
    du1_db=get_du1_db(coefficients, rank);
    du2_da=get_du2_da(coefficients, rank);
    du2_db=get_du2_db(coefficients, rank);
    dv1_da=get_dv1_da(coefficients, rank);
    dv1_db=get_dv1_db(coefficients, rank);
    dv2_da=get_dv2_da(coefficients, rank);
    dv2_db=get_dv2_db(coefficients, rank);
    
    d2logA_da2= get_d2logA_da2( u1, v1, du1_da, dv1_da);
    d2logA_db2= get_d2logA_db2( u2, v2, du2_db, dv2_db);
    d2logA_dbda= get_d2logA_dbda( u1, v1, du1_db, dv1_db);
    d2logA_dadb = get_d2logA_dadb( u2, v2, du2_da, dv2_da);

    I{gene, 1}.I11=-d2logA_da2;
    I{gene, 1}.I12=-d2logA_dadb;
    I{gene, 1}.I21=-d2logA_dbda;
    I{gene, 1}.I22=-d2logA_db2;    
end

end

