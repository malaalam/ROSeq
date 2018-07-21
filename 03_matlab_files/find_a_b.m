function [ results ] =find_a_b( data_set )
% *************************************************************************
% This function 'find_a_b' finds the optimal values of the parameters a and
% b through Maximising the Log Likelihood.
%
% *************************************************************************
% Define a structure that contains handles to all the objects in
% figureHandle
% *************************************************************************

number_of_genes=size(data_set, 1);

% *************************************************************************
% For EACH gene, evaluate the best fitting coefficients
% *************************************************************************
fig_wait=waitbar(0,'Modeling a Sub-Population, one gene at a time...');
for gene=1:number_of_genes
    waitbar(gene/number_of_genes);
    ds=data_set(gene, :);
    load('stats');
    step=0.05;
    rs=[];
% *************************************************************************
% Determine the frequency or number of occurences of read counts within a
% certain bin. Bin Width is set to be equal to 0.05 times sigma. 
% *************************************************************************
    
    for i=floords(gene,1):step:ceilds(gene,1)-step
        LL=meands(gene,1)+i*stdds(gene,1);
        UL=meands(gene,1)+(i+step)*stdds(gene,1);
        c1=find((ds<UL));
        c2=find((ds>=LL));
        rs=[rs length(intersect(c1,c2))];
    end
    
    fds=rs((~isnan(rs)));
    number_of_bins=length(fds);
    rank=1:1:number_of_bins;
    read_count_sorted=sort(fds, 'descend');
    normalized_read_count_sorted=read_count_sorted/sum(read_count_sorted); %actual
    fun = ...
        @(coefficients)minimize_NLL(coefficients, rank, normalized_read_count_sorted);
    coefficients0=[0.25, 3];
    options = optimset('MaxFunEvals',1000);


% *************************************************************************
% Use fminsearch() to determine the best fitting a and b through an
% iterative procedure
% *************************************************************************

    [correctCoefficients,errorValue] = fminsearch(fun,coefficients0,options);
    if isnan(errorValue)
        a=NaN;
        b=NaN;
    else
        a=correctCoefficients(1);
        b=correctCoefficients(2);
        
    end
    
    A=1./sum((number_of_bins+1-rank).^b./rank.^a);
    f=A*((number_of_bins+1-rank).^b)./(rank.^a); % model
    SS_res=sum((normalized_read_count_sorted-f).^2);
    SS_tot=sum((normalized_read_count_sorted-mean(normalized_read_count_sorted)).^2);
    R2=1-SS_res/SS_tot;
% *************************************************************************
% Save the results as a structure variable to be returned from this function
% l************************************************************************
    
    results{gene, 1}.a=a;
    results{gene, 1}.b=b;
    results{gene, 1}.A=A;
    results{gene, 1}.rank=rank;
    results{gene, 1}.number_of_bins=number_of_bins;
    results{gene, 1}.sum_read_counts=sum(normalized_read_count_sorted);
    results{gene, 1}.R2 = R2;
    results{gene, 1}.model=f;
    results{gene, 1}.actual=normalized_read_count_sorted;
end
close(fig_wait)



end

