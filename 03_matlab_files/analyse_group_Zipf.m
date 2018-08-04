function [ results ] = analyse_group_Zipf(subPopulation, stats)
% *************************************************************************
% This function 'analyze_group_Zipf()' makes a call to the function
% "find_a_b" which finds the optimal values of the parameters a and b for a
% given gene and subpopulation read count data.
%
% *************************************************************************
% Import File
% *************************************************************************
real_data_set=subPopulation;
% *************************************************************************
% Normalize Data Set
% *************************************************************************
normalized_data_set=real_data_set;
% *************************************************************************
% Apply Negative Log Likelihood
% *************************************************************************
results=find_a_b(normalized_data_set, stats);


end

