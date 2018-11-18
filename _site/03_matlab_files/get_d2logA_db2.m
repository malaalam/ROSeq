function [ d2logA_db2 ] = get_d2logA_db2( u2, v2, du2_db, dv2_db)

% *************************************************************************
% This function 'get_d2logA_db2' finds the double derivative of A with
% respect to b. This double derivative is evaluated at the optimal (a_hat,
% b_hat).
% *************************************************************************

num1=v2.*du2_db;
num2=u2.*dv2_db;
den1=v2.^2;
d2logA_db2=(num1-num2)./den1;

end