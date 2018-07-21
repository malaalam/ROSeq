function [ d2logA_dbda ] = get_d2logA_dbda( u1, v1, du1_db, dv1_db)

% *************************************************************************
% This function 'get_d2logA_da2' finds the double derivative of A with
% respect to a and b. This double derivative is evaluated at the optimal (a_hat,
% b_hat).
% *************************************************************************

num1=v1.*du1_db;
num2=u1.*dv1_db;
den1=v1.^2;
d2logA_dbda=(num1-num2)./den1;

end