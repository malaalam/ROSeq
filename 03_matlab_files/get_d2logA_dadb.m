function [ d2logA_dadb ] = get_d2logA_dadb( u2, v2, du2_da, dv2_da)

% *************************************************************************
% This function 'get_d2logA_dadb' finds the double derivative of A with
% respect to a and b. This double derivative is evaluated at the optimal (a_hat,
% b_hat).
% *************************************************************************

num1=v2.*du2_da;
num2=u2.*dv2_da;
den1=v2.^2;
d2logA_dadb=(num1-num2)./den1;

end