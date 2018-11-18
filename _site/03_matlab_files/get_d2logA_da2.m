function [ d2logA_da2 ] = get_d2logA_da2( u1, v1, du1_da, dv1_da)
% *************************************************************************
% This function 'get_d2logA_da2' finds the double derivative of A with
% respect to a. This double derivative is evaluated at the optimal (a_hat,
% b_hat).
%
% *************************************************************************
% Define a structure that contains handles to all the objects in
% figureHandle
% *************************************************************************

num1=v1.*du1_da;
num2=u1.*dv1_da;
den1=v1.^2;
d2logA_da2=(num1-num2)./den1;

end

