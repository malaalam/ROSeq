function [ du1_da] = get_du1_da( coefficients, r)
% *************************************************************************
% This function 'get_du1_da' finds the first derivative of u1 with
% respect to a. This first derivative is evaluated at the optimal (a_hat,
% b_hat).
% ************************************************************************

    a=coefficients(1);
    b=coefficients(2);
    N=length(r);
    
    num1=(N+1-r).^b;
    num2=(log(r)).^2;
    den1=r.^a;
    du1_da=...
        -sum(num1.*num2./den1);
end
