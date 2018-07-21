function [ du2_da] = get_du2_da( coefficients, r)

% *************************************************************************
% This function 'get_du2_da' finds the first derivative of u2 with
% respect to a. This first derivative is evaluated at the optimal (a_hat,
% b_hat).
% *************************************************************************

    a=coefficients(1);
    b=coefficients(2);
    N=length(r);
    
    num1=(N+1-r).^b;
    num2=log(r);
    num3=log(N+1-r);
    den1=r.^a;
    du2_da=...
        sum(num1.*num2.*num3./den1);
end