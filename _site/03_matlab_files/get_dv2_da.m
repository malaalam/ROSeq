function [ dv2_da] = get_dv2_da( coefficients, r)

% *************************************************************************
% This function 'get_dv2_da' finds the first derivative of v2 with
% respect to a. This first derivative is evaluated at the optimal (a_hat,
% b_hat).
% *************************************************************************

    a=coefficients(1);
    b=coefficients(2);
    N=length(r);
    
    num1=(N+1-r).^b;
    num2=log(r);
    den1=r.^a;
    dv2_da=...
        -sum(num1.*num2./den1);
end