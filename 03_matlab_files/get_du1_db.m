function [ du1_db] = get_du1_db( coefficients, r)

% *************************************************************************
% This function 'get_du1_db' finds the first derivative of u1 with
% respect to b. This first derivative is evaluated at the optimal (a_hat,
% b_hat).
% *************************************************************************

    a=coefficients(1);
    b=coefficients(2);
    N=length(r);
    
    num1=(N+1-r).^b;
    num2=log(r);
    num3=log(N+1-r);
    den1=r.^a;
    du1_db=...
        sum(num1.*num2.*num3./den1);
end