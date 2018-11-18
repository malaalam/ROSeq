function [ dv2_db] = get_dv2_db( coefficients, r)

% *************************************************************************
% This function 'get_dv2_db' finds the first derivative of v2 with
% respect to b. This double derivative is evaluated at the optimal (a_hat,
% b_hat).
% *************************************************************************

    a=coefficients(1);
    b=coefficients(2);
    N=length(r);
    
    num1=(N+1-r).^b;
    num2=log(N+1-r);
    den1=r.^a;
    dv2_db=...
        sum(num1.*num2./den1);
end