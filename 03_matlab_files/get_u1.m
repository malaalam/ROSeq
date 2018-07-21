function [ u1] = get_u1( coefficients, r)

% *************************************************************************
% This function 'get_u1' finds u1. 
% *************************************************************************

    a=coefficients(1);
    b=coefficients(2);
    N=length(r);
    
    num1=(N+1-r).^b;
    num2=log(r);
    den1=r.^a;
    u1=...
        sum(num1.*num2./den1);
end

