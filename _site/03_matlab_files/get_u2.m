function [ u2] = get_u2( coefficients, r)
% *************************************************************************
% This function 'get_u2' finds u2.
% *************************************************************************

    a=coefficients(1);
    b=coefficients(2);
    N=length(r);
    
    num1=(N+1-r).^b;
    num2=log(N+1-r);
    den1=r.^a;
    u2=...
        -sum(num1.*num2./den1);
end

