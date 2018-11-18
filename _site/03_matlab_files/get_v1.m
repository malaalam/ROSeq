function [ v1] = get_v1( coefficients, r)

% *************************************************************************
% This function 'get_v1' finds v1.
% *************************************************************************

    a=coefficients(1);
    b=coefficients(2);
    N=length(r);
    num1=(N+1-r).^b;
    den1=r.^a;
    v1=...
        sum(num1./den1);
end

