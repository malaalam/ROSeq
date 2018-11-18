function [ v2] = get_v2( coefficients, r)

% *************************************************************************
% This function 'get_v2' finds v2.
% *************************************************************************

    a=coefficients(1);
    b=coefficients(2);
    N=length(r);
    num1=(N+1-r).^b;
    den1=r.^a;
    v2=...
        sum(num1./den1);
end

