function[NLL] = minimize_NLL( coefficients, r, readCount)
% *************************************************************************
% This function 'minimize_NLL' finds the negative log likelihood. This is
% minimized in other functions by iterating over values of a and b.
% *************************************************************************

    a=coefficients(1);
    b=coefficients(2);
    N=length(r);
    sumReadCount=sum(readCount);
    A=1./sum(((N+1-r).^b)./(r.^a));
    NLL=...
        a*sum(readCount.*log(r))...
        -b*sum(readCount.*log(N+1-r))...
        -sumReadCount*log(A);
end