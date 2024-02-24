function pC = pCov(A,B)

[n1, k1]=size(A);
[n2, k2]=size(B);
n=n1+n2;
if(k1~=k2)
    disp('number of columns of A and B must be the same')
else
    xDiff=mean(A)-mean(B);       % mean difference row vector
    cA=covdiag(A);
    cB=covdiag(B);
    pC=n1/n*cA+n2/n*cB;          % pooled covariance matrix
end