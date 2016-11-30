function [kQ,cQ,bQ,fQ]=evals(xq,k,c,b,f)
kQ=k(xq);
cQ=c(xq);
bQ=b(xq);
fQ=f(xq);
%if abs(xq-0.2) <= 0.025
%  fQ = fQ - 100
%end