function [ prodtrace ] = computeTrace( inverseSigma, dInverseSigma )

MAXSIZE = 20000000;
n = size(dInverseSigma,2);
p = symrcm(inverseSigma);
prodtrace = 0;
blockSize = floor(MAXSIZE/n);

for i=1:blockSize:n % loop over columns of dInverseSigma
%     i
    toColumn = min([i+blockSize-1 n]);
    x = inverseSigma(p,p)\dInverseSigma(p,p(i:toColumn));
    prodtrace = prodtrace + trace(x(i:toColumn,:)); 
end
