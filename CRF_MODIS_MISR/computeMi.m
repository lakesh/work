function [ mi ] = computeMi( inverseSigma, b )
% first solution
% mi = inverseSigma\b;

% second solution
p = symrcm(inverseSigma);
mi = inverseSigma(p,p)\b(p);
mi(p) = mi;
end