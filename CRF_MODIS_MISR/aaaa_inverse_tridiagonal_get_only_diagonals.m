function inverse = aaaa_inverse_tridiagonal_get_only_diagonals(A)

% returns sparse inverse matrix, with entries only on main and two adjacent
% diagonals

% returns the inverse of a tridiagonal matrix, according to paper:
% R.-S. Ran, T.-Z Huang - An inversion algorithm for a banded matrix
%
% but of special kind of tridiagonal matrix, where two lower diagonals are
% all -1's

len = length(A);
b = diag(A);
inverse = sparse(len, len);

% calculate theta and phi from page 6 (pp. 1704), the first two for theta
% are -1th and 0th, last two for phi are n+1-st and n+2-nd
theta = zeros(1, len + 2);
phi = zeros(1, len + 2);

theta(1) = 0;
theta(2) = 1;
phi(len + 1) = 1;
phi(len + 2) = 0;
for i = 1 : len
    theta(i + 2) = b(i) * theta(i + 1) - theta(i);
end;
for i = len : -1 : 1
    phi(i) = b(i) * phi(i + 1) - phi(i + 2);
end;

for i = 1 : len
    for j = max([1, i - 1]) : min([len, i + 1])
        inverse(i, j) = theta(i + 1) * phi(j + 1) / theta(len + 2);
        inverse(j, i) = inverse(i, j);
    end;
end;