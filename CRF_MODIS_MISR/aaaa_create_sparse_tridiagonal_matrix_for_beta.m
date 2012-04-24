function tridiag = aaaa_create_sparse_tridiagonal_matrix_for_beta(n)

e = ones(n, 1);
tridiag = spdiags([-e 2*e -e], -1:1, n, n);
tridiag(1, 1) = 1;
tridiag(n, n) = 1;