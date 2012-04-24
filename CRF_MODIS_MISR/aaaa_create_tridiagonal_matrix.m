function tridiag = aaaa_create_tridiagonal_matrix(main_diag, lower_diag, upper_diag)

len = length(main_diag);
if (length(lower_diag) ~= (len - 1))
    lower_diag = lower_diag(1) * ones(1, len - 1);
end;
if (length(upper_diag) ~= (len - 1))
    upper_diag = upper_diag(1) * ones(1, len - 1);
end;

tridiag = diag(main_diag) + diag(upper_diag, 1) + diag(lower_diag, -1);