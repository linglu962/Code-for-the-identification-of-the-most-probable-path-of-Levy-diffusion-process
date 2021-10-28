function A=basis_diff(z,n_basis_function)
A=ones(size(z,2),1);
for i=1:n_basis_function-1
    A=[A z'.^i]; %#ok<AGROW>
end