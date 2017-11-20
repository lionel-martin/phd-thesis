function [ sparse_mat ] = gsp_sparse_mat(nb_rows, nb_cols, density)
%   GSP_SPARSE_MAT computes exactly a sparse matrix of the correct density
%
%   sparse_mat = gsp_sparse_mat(N, M, density) returns a sparse
%   matrix of size N x M of given density. The total number of non-zero
%   element is N * M * density.
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Created by Lionel Martin

nb_elems_tot = nb_rows * nb_cols;
nb_elems = round(nb_elems_tot * density);

indices = randperm(nb_elems_tot, nb_elems);

[Is, Js] = ind2sub([nb_rows, nb_cols], indices);
sparse_mat = sparse(Is, Js, 1, nb_rows, nb_cols);

end