function [ sparse_mat ] = gsp_sparse_mat_sym(nb_rows, density, param)
%   GSP_SPARSE_MAT_SYM computes exactly a symmetrical matrix of the correct density
%
%   sparse_mat = gsp_sparse_mat_sym(N, density) returns a sparse
%   matrix of size N x N of given density. The total number of non-zero
%   element is N * N * density. The matrix is always symmetrical.
%
%   sparse_mat = gsp_sparse_mat_sym(N, density, param) accepts additional
%   parameters:
%   * 'empty_diag' (default true) forces the diagonal elements to be zero.
%   The density is ensured with both behaviors.
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Created by Lionel Martin
if nargin < 3, param = struct; end
if ~isfield(param, 'empty_diag'), param.empty_diag = true; end

if param.empty_diag, diag = 1; else diag = 0; end

triu_idx = find(triu(ones(nb_rows), diag) == 1);

nb_elems_tot = numel(triu_idx);
nb_elems = round(nb_elems_tot * density);

indices = triu_idx(randperm(nb_elems_tot, nb_elems));

[Is, Js] = ind2sub([nb_rows, nb_rows], indices);
sparse_mat = sparse([Is, Js], [Js, Is], 1, nb_rows, nb_rows);
sparse_mat(1:nb_rows+1:end) = sparse_mat(1:nb_rows+1:end) / 2;

end