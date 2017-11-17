function [Uk, ek, timing] = compute_eigen( G, k, params )
%   COMPUTE_EIGEN solves the eigendecomposition for the k smallest
%   eigenpairs. The eigensolver of choice is eigs.
%
%   [Uk, ek] = compute_eigen(G, k) returns the eigenvectors Uk
%   associated with the k smallest eigenvalues ek. k must be strictly
%   positive and less than G.N.
%
%   [Uk, ek, timing] = compute_eigen(..., params) supports additional
%   parameters:
%   * 'timing' (default false) outputs the time required by the method.
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Created by Lionel Martin
assert(k<=G.N && k > 0);

time_eigen = 0;

if nargin < 3
    params = struct('timing', false);
elseif ~isfield(params, 'timing')
    params.timing = false;
end

if params.timing, t_eigen = tic; end
    
[Uk, ek] = eigs(G.L, k, 'sm');
[ek, eid] = sort(diag(ek));
Uk = Uk(:, eid);

if params.timing, time_eigen = toc(t_eigen); end
if nargout == 3, timing = time_eigen; end
