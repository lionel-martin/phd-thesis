function [ graph_filter ] = approx_filter_allenzhu( G, lk, m, param )
%   APPROX_FILTER_ALLENZHU generates the internal function for plotting
%   the approximation of the lowpass filter with Allen-Zhu's method.
%
%   graph_filter = approx_filter_allenzhu(G, lk, m) returns the
%   anonymous function generating the graph filter based on the threshold
%   lk, the maximum eigenvalue G.lmax and the order of the approximation m
%   (default 30).
%
%   graph_filter = approx_filter_allenzhu(G, lk, m, param) supports additional
%   parameters:
%   * 'verbose' (default 1) defines the quantity of outputed information.
%   * 'kap' (default 0.01) specifies how quickly the step function must
%   decay.
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Created by Lionel Martin

if nargin < 4, param = struct; end
if ~isfield(param,'verbose'), param.verbose = 1; end;

if nargin < 3, m = 30; end
param.order = m;

construct_G = @(y) struct('lmax', G.lmax, 'L', spdiags(y(:), 0, speye(numel(y))), 'N', numel(y));

graph_filter = @(x) gsp_filter_new_ideal_lowpass(construct_G(x), lk, ones(size(x(:))), param);

end
