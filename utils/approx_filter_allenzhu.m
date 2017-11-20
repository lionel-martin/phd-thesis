function [ graph_filter ] = approx_filter_allenzhu( G, lk, m, param )


if nargin < 4, param = struct; end
if ~isfield(param,'verbose'), param.verbose = 1; end;

if nargin < 3, m = 30; end
param.order = m;

graph_filter = @(x) gsp_filter_new_ideal_lowpass(G, lk, x', param);

end
