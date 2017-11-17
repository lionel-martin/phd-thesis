%   Experiment of Figure 5.4.
%   Norm Tig on unbalanced SBM results with all the lowest norm in the same cluster.
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Created by Lionel Martin
gsp_start;

N = 700;
k = 5;
z = [ones(1, 300), 2*ones(1, 200), 3*ones(1, 100), 4*ones(1, 50), 5*ones(1, 50)];
gparams = struct('p', 0.16, 'q', 0.04, 'z', z);
G = gsp_stochastic_block_graph(N, k, gparams);

signal = speye(N);
filter = gsp_design_heat(G);

param = struct('method', 'cheby', 'order', 50);
filt_sig = gsp_filter_analysis(G, filter, signal, param);
tig_norm = sum(filt_sig.^2, 1);

gsp_plot_signal(G, tig_norm);
caxis([0.95*min(tig_norm), 1.2*min(tig_norm)]);
