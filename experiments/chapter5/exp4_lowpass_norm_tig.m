%   Experiment of Figure 5.3.
%   Norm Tig using lowpass itersine filter highlights the central nodes (low norm).
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Created by Lionel Martin
gsp_start;

N = 300;
G = gsp_sensor(N);
G = gsp_estimate_lmax(G);

signal = speye(N);
filter = gsp_design_itersine(G, 5);

param = struct('method', 'cheby', 'order', 50);
filt_sig = gsp_filter_analysis(G, filter{1}, signal, param);
tig_norm = sum(filt_sig.^2, 1);

gsp_plot_signal(G, tig_norm);
caxis([0, 0.5*max(tig_norm)]);
