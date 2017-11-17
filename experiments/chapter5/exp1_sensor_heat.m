%   Experiment of Figure 5.1.
%   Diffusion of three Kronecker randomly selected in a Sensor graph.
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Created by Lionel Martin
gsp_start;

N = 400;
k = 3;
G = gsp_sensor(N);
G = gsp_estimate_lmax(G);

indices = randperm(N, k);
signal = sparse(indices, 1, 1, N, 1);

heat = gsp_design_heat(G);

param = struct('method', 'cheby', 'order', 50);
filt_sig = gsp_filter_analysis(G, heat, signal, param);

plotparams.vertex_highlight = indices;
gsp_plot_signal(G, filt_sig, plotparams);

