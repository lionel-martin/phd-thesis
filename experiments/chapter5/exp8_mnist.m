%   Experiment of Table 5.3.
%   Computational benefits tested on the MNIST dataset.
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Created by Lionel Martin
MNIST = load_data('mnist', 'full', {'raw', 'nn'});
G = graph_from_nn(MNIST.nn, MNIST.dd, 'perplexity', 100);
G = gsp_estimate_lmax(G);

%% EBD method
fprintf('Running EBD light clustering...\t');
t_EBD = tic;

filt_param = struct('method', 1, 'order', 100);
[~, G.lk, ~] = estimate_lambda_k(G, k, param);
[~, cheb_coef] = jackson_cheby_poly_coefficients(0, G.lk, [0, G.lmax], param.order);
time_filter(ni, simu) = toc(t_FM);
lowpass_filter = @(x) gsp_cheby_op(G, cheb_coef, x);
assignment = ebd_light_clust(G, k, lowpass_filter);

time_EBD = toc(t_EBD);
ncut_EBD = compute_ncut(G, assignment);
mod_EBD = compute_modularity(G, assignment);
randind_EBD = compute_rand_index(MNIST.labels, assignment);

fprintf('Done!\nResults: time = %d, ncut = %.4f, mod = %.4f, rand index = %.4f.\n', time_EBD, ncut_EBD, mod_EBD, randind_EBD);

%% SC method
fprintf('Running Spectral Clustering (Shi & Malik)...\t');
parpool;

t_SC = tic;

[Uk, ek] = compute_eigen(G, k);
assignment = kmeans(Uk, k, 'Replicate', 150, 'Options', statset('UseParallel', 1));

time_SC = toc(t_SC);
ncut_SC = compute_ncut(G, assignment);
mod_SC = compute_modularity(G, assignment);
randind_SC = compute_rand_index(MNIST.labels, assignment);

fprintf('Done!\nResults: time = %d, ncut = %.4f, mod = %.4f, rand index = %.4f.\n', time_EBD, ncut_EBD, mod_EBD, randind_EBD);
