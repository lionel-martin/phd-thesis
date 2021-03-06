%EXPERIMENT_CLUSTERING_SBM Script SBM clustering experiments (Figure 4.4)

% Copyright (C) 2016 Johan Paratte, Lionel Martin.
% This file is part of the Reproducible Research code base for the Fast
% Eigenspace Approximation using Random Signals (FEARS) method.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

% If you use this code please kindly cite
%     Paratte, Johan and Martin, Lionel. 
%     "Fast Eigenspace Approximation using Random Signals."
%     arXiv preprint arXiv:1611.00938 (2016).
% https://arxiv.org/abs/1611.00938

%% Data init
N = 5000;
k = 20;
deg = 16;

nb_simus = 50;
eps_resolution = 20;

only_our_method = 0;
compute_pm = 0;

params.verbose = 0;
params.lap_type = 'normalized';

if ~only_our_method
    measures_rnd = zeros(nb_simus, 20, 8);
    measures_nmi = zeros(nb_simus, 20, 8);
    times = zeros(nb_simus, 20, 8);
end
%%
for simu=1:nb_simus
    for eps=1:eps_resolution
        
        fprintf('*****\neps : %d/%d, simu: %d/%d\n*****\n', eps, eps_resolution, simu, nb_simus);
        
        paramsG = struct('deg', deg, 'eps_factor', eps/eps_resolution);        
        G = generate_sbm(N, k, paramsG);
        G = gsp_create_laplacian(G, 'normalized');
        G = gsp_estimate_lmax(G);
        G.k = k;
        
        params.order = 50;
        params.poly_order = 50;
        params.R = randn(N, k)/sqrt(k);
        
        % (Shi & Malik) Spectral Clustering
        if ~only_our_method
            disp('SC');
            t_SC = tic;
            Uk = compute_eigen(G, k);
            assign_SC = kmeans(Uk, k, 'Replicates', 10, 'Options', setstats('UseParallel', 1));
            times(simu, eps, 1) = toc(t_SC);
            measures_rnd(simu, eps, 1) = rand_index(G.info.node_com, assign_SC, 'adjusted');
            measures_nmi(simu, eps, 1) = compute_nmi(G.info.node_com, assign_SC);
        end
        
        % Our method
        disp('FEARS 50');
        G = gsp_create_laplacian(G, 'normalized');
        G = gsp_estimate_lmax(G);
        t_US = tic;
        basis = gsp_eigenspace_estimation(G, k, params);
        assign_us = kmeans(basis, k, 'Replicates', 10, 'Options', setstats('UseParallel', 1));
        times(simu, eps, 2) = toc(t_US);
        measures_rnd(simu, eps, 2) = rand_index(G.info.node_com, assign_us, 'adjusted');
        measures_nmi(simu, eps, 2) = compute_nmi(G.info.node_com, assign_us);

        % CSC method
        if ~only_our_method
            disp('CSC 50');
            t_CSC = tic;
            params_CSC = struct('features', 0, 'assignment', 1, 'poly_order', ...
              params.poly_order, 'sampling', 'VD', 'n_factor', 0.15*N/(k*log(k)));
            assign_CSC = my_CSC(G, params_CSC);
            times(simu, eps, 3) = toc(t_CSC);
            measures_rnd(simu, eps, 3) = rand_index(G.info.node_com, assign_CSC, 'adjusted');
            measures_nmi(simu, eps, 3) = compute_nmi(G.info.node_com, assign_CSC);
        end
        
        % Our method
        disp('FEARS 250');
        params.order = 250;
        params.poly_order = 250;
        t_US = tic;
        [basis, ~, out_params] = gsp_eigenspace_estimation(G, k, params);
        t_filt250 = toc(t_US);
        assign_us = kmeans(basis, k, 'Replicates', 10, 'Options', setstats('UseParallel', 1));
        times(simu, eps, 4) = toc(t_US);
        measures_rnd(simu, eps, 4) = rand_index(G.info.node_com, assign_us, 'adjusted');
        measures_nmi(simu, eps, 4) = compute_nmi(G.info.node_com, assign_us);

        t_CUS = tic;
        params_CSC = struct('features', 0, 'assignment', 1, 'poly_order', ...
              params.poly_order, 'sampling', 'uniform', 'n_factor', 0.15*N/(k*log(k)), ...
              'filt_sig', basis, 'lk_est', out_params.lk);
        assign_CFEARS = my_CSC(G, params_CSC);

        times(simu, eps, 5) = toc(t_CUS) + t_filt250;
        measures_rnd(simu, eps, 5) = rand_index(G.info.node_com, assign_CFEARS, 'adjusted');
        measures_nmi(simu, eps, 5) = compute_nmi(G.info.node_com, assign_CFEARS);

        % CSC method
        if ~only_our_method
            disp('CSC 250');
            t_CSC = tic;
            params_CSC = struct('features', 0, 'assignment', 1, 'poly_order', ...
              params.poly_order, 'sampling', 'VD', 'n_factor', 0.15*N/(k*log(k)));
            assign_CSC = my_CSC(G, params_CSC);
            times(simu, eps, 6) = toc(t_CSC);
            measures_rnd(simu, eps, 6) = rand_index(G.info.node_com, assign_CSC, 'adjusted');
            measures_nmi(simu, eps, 6) = compute_nmi(G.info.node_com, assign_CSC);
        end

        % PM method
        if ~only_our_method && compute_pm
            disp('PM p=5');
            t_PM = tic;
            feat = speye(G.N) - G.L;
            feat_pw = (feat^11)*params.R;
            [basis, ~, ~] = svd(feat_pw, 'econ');

            assign_pm = kmeans(basis, k, 'Replicates', 10);
            times(simu, eps, 7) = toc(t_PM);
            measures_rnd(simu, eps, 7) = rand_index(G.info.node_com, assign_pm, 'adjusted');
            measures_nmi(simu, eps, 7) = compute_nmi(G.info.node_com, assign_pm);

            disp('PM p=25');
            t_PM = tic;
            feat = speye(G.N) - G.L;
            feat_pw = (feat^51)*params.R;
            [basis, ~, ~] = svd(feat_pw, 'econ');

            assign_pm = kmeans(basis, k, 'Replicates', 10);
            times(simu, eps, 8) = toc(t_PM);
            measures_rnd(simu, eps, 8) = rand_index(G.info.node_com, assign_pm, 'adjusted');
            measures_nmi(simu, eps, 8) = compute_nmi(G.info.node_com, assign_pm);
        end
    end
end
