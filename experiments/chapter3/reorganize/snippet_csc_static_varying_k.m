clear all;

order = 150;
nb_simus = 50;
sims = {1:0, 1:0, 1:0, 13:50, 1:50, 1:50};

N = 15000;
k = 25;
s = 60;

KEPT = [0, 0.25, 0.5]; n_kept = numel(KEPT);
REG_PARAM = [1e-4]; nb_reg = numel(REG_PARAM);
GAMMA = [2]; nb_gamma = numel(GAMMA);
PARAMS = [0.01, 0.01]; %node removal, edge removal
nb_params_set = size(PARAMS, 1);

D = logspace(log10(25), log10(200), 10); nb_d = numel(D);
load('/mnt/data/thesis_data/simulations/D_small.mat');
nb_d_sm = numel(D_small);
D_full = round([D_small, D]); nb_d_fl = numel(D_full);

K = round(logspace(log10(10), log10(100), 6)); nb_k = numel(K);

measures_var_k = zeros(nb_k, nb_d_fl, nb_simus, nb_gamma, nb_reg, 3); % time, kmcost, ncut

for ki = 1:nb_k
    k = round(K(ki));

    ec = (s - sqrt(s))/(s + sqrt(s)*(k-1));
    eps = ec / 2;
    p = s*k/(N*(k-1)*eps+N-k);
    q = p * eps;
    gparams = struct('p', p, 'q', q);

    for simu = sims{ki}
        fprintf('Simu n°%d/%d with k=%d\n', simu, nb_simus, k);

            fprintf('\tGraph generation...');
            infeasible = 1;
            while infeasible
                G1 = gsp_stochastic_block_graph(N, k, gparams);
                G1 = gsp_create_laplacian(G1, 'normalized');
                G1.k = k;

                try
                    tSC = tic;
                    [U, e] = eigs(G1.L, G1.k, 'sm');
                    [G1.e, feat_idx] = sort(diag(e));
                    G1.U = U(:, feat_idx);
                    infeasible = 0;
                catch ME
                    infeasible = 1;
                    warning('Catched %s', ME.identifier);
                end
                G1.time_spc_dec = toc(tSC);
            end

            fprintf('\tDone!\n');

        for di = 1:nb_d_fl
            d_val = round(D_full(di));

            for gi = 1:nb_gamma
                for regi = 1:nb_reg
                        fprintf('\t\tCSC on G1 (d = %d)...\n', d_val);
                        params_CSC = struct('d_value', d_val, 'features', 1, 'assignment', 1, 'poly_order', order, ...
                            'sampling', 'VD', 'n_factor', GAMMA(gi), 'regu', REG_PARAM(regi));
                        [C_est, F1_est, time_CSC1, comp_terms_CSC1, output_CSC] = my_CSC(G1, params_CSC);


                    save(sprintf('/mnt/data/thesis_data/simulations_varying_k_thesis/CSC_G1_simu_%d_k_%d_d_%d.mat', simu, k, d_val), 'C_est', 'F1_est', 'time_CSC1', 'comp_terms_CSC1', 'output_CSC');

                    measures_var_k(ki, di, simu, gi, regi, 1) = time_CSC1.total;
                    measures_var_k(ki, di, simu, gi, regi, 2) = compute_kmeans_cost(G1.U, C_est);
                    measures_var_k(ki, di, simu, gi, regi, 3) = compute_ncut(G1, C_est);
                    
                    if di == 1
                        t_SC1 = tic;
                        [G1.ass_SC1, ~, kmcost] = kmeans(G1.U, k, 'Replicates', 100, 'MaxIter', 150, 'Options', statset('UseParallel', 1));
                        G1.time_SC_kmeans = toc(t_SC1);
                        G1.time_SC_total = G1.time_SC_kmeans + G1.time_spc_dec;
                        G1.cost_SC1 = sum(kmcost);
                        G1.cost_ncut_SC = compute_ncut(G1, G1.ass_SC1);

                        save(sprintf('/mnt/data/thesis_data/simulations_varying_k_thesis/GT_G1_simu_%d_k_%d.mat', simu, k), 'G1');
                    end
                end
            end
        end
    end
end
