parpool;
n = round(80 * log(40));
N = 15000;
nb_k = 6; K = round(logspace(1, 2, nb_k));

load('/mnt/data/thesis_data/simulations_varying_k_thesis/const_comp_meas_var_k.mat');
load('/mnt/data/thesis_data/simulations/D_small.mat');
D = logspace(log10(25), log10(200), 10);
D = round([D_small, D]); nb_d = numel(D);

for ki = [2, 3, 5, 6]
    for di=1:nb_d
        for simu=1:50
            fprintf('Treating simu=%d, k=%d, d=%d\n', simu, ki, di);
            load(sprintf('/mnt/data/thesis_data/simulations_varying_k_thesis/CSC_G1_simu_%d_k_%d_d_%d.mat', simu,round( K(ki)), round(D(di))));
            G = comp_terms_CSC1.G;
            ind_obs = datasample(1:N, n, 'Replace', false, 'Weights', comp_terms_CSC1.weight_VD);
            X_lk_est_DS = F1_est(ind_obs, :);

            fprintf('Low-dimensional kmeans...')
            tic; 
            IDX_LD = kmeans(X_lk_est_DS , K(ki), 'Replicates', 100, 'MaxIter', 150, 'Options', statset('UseParallel',1));
            time_CSC1.k_means_low_dim=toc;
            fprintf('\t\tDone.\n')

            fprintf('Interpolation of cluster indicators...\n\n')
            tic;
            C_obs_LD = sparse(1:n, IDX_LD, 1, n, K(ki));
            [~,JCH_HP] = jackson_cheby_poly_coefficients(G.lk, G.lmax, [0, G.lmax], output_CSC.poly_order);

            C_est = zeros(N, size(C_obs_LD,2));
            parfor k=1:size(C_obs_LD,2)
               c_obs = C_obs_LD(:,k);
               C_est(:, k) = interpolate_on_complete_graph(c_obs, ind_obs, @(x)gsp_cheby_op(G, JCH_HP, x), output_CSC.regu, G.N, output_CSC.solver);
            end

            [~, assignment] = max(C_est./repmat(sqrt(sum(C_est.^2, 1)), N, 1), [], 2);
            time_CSC1.interpolation=toc;
            fprintf('\t\tDone.\n')
            time_CSC1 = rmfield(time_CSC1, 'total');
            time_CSC1.total = sum(cell2mat(struct2cell(time_CSC1)));
            comp_terms_CSC1.IDX_LD = IDX_LD;
            comp_terms_CSC1.ind_obs = ind_obs;
            
            C_est = assignment;
            kmcost = compute_kmeans_cost(G.U, C_est);
            ncut = compute_ncut(G, C_est);
            const_comp_meas_var_k(ki, di, simu, :) = [time_CSC1.total, kmcost, ncut];

            save(sprintf('/mnt/data/thesis_data/simulations_varying_k_thesis/const_comp_CSC_sim_%d_k_%d_d_%d.mat', simu, round(K(ki)), round(D(di))), 'comp_terms_CSC1', 'time_CSC1', 'C_est');
        end
    end
end
save('/mnt/data/thesis_data/simulations_varying_k_thesis/const_comp_meas_var_k.mat', 'const_comp_meas_var_k');