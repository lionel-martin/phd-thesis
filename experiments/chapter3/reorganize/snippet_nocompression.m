parpool;
N = 15000;
nb_k = 6; K = round(logspace(1, 2, nb_k));

load('/mnt/data/thesis_data/simulations_varying_k_thesis/nocomp_meas_var_k.mat');
load('/mnt/data/thesis_data/simulations/D_small.mat');
D = logspace(log10(25), log10(200), 10);
D = round([D_small, D]); nb_d = numel(D);

for ki = 1:4
    for di=1:nb_d
        for simu=1:50
            fprintf('Treating simu=%d, k=%d, d=%d\n', simu, ki, di);
            load(sprintf('/mnt/data/thesis_data/simulations_varying_k_thesis/CSC_G1_simu_%d_k_%d_d_%d.mat', simu,round( K(ki)), round(D(di))));
            G = comp_terms_CSC1.G;
            
            fprintf('Full kmeans...')
            tic;
            C_est = kmeans(F1_est , K(ki), 'Replicates', 100, 'MaxIter', 150, 'Options', statset('UseParallel',1));
            time_CSC1.k_means_full=toc;
            fprintf('\t\tDone.\n')

            time_CSC1 = rmfield(time_CSC1, 'total');
            time_CSC1.total = sum(cell2mat(struct2cell(time_CSC1)));
            comp_terms_CSC1.IDX_LD = C_est;
            comp_terms_CSC1.ind_obs = 1:N;

            kmcost = compute_kmeans_cost(G.U, C_est);
            ncut = compute_ncut(G, C_est);
            nocomp_meas_var_k(ki, di, simu, :) = [time_CSC1.total, kmcost, ncut];

            save(sprintf('/mnt/data/thesis_data/simulations_varying_k_thesis/nocomp_CSC_sim_%d_k_%d_d_%d.mat', simu, round(K(ki)), round(D(di))), 'comp_terms_CSC1', 'time_CSC1', 'C_est');
        end
    end
end
save('/mnt/data/thesis_data/simulations_varying_k_thesis/nocomp_meas_var_k.mat', 'nocomp_meas_var_k');