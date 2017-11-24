% Generic code for the evaluation of SC, CSC, dCSC and IASC with different
% parameters. This accepts static and dynamic clustering attempts.
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Modified by Lionel Martin

nb_simus = 50;
order = 150;

N = 15000;
k = 25;
s = 60;

ec = (s - sqrt(s))/(s + sqrt(s)*(k-1));
eps = ec / 2;
p = s*k/(N*(k-1)*eps+N-k);
q = p * eps;
gparams = struct('p', p, 'q', q);

% KEPT = [0, 0.25]; n_kept = numel(KEPT);
% D = linspace(10, 200, 5); nb_d = numel(D);
KEPT = [0, 0.25, 0.5]; n_kept = numel(KEPT);
D = logspace(log10(25), log10(200), 10); nb_d = numel(D);

REG_PARAM = [1e-4]; nb_reg = numel(REG_PARAM);

GAMMA = [2]; nb_gamma = numel(GAMMA);

% PARAMS = [1, 1; 1, 3; 3, 1]/100; %node removal, edge removal
PARAMS = [0.01, 0.01]; %node removal, edge removal
nb_params_set = size(PARAMS, 1);

measures = zeros(nb_d, nb_simus, nb_params_set, n_kept, nb_gamma, nb_reg, 9);


%%
for simu=1:nb_simus
    fprintf('Simu n°%d/%d\n', simu, nb_simus);
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
    fprintf('\t Done!\n');
    
    ABBA = 2*speye(G1.N) - G1.L;
    [U1_IASC, e1_IASC] = eigs(ABBA, min(20, G1.N), 'LM');

    for par_idx=1:nb_params_set
        fprintf('\tParam set %d/%d\n', par_idx, nb_params_set);
        fprintf('\t\tGraph Perturbation...');
        nrmi = PARAMS(par_idx, 1); ermi = PARAMS(par_idx, 2);
        param_modif = struct('p', p, 'q', q, 'prop_node_mod', nrmi, 'prop_edge_mod', ermi, 'until_spc_dec', true);
        G2 = modify_graph(G1, param_modif);
        G2.k = k;
        phi = G2.U;
        tic;
        [ass_SC2, ~, kmcost] = kmeans(phi, k, 'Replicates', 100, 'MaxIter', 150, 'Options', statset('UseParallel', 1));
        t_SC2 = toc;
        t_SC2 = t_SC2 + G2.time_spc_dec;
        C_SC2 = sum(kmcost);
        Ncut_SC2 = compute_ncut(G2, ass_SC2);
        save(sprintf('simulations/GT_simu_%d_param_%d.mat', simu, par_idx), 'G1', 'G2', 'C_SC2', 't_SC2', 'ass_SC2', 'param_modif');
        fprintf('\tDone!\n');
        
        tic;
        [e2_IASC, U2_IASC, ~] = eigenUpdate(G2, G1, e1_IASC, U1_IASC);
        [~, e_idx] = sort(diag(e2_IASC), 'descend');
        U2_IASC = U2_IASC(:, e_idx(1:k));
        ass_IASC = kmeans(U2_IASC, k, 'Replicates', 100, 'MaxIter', 150, 'Options', statset('UseParallel', 1));
        t_IASC = toc;
        C_IASC = compute_kmeans_cost(phi, ass_IASC);
        Ncut_IASC = compute_ncut(G2, ass_IASC);

        for di = 1:nb_d
            d_val = round(D(di));

            for gi = 1:nb_gamma
                for regi = 1:nb_reg
                    fprintf('\t\tCSC on G1 (d = %d)...\n', d_val);
                    params_CSC = struct('d_value', d_val, 'features', 1, 'assignment', 1, 'poly_order', order, ...
                        'sampling', 'VD', 'n_factor', GAMMA(gi), 'regu', REG_PARAM(regi));
                    [C_est, F1_est, time_CSC1, comp_terms_CSC1, params_CSC] = my_CSC(G1, params_CSC);
                    save(sprintf('simulations/CSC_G1_simu_%d_param_%d_d_%d.mat', simu, par_idx, d_val), 'G1', 'C_est', 'F1_est', 'time_CSC1', 'comp_terms_CSC1');

                    for pp = 1:numel(KEPT)
                        d = size(F1_est, 2); pi = KEPT(pp);
                        if pi == 0
                            if isfield(params_CSC, 'filt_sig'), params_CSC = rmfield(params_CSC, 'filt_sig'); end
                            if isfield(params_CSC, 'lk_est'), params_CSC = rmfield(params_CSC, 'lk_est'); end
                            if isfield(params_CSC, 'weight_VD'), params_CSC = rmfield(params_CSC, 'weight_VD'); end                            
                        else
                            params_CSC.filt_sig = F1_est(:, randperm(d, round(pi*d)));
                            params_CSC.weight_VD = comp_terms_CSC1.weight_VD;
                            params_CSC.lk_est = comp_terms_CSC1.lk_est;
                        end

                        if isfield(params_CSC, 'filt_sig'),
                           nb_sig = d-size(params_CSC.filt_sig, 2);
                        else
                            nb_sig = d;
                        end
                        fprintf('\t\tCSC on G2 (d_rem = %d)...\n', nb_sig);

                        % CSC on G2 with previous signals from CSC of G1
                        [ass2w1_CSC, F2w1_est, t_C2w1, comp_terms_CSC2, ~] = my_CSC(G2, params_CSC);
                        save(sprintf('simulations/CSC_G2_simu_%d_param_%d_d_%d_p_%f.mat', simu, par_idx, d_val, pi), 'ass2w1_CSC', 'F2w1_est', 't_C2w1', 'comp_terms_CSC2');

                        % CSC without compression on G1
    %                     fprintf('Running full kmeans...');
    %                     tic;
    %                     ass2w1_RF = kmeans(F2w1_est, k, 'Replicates', 100, 'MaxIter', 150, 'Options', statset('UseParallel',1));
    %                     t_RF2w1 = toc;
    %                     fprintf('\t\t\tDone.\n');

                        measures(di, simu, par_idx, pp, gi, regi, 1) = t_C2w1.total;
                        measures(di, simu, par_idx, pp, gi, regi,2) = compute_kmeans_cost(phi, ass2w1_CSC);
                        measures(di, simu, par_idx, pp, gi, regi, 3) = compute_ncut(G2, ass2w1_CSC);
    %                     measures(di, simu, par_idx, pp, 4) = t_RF2w1;
    %                     if isfield(t_C2w1, 'lk_est'), measures(di, simu, par_idx, pp, 4) = measures(di, simu, par_idx, pp, 4) + t_C2w1.lk_est; end
    %                     if isfield(t_C2w1, 'filtering'), measures(di, simu, par_idx, pp, 4) = measures(di, simu, par_idx, pp, 4) + t_C2w1.filtering; end
    %                     measures(di, simu, par_idx, pp, 5) = compute_kmeans_cost(phi, ass2w1_RF);
    %                     measures(di, simu, par_idx, pp, 6) = compute_ncut(G2, ass2w1_RF);
                        measures(di, simu, par_idx, pp, gi, regi, 7) = t_SC2;
                        measures(di, simu, par_idx, pp, gi, regi, 8) = C_SC2;
                        measures(di, simu, par_idx, pp, gi, regi, 9) = Ncut_SC2;
                    end
                end
            end
        end
    end
end

%% Additionnal variables
load('simulations/D_small.mat');
nb_d_sm = numel(D_small);
meas = measures;
meas_dims = size(measures);
meas_dims(1) = meas_dims(1) + nb_d_sm;
measures = zeros(meas_dims);
measures(nb_d_sm+1:end, :, :, :, :, :, :) = meas;
%% Additionnal computations
for simu=71:80
    for par_idx=1:nb_params_set
        fprintf('Simu n°%d/%d\n', simu, nb_simus);
        fprintf('\tGraphs loading...');
        load(sprintf('simulations/GT_simu_%d_param_%d.mat', simu, par_idx));
        fprintf('\tDone!\n');
        
        for di = 1:nb_d_sm
            d_val = round(D_small(di));
            
             for gi = 1:nb_gamma
                for regi = 1:nb_reg
                    fprintf('\t\tCSC on G1 (d = %d)...\n', d_val);
                    params_CSC = struct('d_value', d_val, 'features', 1, 'assignment', 1, 'poly_order', order, ...
                        'sampling', 'VD', 'n_factor', GAMMA(gi), 'regu', REG_PARAM(regi));
                    [C_est, F1_est, time_CSC1, comp_terms_CSC1, params_CSC] = my_CSC(G1, params_CSC);
                    save(sprintf('simulations/CSC_G1_simu_%d_param_%d_d_%d.mat', simu, par_idx, d_val), 'G1', 'C_est', 'F1_est', 'time_CSC1', 'comp_terms_CSC1');

                    for pp = 1:numel(KEPT)
                        d = size(F1_est, 2); pi = KEPT(pp);
                        if pi == 0
                            if isfield(params_CSC, 'filt_sig'), params_CSC = rmfield(params_CSC, 'filt_sig'); end
                            if isfield(params_CSC, 'lk_est'), params_CSC = rmfield(params_CSC, 'lk_est'); end
                            if isfield(params_CSC, 'weight_VD'), params_CSC = rmfield(params_CSC, 'weight_VD'); end                            
                        else
                            params_CSC.filt_sig = F1_est(:, randperm(d, round(pi*d)));
                            params_CSC.weight_VD = comp_terms_CSC1.weight_VD;
                            params_CSC.lk_est = comp_terms_CSC1.lk_est;
                        end

                        if isfield(params_CSC, 'filt_sig'),
                           nb_sig = d-size(params_CSC.filt_sig, 2);
                        else
                            nb_sig = d;
                        end
                        fprintf('\t\tCSC on G2 (d_rem = %d)...\n', nb_sig);

                        % CSC on G2 with previous signals from CSC of G1
                        [ass2w1_CSC, F2w1_est, t_C2w1, comp_terms_CSC2, ~] = my_CSC(G2, params_CSC);
                        save(sprintf('simulations/CSC_G2_simu_%d_param_%d_d_%d_p_%f.mat', simu, par_idx, d_val, pi), 'ass2w1_CSC', 'F2w1_est', 't_C2w1', 'comp_terms_CSC2');

                        % CSC without compression on G1
    %                     fprintf('Running full kmeans...');
    %                     tic;
    %                     ass2w1_RF = kmeans(F2w1_est, k, 'Replicates', 100, 'MaxIter', 150, 'Options', statset('UseParallel',1));
    %                     t_RF2w1 = toc;
    %                     fprintf('\t\t\tDone.\n');

                        measures(di, simu, par_idx, pp, gi, regi, 1) = t_C2w1.total;
                        measures(di, simu, par_idx, pp, gi, regi,2) = compute_kmeans_cost(phi, ass2w1_CSC);
                        measures(di, simu, par_idx, pp, gi, regi, 3) = compute_ncut(G2, ass2w1_CSC);
    %                     measures(di, simu, par_idx, pp, 4) = t_RF2w1;
    %                     if isfield(t_C2w1, 'lk_est'), measures(di, simu, par_idx, pp, 4) = measures(di, simu, par_idx, pp, 4) + t_C2w1.lk_est; end
    %                     if isfield(t_C2w1, 'filtering'), measures(di, simu, par_idx, pp, 4) = measures(di, simu, par_idx, pp, 4) + t_C2w1.filtering; end
    %                     measures(di, simu, par_idx, pp, 5) = compute_kmeans_cost(phi, ass2w1_RF);
    %                     measures(di, simu, par_idx, pp, 6) = compute_ncut(G2, ass2w1_RF);
                        measures(di, simu, par_idx, pp, gi, regi, 7) = t_SC2;
                        measures(di, simu, par_idx, pp, gi, regi, 8) = C_SC2;
                        measures(di, simu, par_idx, pp, gi, regi, 9) = Ncut_SC2;
                    end
                end
            end
        end
    end
end

%% Varying k - variables
K = round(logspace(log10(10), log10(100), 6)); nb_k = numel(K);
D_full = round([D_small, D]); nb_d_fl = numel(D_full);
measures_var_k = zeros(nb_k, nb_d_fl, nb_simus, nb_gamma, nb_reg, 3); % time, kmcost, ncut
%% Varying k - run
for ki = 1:nb_k
    k = round(K(ki));

    ec = (s - sqrt(s))/(s + sqrt(s)*(k-1));
    eps = ec / 2;
    p = s*k/(N*(k-1)*eps+N-k);
    q = p * eps;
    gparams = struct('p', p, 'q', q);

    for simu = sims{ki}
        fprintf('Simu n°%d/%d with k=%d\n', simu, nb_simus, k);

%         if k ~= 25
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
%         end

        for di = 1:nb_d_fl
            d_val = round(D_full(di));

            for gi = 1:nb_gamma
                for regi = 1:nb_reg

%                     if k == 25
%                         load(sprintf('simulations/CSC_G2_simu_%d_param_1_d_%d_p_%f.mat', simu, d_val, 0), 't_C2w1', 'ass2w1_CSC', 'comp_terms_CSC2');
%                         G1 = comp_terms_CSC2.G;
%                         time_CSC1 = t_C2w1;
%                         C_est = ass2w1_CSC;
%                         comp_terms_CSC1 = comp_terms_CSC2;
%                     else
                        fprintf('\t\tCSC on G1 (d = %d)...\n', d_val);
                        params_CSC = struct('d_value', d_val, 'features', 1, 'assignment', 1, 'poly_order', order, ...
                            'sampling', 'VD', 'n_factor', GAMMA(gi), 'regu', REG_PARAM(regi));
                        [C_est, F1_est, time_CSC1, comp_terms_CSC1, output_CSC] = my_CSC(G1, params_CSC);
%                     end

                    save(sprintf('simulations_varying_k_thesis/CSC_G1_simu_%d_k_%d_d_%d.mat', simu, k, d_val), 'C_est', 'F1_est', 'time_CSC1', 'comp_terms_CSC1', 'output_CSC');

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

                        save(sprintf('simulations_varying_k_thesis/GT_G1_simu_%d_k_%d.mat', simu, k), 'G1');
                    end
                end
            end
        end
    end
end

%% Varying k - postprocessing

SC_measures_var_k = zeros(nb_k, nb_simus, 3);
load('SC_measures_var_k.mat');

for ki = 1:nb_k
    k = round(K(ki));
    for ii = 1:nb_simus
        fprintf('Computing Graph #%d\n', ii);
        clear G1;
        load(sprintf('GT_G1_simu_%d_k_%d.mat', ii, k));
        if ~isfield(G1, 'cost_ncut_SC')
            G1.cost_ncut_SC = compute_ncut(G1, G1.ass_SC1);
            save(sprintf('GT_G1_simu_%d_k_%d.mat', ii, k), 'G1');
        end
        SC_measures_var_k(ki, ii, :) = [G1.time_SC_total, G1.cost_SC1, G1.cost_ncut_SC];
    end
end

save('SC_measures_var_k.mat', 'SC_measures_var_k');

%% Varying N - variables

NN = [1000, 10000]; %linspace(5000, 50000, 10); nb_N = numel(NN);
nb_N = numel(NN);
d_val = 50;

measures_var_n = zeros(nb_N, nb_simus, nb_params_set, n_kept, nb_gamma, nb_reg, 9);

%% Varying N - run

for simu=1:nb_simus
    fprintf('Simu n°%d/%d\n', simu, nb_simus);

    for ni = 1:nb_N
        N = NN(ni);
        p = s*k/(N*(k-1)*eps+N-k);
        q = p * eps;
        gparams = struct('p', p, 'q', q);

        fprintf('\tGraph generation (N=%d)...', N);

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
        fprintf('\t Done!\n');
        
%         ABBA = 2*speye(G1.N) - G1.L;
%         G1.ABBA = ABBA;
%         [U1_IASC, e1_IASC] = eigs(ABBA, min(20, G1.N), 'LM');
%         e1_IASC = diag(e1_IASC);

        for par_idx=1:nb_params_set
            fprintf('\tParam set %d/%d\n', par_idx, nb_params_set);
            fprintf('\t\tGraph Perturbation...');
            nrmi = PARAMS(par_idx, 1); ermi = PARAMS(par_idx, 2);
            param_modif = struct('p', p, 'q', q, 'prop_node_mod', nrmi, 'prop_edge_mod', ermi, 'until_spc_dec', true);
            G2 = modify_graph(G1, param_modif);
            G2.k = k;
            phi = G2.U;
            tic;
            [ass_SC2, ~, kmcost] = kmeans(phi, k, 'Replicates', 100, 'MaxIter', 150, 'Options', statset('UseParallel', 1));
            t_SC2 = toc;
            t_SC2 = t_SC2 + G2.time_spc_dec;
            C_SC2 = sum(kmcost);
            Ncut_SC2 = compute_ncut(G2, ass_SC2);
            save(sprintf('simulations_varying_n/GT_simu_%d_param_%d_N_%d.mat', simu, par_idx, N), 'G1', 'G2', 'C_SC2', 'Ncut_SC2', 't_SC2', 'ass_SC2', 'param_modif'); %, 'e1_IASC', 'U1_IASC');
            fprintf('\tDone!\n');

%             fprintf('\t\tIASC on G2...\n');
%             tic;
%             [e2_IASC, U2_IASC, ~] = eigenUpdate(G2, G1, e1_IASC, U1_IASC);
%             [~, e_idx] = sort(diag(e2_IASC), 'descend');
%             U2_IASC = U2_IASC(:, e_idx(1:k));
%             ass_IASC = kmeans(U2_IASC, k, 'Replicates', 100, 'MaxIter', 150, 'Options', statset('UseParallel', 1));
%             t_IASC = toc;
%             C_IASC = compute_kmeans_cost(phi, ass_IASC);
%             Ncut_IASC = compute_ncut(G2, ass_IASC);
%             save(sprintf('simulations_varying_n/IASC_simu_%d_param_%d_N_%d.mat', simu, par_idx, N), 'C_IASC', 'Ncut_IASC', 't_IASC', 'ass_IASC', 'param_modif');
%             fprintf('\tDone!\n');

            for gi = 1:nb_gamma
                for regi = 1:nb_reg
                    fprintf('\t\tCSC on G1 (d = %d)...\n', d_val);
                    params_CSC = struct('d_value', d_val, 'features', 1, 'assignment', 1, 'poly_order', order, ...
                        'sampling', 'VD', 'n_factor', GAMMA(gi), 'regu', REG_PARAM(regi));
                    [C_est, F1_est, time_CSC1, comp_terms_CSC1, params_CSC] = my_CSC(G1, params_CSC);
                    save(sprintf('simulations_varying_n/CSC_G1_simu_%d_param_%d_d_%d_N_%d.mat', simu, par_idx, d_val, N), 'C_est', 'F1_est', 'time_CSC1', 'comp_terms_CSC1');

                    for pp = 1:numel(KEPT)
                        d = size(F1_est, 2); pi = KEPT(pp);
                        if pi == 0
                            if isfield(params_CSC, 'filt_sig'), params_CSC = rmfield(params_CSC, 'filt_sig'); end
                            if isfield(params_CSC, 'lk_est'), params_CSC = rmfield(params_CSC, 'lk_est'); end
                            if isfield(params_CSC, 'weight_VD'), params_CSC = rmfield(params_CSC, 'weight_VD'); end                            
                        else
                            params_CSC.filt_sig = F1_est(:, randperm(d, round(pi*d)));
                            params_CSC.weight_VD = comp_terms_CSC1.weight_VD;
                            params_CSC.lk_est = comp_terms_CSC1.lk_est;
                        end

                        if isfield(params_CSC, 'filt_sig'),
                           nb_sig = d-size(params_CSC.filt_sig, 2);
                        else
                            nb_sig = d;
                        end
                        fprintf('\t\tCSC on G2 (d_rem = %d)...\n', nb_sig);

                        % CSC on G2 with previous signals from CSC of G1
                        [ass2w1_CSC, F2w1_est, t_C2w1, comp_terms_CSC2, ~] = my_CSC(G2, params_CSC);
                        C_CSC2 = compute_kmeans_cost(phi, ass2w1_CSC);
                        Ncut_CSC2 = compute_ncut(G2, ass2w1_CSC);
                        save(sprintf('simulations_varying_n/CSC_G2_simu_%d_param_%d_d_%d_N_%d_p_%f.mat', simu, par_idx, d_val, N, pi), 'C_CSC2', 'Ncut_CSC2', 'ass2w1_CSC', 'F2w1_est', 't_C2w1', 'comp_terms_CSC2');

                        % CSC without compression on G1
    %                     fprintf('Running full kmeans...');
    %                     tic;
    %                     ass2w1_RF = kmeans(F2w1_est, k, 'Replicates', 100, 'MaxIter', 150, 'Options', statset('UseParallel',1));
    %                     t_RF2w1 = toc;
    %                     fprintf('\t\t\tDone.\n');

                        measures_var_n(ni, simu, par_idx, pp, gi, regi, 1) = t_C2w1.total;
                        measures_var_n(ni, simu, par_idx, pp, gi, regi, 2) = C_CSC2;
                        measures_var_n(ni, simu, par_idx, pp, gi, regi, 3) = Ncut_CSC2;
    
%                         measures_var_n(ni, simu, par_idx, pp, gi, regi, 4) = t_IASC;
%                         measures_var_n(ni, simu, par_idx, pp, gi, regi, 5) = C_IASC;
%                         measures_var_n(ni, simu, par_idx, pp, gi, regi, 6) = Ncut_IASC;

                        measures_var_n(ni, simu, par_idx, pp, gi, regi, 7) = t_SC2;
                        measures_var_n(ni, simu, par_idx, pp, gi, regi, 8) = C_SC2;
                        measures_var_n(ni, simu, par_idx, pp, gi, regi, 9) = Ncut_SC2;
                    end
                end
            end
        end
    end
end
%% (TMP) Compute error difference
list_dir = dir;
nb_files = numel(list_dir);
SC_costs = zeros(nb_k, nb_simus);
for fi = 1:nb_files
    if strfind(list_dir(fi).name, 'GT_G1')
        fprintf('Loading %s\n', list_dir(fi).name);
        load(list_dir(fi).name);
        
        filename = strrep(list_dir(fi).name, '.mat', '');
        blocks = strsplit(filename, '_');
        SC_costs(kr==str2num(cell2mat(blocks(6))), str2num(cell2mat(blocks(4)))) = cost_SC1;
        
    end
end

%% (TMP) Recompute Kmeans for K=40
list_dir = dir;
nb_files = numel(list_dir);
SC_costs = zeros(nb_k, nb_simus);
for fi = 1:nb_files
    if strfind(list_dir(fi).name, 'GT_G1') & strfind(list_dir(fi).name, 'k_40')
        fprintf('Loading %s\n', list_dir(fi).name);
        load(list_dir(fi).name);

        t_SC1 = tic;
        [ass_SC1, ~, kmc] = kmeans(G1.U, 40, 'Replicates', 200, 'MaxIter', 150, 'Options', statset('UseParallel', 1));
        time_SC1 = toc(t_SC1) + G1.time_spc_dec;
        cost_SC1 = sum(kmcost);

        save(list_dir(fi).name, 'G1', 'ass_SC1', 'time_SC1', 'cost_SC1');
    end
end


%% Plotting dynamic
% close all;
min_perc = 40;
max_perc = 60;
for pari = 1:nb_params_set % 1: [1,1], 2: [0.5, 0.5], 3: [1, 0]
    time_SC = vec(measures(:, :, pari, :, :, :, 7));
    for meas = 2:3 % 2: kmeans cost, 3: ncut
        if meas == 2, meas_name = 'kmeans cost'; else meas_name = 'NCut'; end
        cost_SC = vec(measures(:, :, pari, :, :, :, 6+meas));
        %for csc_mode = 0:1 % 0: CSC, 1: RF
        csc_mode = 0;
            if csc_mode == 0, csc_mode_name = 'sampling'; else csc_mode_name = 'filtering'; end
%             data = permute(squeeze(measures(:, :, pari, :, csc_mode*3+meas)), [1, 3, 2]);
%             times = permute(squeeze(measures(:, :, pari, :, csc_mode*3+1)), [1, 3, 2]);
            data = permute(squeeze(measures(:, :, pari, :, :, :, csc_mode*3+meas)), [1, 3, 4, 2]);
            times = permute(squeeze(measures(:, :, pari, :, :, :, csc_mode*3+1)), [1, 3, 4, 2]);
            
            figure; hold on;
            low_tSC = prctile(time_SC, min_perc); med_tSC = median(time_SC); hig_tSC = prctile(time_SC, max_perc);
            low_cSC = prctile(cost_SC, min_perc); med_cSC = median(cost_SC); hig_cSC = prctile(cost_SC, max_perc);
            errorbarxy(med_cSC, med_tSC, med_cSC - low_cSC, hig_cSC - med_cSC, med_tSC - low_tSC, hig_tSC - med_tSC, {'gx', 'g', 'g'});
            scatter(measures(1, :, pari, 1, 1, 1, 6+meas), measures(1, :, pari, 1, 1, 1, 7), 15, 'g', 'filled');

            for ii = 1:size(data, 2)
                for jj = 1:size(data, 3)
                    sub_data = squeeze(data(:, ii, jj, :));
                    sub_times = squeeze(times(:, ii, jj, :));
                    low_meas = prctile(sub_data, min_perc, 2); med_meas = median(sub_data, 2); hig_meas = prctile(sub_data, max_perc, 2);
                    low_time = prctile(sub_times, min_perc, 2); med_time = median(sub_times, 2); hig_time = prctile(sub_times, max_perc, 2);
                    switch ii
                        case 1
                            color = 'b';
                        case 2
                            color = 'r';
                        case 3
                            color = 'k';
                        case 4
                            color = 'm';
                    end
                    switch jj
                        case 1
                            form = 'o';
                            line = '-';
                        case 2
                            form = 's';
                            line = '-.';
                        case 3
                            form = 'v';
                            line = ':';
                        case 4
                            form = '+--';
                    end
                    plot_cell = {strcat(color, form, line), color, color};
%                     errorbarxy(med_meas, med_time, med_meas - low_meas, hig_meas - med_meas, med_time - low_time, hig_time - med_time, plot_cell);
                    scatter(vec(sub_data), vec(sub_times), 15, color, 'filled', form);
                end
            end
            xlabel('Error');
            ylabel('Time [s]');
%             legend('exact', 'p = 0', 'p = 0.25');
%             legend('exact', '', 'n = 10%', '', 'n = 30%', '', 'n = 40%', '', 'n = 70%', '');
            hold off;

            title(strcat('Measure: ', meas_name, ', csc_mode: ', csc_mode_name, ', modifs (node, edge): ', num2str(PARAMS(pari, :))));
        %end
    end
end

%% Plotting varying k
% close all;

min_perc = 30;
max_perc = 70;

for meas = 2:3 % 2: kmeans cost, 3: ncut
    if meas == 2, meas_name = 'kmeans cost'; else meas_name = 'NCut'; end
        data = squeeze(measures_var_k(:, :, :, :, :, meas));
        times = squeeze(measures_var_k(:, :, :, :, :, 1));

        figure; hold on;
        for ii = 1:size(data, 1)
                sub_data = squeeze(data(ii, :, :));
                sub_times = squeeze(times(ii, :, :));
                low_meas = prctile(sub_data, min_perc, 2); med_meas = median(sub_data, 2); hig_meas = prctile(sub_data, max_perc, 2);
                low_time = prctile(sub_times, min_perc, 2); med_time = median(sub_times, 2); hig_time = prctile(sub_times, max_perc, 2);
                switch ii
                    case 1
                        color = 'b';
                    case 2
                        color = 'r';
                    case 3
                        color = 'k';
                    case 4
                        color = 'm';
                end
                plot_cell = {strcat(color, 'o-'), color, color};
                errorbarxy(med_meas, med_time, med_meas - low_meas, hig_meas - med_meas, 0, 0, plot_cell);
%                 scatter(vec(sub_data), vec(sub_times), 15, color, 'filled', 'o');
        end
        xlabel('Error');
        ylabel('Time [s]');
        legend(cellstr(num2str(round(K)', 'K = %-d')));
%             legend('exact', 'p = 0', 'p = 0.25');
%             legend('exact', '', 'n = 10%', '', 'n = 30%', '', 'n = 40%', '', 'n = 70%', '');
        hold off;

        title(strcat('Measure: ', meas_name));
    %end
end


figure; hold on;
data = squeeze(measures_diff(:, :, :, :, :, 2));
for ii = 1:size(data, 1)
    sub_data = squeeze(data(ii, :, :));
    low_meas = prctile(sub_data, min_perc, 2); med_meas = median(sub_data, 2); hig_meas = prctile(sub_data, max_perc, 2);
    switch ii
        case 1
            color = 'b';
        case 2
            color = 'r';
        case 3
            color = 'k';
        case 4
            color = 'm';
    end
    plot_cell = {strcat(color, 'o-'), color, color};
%     errorbarxy(D(1:10), med_meas, 0, 0, med_meas - low_meas, hig_meas - med_meas, plot_cell);
    scatter(repmat(D(1:10), 1, nb_simus), vec(sub_data), 15, color, 'filled', 'o');
end
xlabel('Number of signals');
ylabel('Error diff.');
legend(cellstr(num2str(round(K)', 'K = %-d')));
hold off;


%% Plotting varying n
% close all;

min_perc = 50;
max_perc = 50;

for meas = 1:3 % 1: time, 2: kmeans cost, 3: ncut
    figure; hold on;
    data = permute(squeeze(measures_var_n(:, 1:nb_simus, :, :, :, :, meas)), [3, 1, 2]);

    for ii = 1:size(data, 1)  % values of p
        sub_data = squeeze(data(ii, :, :));
        low_meas = prctile(sub_data, min_perc, 2); med_meas = median(sub_data, 2); hig_meas = prctile(sub_data, max_perc, 2);
        switch ii
            case 1
                color = 'b';
            case 2
                color = 'r';
            case 3
                color = 'k';
        end
        plot_cell = {strcat(color, 'o-'), color, color};
        errorbar(NN, med_meas, med_meas - low_meas, hig_meas - med_meas, plot_cell{1});

    end

    SC_data = squeeze(measures_var_n(:, 1:nb_simus, :, 1, :, :, 6+meas));
    low_meas = prctile(SC_data, min_perc, 2); med_meas = median(SC_data, 2); hig_meas = prctile(SC_data, max_perc, 2);
    errorbar(NN, med_meas, med_meas - low_meas, hig_meas - med_meas, 'go-');

    hold off;
    leg = cellstr(num2str(KEPT', 'p = %-.2f'));
    leg{4} = 'SC';
    legend(leg);

    xlabel('Number of nodes');
    switch meas
        case 1
            ylabel('Time [s]');
        case 2
            ylabel('KM Cost');
        case 3
            ylabel('NCut');
    end
end
