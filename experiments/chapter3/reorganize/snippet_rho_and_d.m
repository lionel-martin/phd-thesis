n = 10000;
k = 25;
s = 60;
ec = (s - sqrt(s))/(s + sqrt(s)*(k-1));
eps = ec / 2;
p = s*k/(n*(k-1)*eps+n-k);
q = p * eps;
gparams = struct('p', p, 'q', q);

order = 100;

D = [50, 100]; nb_d = numel(D);
P = 0:0.05:1; nb_p = numel(P);
Pert = [0.005, 0.015]; nb_pert = numel(Pert);
nb_simus = 50;

rhos = zeros(nb_pert, nb_simus);
SC_costs = zeros(nb_pert, nb_simus, 2);
d_CSC_costs = zeros(nb_d, nb_pert, nb_p, nb_simus, 2);

%%
for simu = 11:nb_simus
    fprintf('Simu #%d: SC...', simu);
    G1 = gsp_stochastic_block_graph(n, k, gparams);
    G1 = gsp_create_laplacian(G1, 'normalized');
    G1.k = k;

    [U, e] = eigs(G1.L, k, 'sm');
    [e, eid] = sort(diag(e));
    U1 = U(:, eid);
    fprintf('\t\tDone.\n');
    
    for perti = 1:nb_pert
        fprintf('Constructing G2 and SC2 with perti=%d...', perti);
        param_modif = struct('p', p, 'q', q, 'prop_node_mod', Pert(perti), 'prop_edge_mod', Pert(perti), 'until_spc_dec', true);
        G2 = modify_graph(G1, param_modif);
        G2.k= k;
        U2 = G2.U;

        rhos(perti, simu) = norm(U1*U1'-U2*U2', 'fro');
        [ass_SC2, ~, km_SC2] = kmeans(U2, k, 'Replicates', 100, 'Options', statset('UseParallel',1));
        km_SC2 = sum(km_SC2);
        ncut_SC2 = compute_ncut(G2, ass_SC2);
        SC_costs(perti, simu, :) = [km_SC2, ncut_SC2];
        fprintf('\t\tDone.\n');

        for di=1:nb_d
            fprintf('\tsimu=%d, di=%d: CSC G1...', simu, di);
            params_CSC = struct('d_value', D(di), 'features', 1, 'assignment', 1, 'poly_order', order, ...
                                'sampling', 'VD', 'n_factor', 25, 'regu', 1e-4);
            [C_est, F1_est, time_CSC1, comp_terms_CSC1, params_CSC] = my_CSC(G1, params_CSC);
            d = size(F1_est, 2);

            for pi = 1:nb_p
                pp = P(pi);
                if pp == 0
                    if isfield(params_CSC, 'filt_sig'), params_CSC = rmfield(params_CSC, 'filt_sig'); end
                    if isfield(params_CSC, 'lk_est'), params_CSC = rmfield(params_CSC, 'lk_est'); end
                    if isfield(params_CSC, 'weight_VD'), params_CSC = rmfield(params_CSC, 'weight_VD'); end                            
                else
                    params_CSC.filt_sig = F1_est(:, randperm(d, round(pp*d)));
                    params_CSC.weight_VD = comp_terms_CSC1.weight_VD;
                    params_CSC.lk_est = comp_terms_CSC1.lk_est;
                end

                fprintf('\t\tCSC on G2 (pp = %.2f)...\n', pp);
                [ass2w1_CSC, F2w1_est, t_C2w1, comp_terms_CSC2, ~] = my_CSC(G2, params_CSC);

                km_dCSC = compute_kmeans_cost(U2, ass2w1_CSC);
                ncut_dCSC = compute_ncut(G2, ass2w1_CSC);
                d_CSC_costs(di, perti, pi, simu, :) = [km_dCSC, ncut_dCSC];
            end
        end
    end
end
