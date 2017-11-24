% Generic code for the evaluation of SC, CSC, dCSC and IASC with different
% parameters. This accepts static and dynamic clustering attempts.
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Modified by Lionel Martin

%% Parametrization
if ~isfield(params, 'NN')
    if ~isfield(params, 'N')
        params.NN = 15000;
    else
        params.NN = [params.N];
    end
end

if ~isfield(params, 'K')
    if ~isfield(params, 'k')
        params.K = 25;
    else
        params.K = [params.k];
    end
end

if ~isfield(params, 'D')
    if ~isfield(params, 'd')
        params.D = 50;
    else
        params.D = [params.d];
    end
end

if ~isfield(params, 'nb_simus'), params.nb_simus = 50; end
if ~isfield(params, 'order'), params.order = 100; end
if ~isfield(params, 'iasc_rank'), params.iasc_rank = 20; end
if ~isfield(params, 'n_prop'), params.n_prop = 0.1; end

if ~isfield(params, 'deg'), params.deg = 60; end
if ~isfield(params, 'EPS_FACTORS'), params.EPS_FACTORS = 0.5; end

if ~isfield(params, 'dynamic'), params.dynamic = true; end
if ~isfield(params, 'METHODS'), params.METHODS = {'SC', 'CSC', 'IASC', 'dCSC'}; end
if ~isfield(params, 'MEASURES'), params.MEASURES = {'kmeans', 'ncut', 'time'}; end

if ~params.dynamic
    params.PNR = 0;
    params.PER = 0;
    params.P = 0;
else
    if ~isfield(params, 'PNR'), params.PNR = 0.01; end
    if ~isfield(params, 'PER'), params.PER = 0.01; end
    if ~isfield(params, 'P'), params.P = 0.1; end

end

nb_N = numel(params.NN);
nb_k = numel(params.K);
nb_d = numel(params.D);
nb_eps = numel(params.EPS_FACTORS);
nb_p = numel(params.P);
nb_pnr = numel(params.PNR);
nb_per = numel(params.PER);
nb_methods = numel(params.METHODS);
nb_measures = numel(params.MEASURES);

measures = zeros(nb_N, nb_k, nb_d, nb_eps, nb_pnr, nb_per, nb_p, nb_methods, nb_measures, nb_simus);
gcp; % start the parallel pool if none is running

%% Running

for simu = 1:nb_simus
    fprintf('Simu n°%d/%d\n', simu, nb_simus);

    for ni = 1:nb_N
        N = NN(ni);

        for ki = 1:nb_k
            k = K(ki);

            for epsi = 1:nb_eps
                eps_factor = EPS_FACTORS(epsi);
                
                fprintf('\tGraph generation...');
                gparams = struct('deg', deg, 'eps_factor', eps_factor, 'lap_type', 'normalized');
                
                if ~params.dynamic && any(strcmp(METHODS, 'SC'))
                    gparams.until_spec_dec = true;
                    gparams.time_spec_dec = true;
                end

                G1 = generate_sbm(N, k, gparams);
                G1.k = k;
                fprintf('\t Done!\n');

                if ~params.dynamic && any(strcmp(METHODS, 'SC'))
                    fprintf('\tSpectral Clustering on G...');
                    if any(strcmp(MEASURES, 'time')), t_SC_KM = tic; end
                    [G1.SC.assignment, ~, kmcost] = kmeans(G1.Uk, k, 'Replicates', 100, 'Options', setstats('UseParallel', 1));
                    fprintf('\tDone!\n');

                    if any(strcmp(MEASURES, 'time'))
                        G1.SC.time.eigen = G1.time_spec_dec;
                        G1.SC.time.kmeans = toc(t_SC_KM);
                        G1.SC.time.total = sum(cell2mat(struct2cell(G1.SC.time)));
                    end
                    
                    if any(strcmp(MEASURES, 'kmeans')), G1.SC.kmeans = sum(kmcost); end
                    if any(strcmp(MEASURES, 'ncut')), G1.SC.ncut = compute_ncut(G1, G1.SC.assignment); end
                end

                if params.dynamic && any(strcmp(METHODS, 'IASC'))
                    fprintf('\tIncremental Approximative Spectral Clustering on G1...');
                    if any(strcmp(MEASURES, 'time')), t_IASC_eigs = tic; end
                    G1.IASC.ABBA = 2*speye(G1.N) - G1.L;
                    [G1.IASC.U, G1.IASC.e] = eigs(G1.IASC.ABBA, max(k, min(params.iasc_rank, G1.N)), 'LM');
                    fprintf('\tDone!\n');

                    if any(strcmp(MEASURES, 'time')), G1.IASC.time.total = toc(t_IASC_eigs); end
                end

                for pnri = 1:nb_pnr
                    pnr = PNR(pnri);

                    for peri = 1:nb_per
                        per = PER(peri);

                        if params.dynamic && pnr == 0 && per == 0
                            warning('Skipping graph modification pnr = 0, per = 0.');
                            continue;
                        end
                        
                        if params.dynamic
                            fprintf('Constructing perturbed graph...');
                            param_modif = struct('p', G1.info.params.p, 'q', G1.info.params.q, ...
                                'prop_node_mod', pnr, 'prop_edge_mod', per);
                            if any(strcmp(METHODS, 'SC')), param_modif.until_spc_dec = true; end
                            G2 = modify_graph(G1, param_modif);
                            G2.k = k;
                            fprintf('\tDone!\n');

                            if any(strcmp(METHODS, 'SC'))
                                fprintf('\tSpectral Clustering on G2...');
                                if any(strcmp(MEASURES, 'time')), t_SC_KM = tic; end
                                [G2.SC.assignment, ~, kmcost] = kmeans(G2.Uk, k, 'Replicates', 100, 'Options', setstats('UseParallel', 1));
                                fprintf('\tDone!\n');

                                if any(strcmp(MEASURES, 'time'))
                                    G2.SC.time.eigen = G2.time_spec_dec;
                                    G2.SC.time.kmeans = toc(t_SC_KM);
                                    G2.SC.time.total = sum(cell2mat(struct2cell(G2.SC.time)));
                                end

                                if any(strcmp(MEASURES, 'kmeans')), G2.SC.kmeans = sum(kmcost); end
                                if any(strcmp(MEASURES, 'ncut')), G2.SC.ncut = compute_ncut(G2, G2.SC.assignment); end
                            end

                            if any(strcmp(METHODS, 'IASC'))
                                fprintf('\tIncremental Approximative Spectral Clustering on G2...');
                                if any(strcmp(MEASURES, 'time')), t_IASC_eigs = tic; end
                                [G2.IASC.e, G2.IASC.U, ~] = eigenUpdate(G2, G1, G1.IASC.e, G1.IASC.U);
                                [G2.IASC.e, e_idx] = sort(diag(G2.IASC.e), 'descend');
                                G2.IASC.U = G2.IASC.U(:, e_idx(1:k));
                                G2.IASC.assignment = kmeans(G2.IASC.U, k, 'Replicates', 100, 'Options', statset('UseParallel', 1));
                                fprintf('\tDone!\n');

                                if any(strcmp(MEASURES, 'time')), G2.IASC.time.total = toc(t_IASC_eigs); end
                                if any(strcmp(MEASURES, 'kmeans')), G2.IASC.kmeans = compute_kmeans_cost(G2.Uk, G2.IASC.assignment); end
                                if any(strcmp(MEASURES, 'ncut')), G2.IASC.ncut = compute_ncut(G2, G2.IASC.assignment); end
                            end
                            
                            G = G2;
                        else
                            G = G1;
                        end

                        for di = 1:nb_d
                            d_val = D(di);
                            
                            if params.dynamic && any(strcmp(METHODS, 'dCSC')) && pnri == 1 && peri == 1  % compute CSC on G1 for dCSC only once.
                                fprintf('\t\tCSC on G1 (d = %d)...\n', d_val);
                                params_CSC = struct('d_value', d_val, 'features', 1, 'assignment', 1, 'poly_order', order, ...
                                    'sampling', 'VD', 'n_factor', params.n_prop * N /(k*log(k)), 'regu', 1e-3);
                                [G1.CSC.assignment, G1.CSC.features, G1.CSC.time, G1.CSC.outputs, G1.CSC.params] = my_CSC(G1, params_CSC);
                            end
                            
                            if any(strcmp(METHODS, 'CSC'))
                                if params.dynamic
                                    fprintf('\t\tCSC on G2 (d = %d)...\n', d_val);
                                else
                                    fprintf('\t\tCSC on G (d = %d)...\n', d_val);
                                end

                                if any(strcmp(MEASURES, 'time')), t_CSC = tic; end

                                params_CSC = struct('d_value', d_val, 'features', 1, 'assignment', 1, 'poly_order', order, ...
                                    'sampling', 'VD', 'n_factor', params.n_prop * N /(k*log(k)), 'regu', 1e-3);
                                [G.CSC.assignment, G.CSC.features, G.CSC.time, G.CSC.outputs, G.CSC.params] = my_CSC(G, params_CSC);
                                
                                if any(strcmp(MEASURES, 'time')), G.CSC.time.total = toc(t_CSC); end
                                if any(strcmp(MEASURES, 'kmeans')), G.CSC.kmeans = compute_kmeans_cost(G.Uk, G.CSC.assignment); end
                                if any(strcmp(MEASURES, 'ncut')), G.CSC.ncut = compute_ncut(G, G.CSC.assignment); end
                            end
                            
                            
                            for pi = 1:nb_P
                                if any(strcmp(METHODS, 'dCSC'))
                                    p = P(pi);
                                    d = size(G1.CSC.features, 2);
                                    params_CSC = G1.CSC.params;
                                    if p == 0
                                        if isfield(params_CSC, 'filt_sig'), params_CSC = rmfield(params_CSC, 'filt_sig'); end
                                        if isfield(params_CSC, 'lk_est'), params_CSC = rmfield(params_CSC, 'lk_est'); end
                                        if isfield(params_CSC, 'weight_VD'), params_CSC = rmfield(params_CSC, 'weight_VD'); end                            
                                    else
                                        params_CSC.filt_sig = G1.CSC.features(:, randperm(d, round(p*d)));
                                        params_CSC.weight_VD = G1.CSC.outputs.weight_VD;
                                        params_CSC.lk_est = G1.CSC.outputs.lk_est;
                                    end

                                    if isfield(params_CSC, 'filt_sig'),
                                       nb_sig = d-size(params_CSC.filt_sig, 2);
                                    else
                                        nb_sig = d;
                                    end
                                    fprintf('\t\tCSC on G2 (d_rem = %d)...\n', nb_sig);

                                    % CSC on G2 with previous signals from CSC of G1
                                    [G.dCSC.assignment, G.dCSC.features, G.dCSC.time, G.dCSC.outputs, ~] = my_CSC(G, params_CSC);
                                end

                                for methi = 1:nb_methods
                                    method = METHODS{methi};
                                    
                                    for measi = 1:nb_measures
                                        measure = MEASURES{measi};
                                        
                                        measures(ni, ki, di, epsi, pnri, peri, pi, methi, measi, simu) = getfield(getfield(G, method), measure);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
