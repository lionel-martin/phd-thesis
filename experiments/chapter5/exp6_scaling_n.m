%   Experiment of Table 5.2.
%   Scalability of our method on sensors of increasing size and comparison.
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Created by Lionel Martin
NN = [1000, 5000, 10000, 20000, 30000, 50000]; nb_N = numel(NN);
k = 10;
nb_simus = 25;

storage_path = '/home/lmartin/Programming/thesis/chap5/measures_exp6.mat';

METHODS = {'EBD', 'SC', 'Newman'};
nb_methods = numel(METHODS);
nb_metrics = 3; % time, Mod, Ncut

measures = zeros(nb_N, nb_methods, nb_metrics, nb_simus);
time_filter = zeros(nb_N, nb_simus);

for simu = 1:nb_simus
    fprintf('Simu #%d\n', simu);

    for ni = 1:nb_N
        N = NN(ni);
        fprintf('\tGraph size: N = %d\n', N);
        
        G = gsp_sensor(N);
        G = gsp_create_laplacian(G, 'normalized');
        G = gsp_estimate_lmax(G);

        param.method = 1;
        param.order = 150;
        t_FM = tic;
        [~, G.lk, ~] = estimate_lambda_k(G, k, param);
        [~, cheb_coef] = jackson_cheby_poly_coefficients(0, G.lk, [0, G.lmax], param.order);
        time_filter(ni, simu) = toc(t_FM);

        for method = 1:nb_methods
            fprintf('\t\t\tMethod: %s\n', METHODS{method});

            t_met = tic;
            switch method
                case 1
                    lowpass_filter = @(x) gsp_cheby_op(G, cheb_coef, x);
                    assignment = ebd_light_clust(G, k, lowpass_filter);
                case 2
                    [Uk, ek] = compute_eigen(G, k);
                    assignment = kmeans(Uk, k, 'Replicate', 150, 'Options', statset('UseParallel', 1));
                case 3
                    assignment = fast_newman(G.W);
            end
            time_method = toc(t_met);
            measures(ni, method, :, simu) = [time_method, compute_modularity(G, assignment), compute_ncut(G, assignment)];
        end
    end
end

save(storage_path, 'measures');

mean(measures, 4)
std(measures, 4)

mean(time_filter, 2)
std(time_filter, 2)