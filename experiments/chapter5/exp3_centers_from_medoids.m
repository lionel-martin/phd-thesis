%   Experiment of Table 5.1.
%   Replacing the assignment of k-medoids with that of the diffusion of
%   Kronecker positioned at the centers of k-medoids.
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Created by Lionel Martin
N = 10000;
k = 15;
deg = 60;

verbose = true;
storage_path = '';
nb_simus = 100;
nb_graphs = 2;      % sensor(N), sbm(N, k, deg, ec/2)
nb_measures = 4;    % mod exact, ncut exact, mod approx, ncut approx
measures = zeros(nb_simus, nb_graphs, nb_measures);

for simu = 1:nb_simus
    if verbose, fprintf('Simu #%d\n', simu); end
    for gid = 1:nb_graphs
        if gid == 1
            G = gsp_sensor(N);
        else
            sbmparams = struct('deg', deg, 'eps_factor', 0.5);
            generate_sbm(N, k, sbmparams)
        end

        G = gsp_create_laplacian(G, 'normalized');
        [Uk, ek] = compute_eigen(G, k);
        [ass, ~, ~, dist_ex, centers] = kmedoids(Uk, k, 'Replicates', 50);
        cent_sigs = zeros(G.N, k);
        cent_sigs(sub2ind([G.N, k], centers', 1:k)) = 1;
        dist_apx = Uk*Uk'*cent_sigs; % exact filtering with the ideal low-pass filter
        [~, ass_apx] = max(dist_apx, [], 2);

        measures(simu, gid, :) = [compute_modularity(G, ass), compute_ncut(G, ass), compute_modularity(G, ass_apx), compute_ncut(G, ass_apx)];
    end
end

if ~strcmp(storage_path, '')
    save(sprintf('%s/exp3_measures.mat', storage_path), 'measures');
end

if verbose, fprintf('Done!\n'); end

% Data analysis
mean(measures, 1)
std(measures, 1)
