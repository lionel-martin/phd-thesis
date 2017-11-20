%   Experiment of Figure 5.5.
%   Clustering quality using the light clustering method tested on SBM
%   of changing clusterability.
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Created by Lionel Martinnb_simus = 25;

N = 5000;
k = 15;
deg = 40;
nb_simus = 50;

storage_path = '/home/lmartin/Programming/thesis/chap5/measures_exp7.mat';

ec = (deg - sqrt(deg)) / (deg + sqrt(deg) * (k-1));
EPS = 0:0.01:ec; nb_eps = numel(EPS);

MEASURES = {'Time(sec)', 'Modularity', 'Normalized cut'};
nb_measures = numel(MEASURES);

measures = zeros(nb_eps, 2, nb_measures, nb_simus); %third dim for (EBD and SC)

for simu = 1:nb_simus
    fprintf('Simu #%d\n', simu);
    for epsi = 1:nb_eps
        eps = EPS(epsi);
        fprintf('\tGraph complexity %.2f\n', eps);

        sbmparams = struct('deg', deg, 'eps_value', eps, 'lap_type', ...
            'normalized', 'until_spec_dec', true, 'time_spec_dec', true);
        G = generate_sbm(N, k, sbmparams);
        G = gsp_estimate_lmax(G);

        % Construction of the approximated ideal low-pass filter
        param = struct('method', 1, 'order', 100);
        [~, G.lk, ~] = estimate_lambda_k(G, k, param);
        [~, cheb_coef] = jackson_cheby_poly_coefficients(0, G.lk, [0, G.lmax], param.order);
        lowpass_filter = @(x) gsp_cheby_op(G, cheb_coef, x);
        
        fprintf('\t\tLight clustering with EBD\n');
        t_EBD = tic;
        assignment_EBD = ebd_light_clust(G, k, lowpass_filter);
        time_EBD = toc(t_EBD);
        measures(epsi, 1, :, simu) = [time_EBD, compute_modularity(G, assignment_EBD), compute_ncut(G, assignment_EBD)];

        fprintf('\t\tSpectral Clustering\n');

        t_SC = tic;
        assignment_SC = kmeans(G.Uk, k, 'Replicate', 150, 'Options', statset('UseParallel', 1));
        time_SC = toc(t_SC) + G.time_spec_dec;
        measures(epsi, 2, :, simu) = [time_SC, compute_modularity(G, assignment_SC), compute_ncut(G, assignment_SC)];
    end
end

save(storage_path, 'measures');

%% Plotting

% Ploting the +/- 1 sigma deviation
med = prctile(measures, 50, 4);
low = prctile(measures, 16, 4);
hig = prctile(measures, 84, 4);

for measure = 1:3
    figure(measure); hold on;
    for method = 1:2
        errorbar(EPS, med(:, method, measure), med(:, method, measure) - low(:, method, measure), hig(:, method, measure) - med(:, method, measure));
    end
    legend('EBD', 'SC');
    ylabel(MEASURES{measure});
    xlabel('$\varepsilon$', 'Interpreter', 'LaTex');
    set(measure, 'Color', 'w');
    hold off;
end