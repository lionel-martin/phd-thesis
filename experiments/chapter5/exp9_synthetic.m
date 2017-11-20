%   Experiment of Tables 5.4 - 5.9.
%   Comparison of the two heuristics for centroid determination and of the 
%   different graph filters (each with two parameters) on various graphs.
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Created by Lionel Martin

nb_simus = 50;
k = 8;

storage_path = '/home/lmartin/Programming/thesis/chap5/';

GRAPHS = {'Sensor 1k', 'Sensor 10k', 'SBM 1k', 'SBM 10k', 'Bunny', 'Barbara'};
nb_graphs = numel(GRAPHS);

FILTERS = {'Tik 1e-2', 'Tik 1e-1', 'Heat 10', 'Heat 50', 'Step 50', 'Step 150'};
nb_filters = numel(FILTERS);

GTS = {'Spectral Clustering'};
nb_gt = numel(GTS);

METHODS = {'Energy Based Determination', 'Agllomerative Selection'};
nb_methods = numel(METHODS);

MEASURES = {'Time(sec)', 'Modularity', 'Normalized cut'};
nb_measures = numel(MEASURES);

filt_tik = @(x, gamma) 1 ./ (1+gamma * x);

groundtruths = zeros(nb_graphs, nb_gt, nb_measures, nb_simus);
measures = zeros(nb_graphs, nb_filters, nb_methods, nb_measures, nb_simus);

gcp; % start parpool if not active

for simu = 1:nb_simus
    fprintf('Simu #%d\n', simu);

    for gri = [1, 3, 5, 6]
        fprintf('\tGraph: %s\n', GRAPHS{gri});
        switch gri
            case 1
                N = 1000;
                G = gsp_sensor(N);
            case 2
                N = 10000;
                G = gsp_sensor(N);
            case 3
                N = 1000; deg = 30;
                gparams = struct('deg', deg, 'eps_factor', 0.25);
                G = generate_sbm(N, k, gparams);
            case 4
                N = 10000; deg = 60;
                gparams = struct('deg', deg, 'eps_factor', 0.25);
                G = generate_sbm(N, k, gparams);
            case 5
                G = gsp_bunny();
                N = G.N;
            case 6
                img = imread('/../data/barbara.png');
                img = imresize(img, 0.125);
                parampatch.rho = 50;
                G = gsp_patch_graph(img, parampatch);
                N = G.N;
        end
        
        G = gsp_create_laplacian(G, 'normalized');
        G = gsp_estimate_lmax(G);

        for filt = 1:nb_filters
            fprintf('\t\tFilter: %s\n', FILTERS{filt});
            clear g cheb_coef;
            tic;
            switch filt
                case 1
                    g = @(x) filt_tik(x, 1e-2);
                case 2
                    g = @(x) filt_tik(x, 1e-1);
                case 3
                    g = gsp_design_heat(G, 10);
                case 4
                    g = gsp_design_heat(G, 50);
                case 5
                    param = struct('method', 1, 'order', 50);
                    [~, G.lk, ~] = estimate_lambda_k(G, k, param);
                    [~, cheb_coef] = jackson_cheby_poly_coefficients(0, G.lk, [0, G.lmax], param.order);
                case 6
                    param = struct('method', 1, 'order', 150);
                    [~, G.lk, ~] = estimate_lambda_k(G, k, param);
                    [~, cheb_coef] = jackson_cheby_poly_coefficients(0, G.lk, [0, G.lmax], param.order);
            end
            
            if ~exist('cheb_coef', 'var'), cheb_coef = gsp_cheby_coeff(G, g, 100, 101); end
            filter = @(x) gsp_cheby_op(G, cheb_coef, x);
            
            for method = 1:nb_methods
                fprintf('\t\t\tMethod: %s\n', METHODS{method});

                tic;
                switch method
                    case 1
                        assignment = ebd_light_clust(G, k, filter);
                    case 2
                        assignment = as_light_clust(G, k, filter);
                end
                time_method = toc;
                
                measures(gri, filt, method, :, simu) = [time_method, compute_modularity(G, assignment), compute_ncut(G, assignment)];
            end
        end
        
        for gti = 1:nb_gt
            fprintf('\t\t %s\n', GTS{gti});
            tic;
            switch gti
                case 1
                    [Uk, ek] = compute_eigen(G, k);
                    assignment = kmeans(Uk, k, 'Replicate', 150, 'Options', statset('UseParallel', 1));
                case 2
                    assignment = fast_newman(G.W);
            end
            time_gt = toc;
        end

        groundtruths(gri, gti, :, simu) = [time_gt, compute_modularity(G, assignment), compute_ncut(G, assignment)];
    end
end

save(sprintf('%s/exp8_measures.mat', storage_path), 'measures');
save(sprintf('%s/exp8_groundtruths.mat', storage_path), 'groundtruths');

mean(measures, 5)
std(measures, 5)