% Modified version of CSC for personal use.
% This version supports reusing signals as required by the dynCSC algorithm.
% Original work by Tremblay et al. is available in the GSPBox.
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Modified by Lionel Martin
function [assignment, features, time, computed_terms, param_CSC] = my_CSC(G, param_CSC)

%%%
%% ====================== check parameter list ================================
%%%

if nargin == 1, param_CSC = struct; end

if ~isfield(param_CSC, 'poly_order'), param_CSC.poly_order = 50; end
if ~isfield(param_CSC, 'regu'), param_CSC.regu = 1e-3; end
if ~isfield(param_CSC, 'sampling'), param_CSC.sampling = 'uniform'; end

if ~isfield(param_CSC, 'n_factor'),
    % n = param_CSC.n_factor * G.k * log(G.k);
    param_CSC.n_factor = 2;
end

if ~isfield(param_CSC, 'd_factor'),
    % d = param_CSC.d_factor * log(n);
    param_CSC.d_factor = 4;
end

if ~isfield(param_CSC, 'lap_type'), param_CSC.lap_type = 'normalized'; end
if ~isfield(param_CSC, 'normBOOL'), param_CSC.normBOOL = strcmp(param_CSC.lap_type,'normalized'); end
if ~isfield(param_CSC, 'solver'), param_CSC.solver = 'gmres'; end
if ~isfield(param_CSC, 'features'), param_CSC.features = 0; end
if ~isfield(param_CSC, 'assignment'), param_CSC.assignment = 1; end
if ~isfield(param_CSC, 'filt_sig'), param_CSC.filt_sig = zeros(G.N, 0); end
if size(param_CSC.filt_sig, 1) ~= G.N, param_CSC.filt_sig = zeros(G.N, 0); end

n=round(param_CSC.n_factor*G.k*log(G.k));

if isfield(param_CSC, 'd_value'),
    d = param_CSC.d_value;
elseif param_CSC.assignment == 1,
    d=round(param_CSC.d_factor*log(n));
else
    d=round(param_CSC.d_factor*log(G.N));
end

%%%
%% ==================== compute required Laplacian matrix ===========================
%%%

if strcmp(G.lap_type, param_CSC.lap_type) == 0
    G = gsp_create_laplacian(G,param_CSC.lap_type);
end

if size(param_CSC.filt_sig, 2) > d,
    param_CSC.filt_sig = param_CSC.filt_sig(:, 1:d);
end

d_rem = d - size(param_CSC.filt_sig, 2);

X_lk_est = [];
%%%
%% ====================== estimate lambda_k and weight VD ===========================
%%%
if d_rem > 0
    if param_CSC.normBOOL
        G.lmax=2;
        time.lmax=0;
    else 
        tic;
        opts.isreal=1;opts.issym=1;opts.maxit=10000;
        G.lmax=eigs(G.L,1,'LA',opts);
        time.lmax=toc;
    end

    fprintf('\n\nEstimating lambda_k...')
    tic;
    param.order=param_CSC.poly_order;
    [~, lk_est, cum_coh_k, ~] = estimate_lambda_k(G, G.k, param);

    param.hint_lambda_max=lk_est*2;   
    [~, lk_estp1, cum_coh_k_p1, ~] = estimate_lambda_k(G, G.k, param);

    lk_est=(lk_est+lk_estp1)/2;
    time.lk_est=toc;

    mean_num_coh=mean([cum_coh_k,cum_coh_k_p1],2);
    weight_VD = sum(mean_num_coh,2); 
    weight_VD=weight_VD./sum(weight_VD); 
    fprintf('\t\t\tDone.\n')

%%%
%% ====================== filter d random vectors ===========================
%%%

    G.lk=lk_est;

    fprintf('Filtering random signals...')

    tic;
    R=(1/sqrt(d)).*randn(G.N,d_rem);

    [~,JCH] = jackson_cheby_poly_coefficients(0,G.lk,[0,G.lmax],param_CSC.poly_order);
    X_lk_est = gsp_cheby_op(G, JCH, R);

    % if normBOOL
    %     X_lk_est=X_lk_est./repmat(sqrt(sum(X_lk_est.^2,2)),1,d);
    % end

    time.filtering=toc;
    fprintf('\t\tDone.\n')

else
    param_CSC.sampling = 'uniform';
    lk_est = param_CSC.lk_est;
end

if d_rem < d
    warning('Combining with input filtered signals');
    X_lk_est = [X_lk_est, param_CSC.filt_sig];
end

features = [];
if param_CSC.features
    features = X_lk_est;
    assignment = [];

    if d_rem > 0
        computed_terms = struct('lk_est', lk_est, 'weight_VD', weight_VD, 'R', R, 'JCH', JCH, 'G', G);
    else
        computed_terms = struct();
    end
end


%%%
%% ====================== downsample n nodes ===========================
%%%
if param_CSC.assignment == 1

    if strcmp(param_CSC.sampling, 'uniform')
        weight = ones(G.N, 1)/(G.N); % Uniform density
    elseif strcmp(param_CSC.sampling, 'VD')
        weight = weight_VD; % Variable density
    else
        error(['CSC: param_CSC.sampling must be either set to ''uniform'' or to ''VD''']);
    end

    ind_obs = datasample(1:G.N, n, 'Replace', false, 'Weights', weight);
    X_lk_est_DS = X_lk_est(ind_obs, :);

    %%%
    %% ====================== do k-means in low dimension ===========================
    %%%

    fprintf('Low-dimensional kmeans...')
    tic; 
    IDX_LD = kmeans(X_lk_est_DS , G.k, 'Replicates', 100, 'Options', setstats('UseParallel', 1));
    time.k_means_low_dim=toc;
    fprintf('\t\tDone.\n')

    %%%
    %% ====================== Interpolate in high dimensions: ===========================
    %%%

    fprintf('Interpolation of cluster indicators...\n\n')
    tic;
    C_obs_LD = sparse(1:n, IDX_LD, 1, n, G.k);
    [~,JCH_HP] = jackson_cheby_poly_coefficients(lk_est, G.lmax, [0, G.lmax], param_CSC.poly_order);

    C_est = zeros(G.N, size(C_obs_LD, 2));
    parfor k=1:size(C_obs_LD, 2)
       c_obs = C_obs_LD(:, k);
       C_est(:, k) = interpolate_on_complete_graph(c_obs, ind_obs, @(x)gsp_cheby_op(G, JCH_HP, x), param_CSC.regu, G.N, param_CSC.solver);
    end

    [~, assignment] = max(C_est./repmat(sqrt(sum(C_est.^2, 1)), G.N, 1), [], 2);
    time.interpolation=toc;

    time.total = sum(cell2mat(struct2cell(time)));

    if d_rem > 0
        computed_terms = struct('lk_est', lk_est, 'IDX_LD', IDX_LD, 'ind_obs', ind_obs, 'weight_VD', weight_VD, 'R', R, 'JCH', JCH, 'G', G);
    else
        computed_terms = struct('IDX_LD', IDX_LD, 'ind_obs', ind_obs, 'G', G);
    end
end
