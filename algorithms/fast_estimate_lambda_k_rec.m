function [ G, lambda_k, nb_iter, k_est ] = fast_estimate_lambda_k_rec( G, k, param )
% 
% This function estimates the k-th eigenvalue of the Laplacian matrix of a graph G.
% Inputs:
% - G.L should be the sparse symmetrical Laplacian matrix of the graph.
% - G.lmax should be the maximal eigenvalue of G.L.
% - k should be the index of the eigenvalue you want to estimate. 
% - param are parameters. 
% *param.nb_estimation is the number of estimation. Default is 1. 
% *param.nb_features is the number of random features for estimation. Default is 2*round(log(G.N)). 
% *param.epsilon is the precision one desires. Default is 1e-1. 
% *param.hint_lambda_max is a hint given to the algorithm to start the search 
% with param.hint_lambda_max as upper bound. Default is G.lmax.
% *param.hint_c_max is a hint to start the search that indicates the number
% of eigenvalues between 0 and hint_lambda_max. Default is G.N.
% *param.hint_lambda_min is a hint to start the search with that lower
% bound. Default is 0.
% *param.hint_c_min is a hint to start the search with that number of
% eigenvalues between 0 and hint_lambda_min. Default is 0.
% *param.jackson=0 if one wants to use Chebychev polynomials. param.jackson=1 if 
% one wants to use Jackson-Chebychev polynomials. Default is 1.
% *param.order is the order of the polynomial filters used. Default is 50.
% *param.max_calls is the maximum number of iterations of the recursive loop. Default is 10.
% *param.verbose is the verbosity level used. Default is 0.

% Copyright (C) 2016 Johan Paratte, Lionel Martin.
% This file is part of the Reproducible Research code base for the Fast
% Eigenspace Approximation using Random Signals (FEARS) method.


% If you use this code please kindly cite
%     Paratte, Johan, and Martin, Lionel. 
%     "Fast Eigenspace Approximation using Random Signals."
%     arXiv preprint arXiv:1611.00938 (2016).
% https://arxiv.org/abs/1611.00938
%% Parameters
if nargin < 3, param = struct; end
if ~isfield(G, 'lmax'), G = gsp_estimate_lmax(G); end
if ~isfield(param, 'nb_estimation'), param.nb_estimation = 1; end
if ~isfield(param, 'nb_features'), param.nb_features = 2*round(log(G.N)); end
if ~isfield(param, 'epsilon'), param.epsilon = 1e-1; end
if ~isfield(param, 'jackson'), param.jackson = 1; end
if ~isfield(param, 'order'), param.order = 50; end
if ~isfield(param, 'max_calls'), param.max_calls = 10; end
if ~isfield(param, 'verbose'), param.verbose = 0; end

if ~isfield(param, 'hint_lambda_max') || ~isfield(param, 'hint_c_max'),
    param.hint_lambda_max = G.lmax;
    param.hint_c_max = G.N;
end

if ~isfield(param, 'hint_lambda_min') || ~isfield(param, 'hint_c_min'),
    param.hint_lambda_min = 0;
    param.hint_c_min = 0;
end

if ~isfield(param, 'filter')
    if ~isfield(param, 'lk_filter')
        param.filter = 'lp-jch';
    else
        param.filter = param.filt_filter;
    end
end

% List of estimations for lambda_k
lambda_k_est = zeros(param.nb_estimation, 1);

% Perform nb_estimation on different of set feature vectors
for ind_est = 1:param.nb_estimation
    
    % Random signals (fixed for one estimation)
    Sig = randn(G.N, param.nb_features)*1/sqrt(param.nb_features);
    
    % Initial values
    lambda_min = param.hint_lambda_min;
    lambda_max = param.hint_lambda_max;
    cmin = param.hint_c_min;
    cmax = param.hint_c_max;
    nb_calls = 0;
    lambda_est = lambda_min + (k-cmin)*(lambda_max - lambda_min)/(cmax-cmin);

    [G, lk, k_est, nb_iter] = recursive_call(G, k, Sig, lambda_min, lambda_max, lambda_est, lambda_max, cmin, cmax, cmax, nb_calls, param);

    % Store result
    lambda_k_est(ind_est) = lk;
end

% Final estimation
G.lk = mean(lambda_k_est);
lambda_k = G.lk;

end

function [G, lk, count, nb_calls] = recursive_call(G, k, Sig, lmin, lmax, lest, lbest, cmin, cmax, cbest, nb_calls, param)

    nb_calls = nb_calls + 1;

    [c_abs, ~] = compute_eigenvalues(G, lest, Sig, param);

    ck = k - cmin; %local to the recursion
    
    cl = c_abs - cmin; %local to the recursion
    cr = cmax - c_abs; %local to the recursion
    
    if abs(k - c_abs) <= abs(k - cbest)
        cbest = c_abs;
        lbest = lest;
    end
    
    if param.verbose
        fprintf('Iter #%d - lest: %f [%d] - lmin: %f [%d] - lmax: %f [%d] || best: %f [%d]', nb_calls, lest, c_abs, lmin, cmin, lmax, cmax, lbest, cbest);
        if (cl == 0 || cr == 0) fprintf('   [DICHOTOMY]\n'); else fprintf('\n'); end
    end
    
    if cl == ck || cbest == k
       lk = lbest;
       count = cbest;
       return;
    end
    
    %Special cases : c_abs = c_min or c_abs = c_max
    if cl == 0 || cr == 0
        if cl == 0
            lmin = lest;
        elseif cr == 0
            lmax = lest;
        end

        lest = (lmax + lmin) / 2;
    else
        if cl > ck
            lmax = lest;
            lest = lmin + ck * (lest - lmin)/cl;
            cmax = c_abs;
        end

        if cl < ck
            lmin = lest;
            lest = lest + (ck - cl)*(lmax - lest)/cr;
            cmin = c_abs;
        end
    end
    
    if nb_calls == param.max_calls
       lk = lbest;
       count = cbest;
       return;
    end
    
    [G, lk, count, nb_calls] = recursive_call(G, k, Sig, lmin, lmax, lest, lbest, cmin, cmax, cbest, nb_calls, param);
end