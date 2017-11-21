function [ measures ] = compute_proj(create_graph, in_params)
%COMPUTE_PROJ computes the statistics for the projection over Uk.
% This method is used to generates the results of Tables 4.1 and 4.3.

% Copyright (C) 2016 Johan Paratte, Lionel Martin.
% This file is part of the Reproducible Research code base for the Fast
% Eigenspace Approximation using Random Signals (FEARS) method.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

% If you use this code please kindly cite
%     Paratte, Johan and Martin, Lionel. 
%     "Fast Eigenspace Approximation using Random Signals."
%     arXiv preprint arXiv:1611.00938 (2016).
% https://arxiv.org/abs/1611.00938

    measures = {};

    nb_avg = in_params.nb_avg;
    
    mean_proj_exact = zeros(nb_avg, 1);
    mean_proj_fast = zeros(nb_avg, 1);
    mean_proj_std = zeros(nb_avg, 1);
    mean_proj_gt = zeros(nb_avg, 1);
    mean_proj_nosvd = zeros(nb_avg, 1);
    nb_iter_fast = zeros(nb_avg, 1);
    nb_iter_std = zeros(nb_avg, 1);
    k_est_fast = zeros(nb_avg, 1);
    k_est_std = zeros(nb_avg, 1);
    lk_exact = zeros(nb_avg, 1);
    lk_std = zeros(nb_avg, 1);
    lk_fast = zeros(nb_avg, 1);
    
    if ~isfield(in_params, 'verbose')
        verbose = 0;
    else
        verbose = in_params.verbose;
    end
    
    k = in_params.k;

    for nn = 1:in_params.nb_avg

        clear params;
        params.order = in_params.order;
        params.max_calls = 10;
        params.estimate_lk = 0;
        if isfield(in_params, 'filt_filter'), params.filt_filter = in_params.filt_filter; end
        if isfield(in_params, 'lk_filter'), params.lk_filter = in_params.lk_filter; end
        
        ts = tic;
        G = create_graph(in_params.N);
        G = gsp_estimate_lmax(G);
        te = toc(ts);
        if verbose
            fprintf('graph creation : %0.2f s\n',te); 
        end
        
        ts = tic;
        [Uk, ek] = eigs(G.L, in_params.k, 'sm');
        [ek, e_id] =sort(diag(ek));
        Uk = Uk(:, e_id);
        lk = ek(k);
        lk_exact(nn) = lk;
        te = toc(ts);
        if verbose
            fprintf('eigs : %0.2f s\n',te); 
        end
        
        ts = tic;
        [~, flk, ~, ~, nb_iter_lk, k_est_lk] = estimate_lambda_k(G, in_params.k, params);
        te = toc(ts);
        if verbose
            fprintf('est-lk-s : %0.2f s\n',te); 
        end
        nb_iter_std(nn) = nb_iter_lk;
        k_est_std(nn) = k_est_lk;
        lk_std(nn) = flk;
        
        params.lk = flk;

        ts = tic;
        [~, approx_U, ~] = filtering(G, in_params.k, in_params.k, params);
        [Bk, ~, ~] = svd(approx_U, 'econ');
        te = toc(ts);
        if verbose
            fprintf('approx-s : %0.2f s\n',te); 
        end
        mean_proj_std(nn) = mean(sum((Bk(:, 1:k)'*Uk).^2, 2));


        ts = tic;
        [~, flk, nb_iter_lk, k_est_lk] = fast_estimate_lambda_k_rec(G, in_params.k, params);
        te = toc(ts);
        if verbose
            fprintf('est-lk-f : %0.2f s\n',te); 
        end
        nb_iter_fast(nn) = nb_iter_lk;
        k_est_fast(nn) = k_est_lk;
        lk_fast(nn) = flk;
        params.lk = flk;

        ts = tic;
        [~, approx_U, ~] = filtering(G, in_params.k, in_params.k, params);
        [Bk, ~, ~] = svd(approx_U, 'econ');
        te = toc(ts);
        if verbose
            fprintf('approx-f : %0.2f s\n',te); 
        end
        mean_proj_fast(nn) = mean(sum((Bk(:, 1:k)'*Uk).^2, 2));

        
        params.lk = lk;

        ts = tic;
        [~, approx_U, ~] = filtering(G, in_params.k, in_params.k, params);
        [Bk, ~, ~] = svd(approx_U, 'econ');
        te = toc(ts);
        if verbose
            fprintf('approx-ex : %0.2f s\n',te); 
        end
        mean_proj_exact(nn) = mean(sum((Bk(:, 1:k)'*Uk).^2, 2));
        
        ts = tic;
        R = randn(G.N, in_params.k)/sqrt(in_params.k); % the random matrix
        [Bk, ~, ~] = svd(Uk * Uk' * R, 'econ');
        te = toc(ts);
        if verbose
            fprintf('approx-gt : %0.2f s\n',te); 
        end
        mean_proj_gt(nn) = mean(sum((Bk(:, 1:k)'*Uk).^2, 2));


        ts = tic;
        HR = Uk * Uk' * R;
        Bk = HR(:, 1:k);
        Bk = Bk ./ repmat(sqrt(sum(Bk.^2)), G.N, 1);
        te = toc(ts);
        if verbose
            fprintf('approx-nosvd : %0.2f s\n',te); 
        end
        mean_proj_nosvd(nn) = mean(sum((Bk(:, 1:k)'*Uk).^2, 2));

    end
   
    measures.mp = mean(mean_proj_gt);
    measures.mpstd = std(mean_proj_gt);

    measures.mpe = mean(mean_proj_exact);
    measures.mpestd = std(mean_proj_exact);

    measures.mpnosvd = mean(mean_proj_nosvd);
    measures.mpnosvdstd = std(mean_proj_nosvd);

    measures.mps = mean(mean_proj_std);
    measures.mpsstd = std(mean_proj_std);
    
    measures.mpf = mean(mean_proj_fast);
    measures.mpfstd = std(mean_proj_fast);
    
    measures.mif = mean(nb_iter_fast);
    measures.mifstd = std(nb_iter_fast);
    
    %measures.mkf = mean(k_est_fast);
    
    measures.mis = mean(nb_iter_std);
    measures.misstd = std(nb_iter_std);
    
    %measures.mks = mean(k_est_std);
    
    measures.mkf = mean(sqrt((k_est_fast - k).^2));
    measures.mkfstd = std(sqrt((k_est_fast - k).^2));
    
    measures.mks = mean(sqrt((k_est_std - k).^2));
    measures.mksstd = std(sqrt((k_est_std - k).^2));
    
    measures.lkf = mean(sqrt((lk_exact - lk_fast).^2));
    measures.lkfstd = std(sqrt((lk_exact - lk_fast).^2));
    
    measures.lks = mean(sqrt((lk_exact - lk_std).^2));
    measures.lksstd = std(sqrt((lk_exact - lk_std).^2));

    mm = measures;
    
    fprintf('mp=%0.2f, mpe=%0.2f, mpnosvd=%0.2f, mp_f=%0.2f, mp_s=%0.2f, iter_fast=%0.2f, iter_std=%0.2f, dev_k_fast=%0.2f, dev_k_std=%0.2f, dev_lk_fast=%0.2f, dev_lk_std=%0.2f\n', mm.mp, mm.mpe, mm.mpnosvd, mm.mpf, mm.mps, mm.mif, mm.mis, mm.mkf, mm.mks, mm.lkf, mm.lks);

end

