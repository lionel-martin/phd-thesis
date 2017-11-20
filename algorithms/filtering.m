function [ G, approx_U, basis, params ] = filtering( G, k, d, params )
%   FILTERING approximates U_k
%   G: the graph
%   k: the dimension of the subspace
%   d: the number of signals to filter
%   params: a struct of parameters
%       - is_exact: boolean for exact low-pass filtering or cheby approx
%       - order: the order for filtering with jackson-cheby
%       - jch: the jackson-chebyshev coefficients if they are already known
%       - lk: the estimate of lambda_k if already known
%       - verbose: the quantity of log displayed

    if nargin < 4, params = struct; end
    if ~isfield(params, 'is_exact'), params.is_exact = 0; end
    if ~isfield(params, 'order'), params.order = 50; end
    if ~isfield(params, 'verbose'), params.verbose = 0; end
    if ~isfield(params, 'fast_lk'), params.fast_lk = 1; end
    if ~isfield(params, 'estimate_lk'), params.estimate_lk = 1; end
    if ~isfield(params, 'filter')
        if ~isfield(params, 'filt_filter')
            params.filter = 'lp-jch';
        else
            params.filter = params.filt_filter;
        end
    end

    if ~isfield(G, 'U') && params.is_exact
        tic;
        G = gsp_compute_fourier_basis(G);
        t = toc;
        if params.verbose, disp(['* Time to compute fourier basis: ', num2str(t)]); end
    end



    if ~isfield(params, 'jch') && ~params.is_exact
        if ~isfield(G, 'lmax')
            G = gsp_estimate_lmax(G);
        end

        if params.estimate_lk || ~isfield(params, 'lk')
            tic;
            if params.fast_lk
                [~, params.lk, nb_iter_lk, k_est_lk] = fast_estimate_lambda_k_rec(G, k, params);
            else
                [~, params.lk, ~, ~, nb_iter_lk, k_est_lk] = estimate_lambda_k(G, k, params);
            end

            t = toc;
            if params.verbose, disp(['* Estimated lk: ', num2str(params.lk)]); end
            if params.verbose, disp(['* Time to estimate lk: ', num2str(t)]); end
            if params.verbose, fprintf('* in %d iterations with k_est=%d (target=%d) \n', nb_iter_lk, k_est_lk, k); end
        end

        tic;
        switch params.filter
            case 'lp-ch'
                [params.pcoefs, ~] = jackson_cheby_poly_coefficients(0, params.lk, [0, G.lmax], params.order);

            case 'lp-jch'
                [~, params.pcoefs] = jackson_cheby_poly_coefficients(0, params.lk, [0, G.lmax], params.order);

            case 'expwin'
                ew = gsp_design_expwin(G, params.lk/G.lmax);
                params.pcoefs = gsp_cheby_coeff(G, ew, params.order);

            case 'linear'
                lin = @(x) 1000*k*(1-x/params.lk) .* (x < params.lk);
                params.pcoefs = gsp_cheby_coeff(G, lin, params.order);

            otherwise
                disp('unknown filter type');     
        end

        t = toc;
        if params.verbose, gsp_plot_filter(G, @(x) gsp_cheby_eval(x, params.pcoefs, [0, G.lmax])); end
        if params.verbose, disp(['* Time to compute jch coeffs: ', num2str(t)]); end
    end


    tic;
    R = norm(G.N, d)/sqrt(d);

    if params.is_exact

        %Vcum = zeros(G.N,1);

        %for i=1:G.N
        %    delta = zeros(G.N, 1);
        %    delta(i) = 1;
        %    Vcum(i) = 1/sqrt(norm(G.U(:, 1:k)'*delta));
        %end

        approx_U = G.U(:, 1:k) * G.U(:, 1:k)' * R;
    else
        approx_U = gsp_cheby_op(G, params.pcoefs, R);
    end
    t = toc;
    if params.verbose, disp(['* Time to filter random signals: ', num2str(t)]); end

    tic;
    [basis, ~, ~] = svd(approx_U, 'econ');

    t = toc;
    if params.verbose, disp(['* Time to compute SVD: ', num2str(t)]); end
end
