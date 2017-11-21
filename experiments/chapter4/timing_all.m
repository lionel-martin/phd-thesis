function [ data ] = timing_all(G, k, fast)
%TIMING_ALL Function which compute timing for different methods.

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

    %Eigs
    t1 = tic;
    Uk = compute_eigen(G, k);
    tim = toc(t1);
    data.eigs = tim;
    
    %FEARS
    clear params_filt;
    params_filt.order = 50; % order of the polynomial approx
    t1 = tic;
    Bk = gsp_eigenspace_estimation(G, k, params_filt);
    kmeans(Bk, k, 'Replicates', 100, 'Options', setstats('UseParallel', 1));
    tim = toc(t1);
    
    data.our = tim;
    
    %Power method
    if ~fast
        t1 = tic;
        Wk = (speye(G.N) - G.L);
        approx_U = (Wk^k)*randn(G.N, k);
        svd(approx_U, 'econ');
        tim = toc(t1);
        data.power = tim;
    end
    
    % CSC
    param_CSC = struct('poly_order', params_filt.order, 'features', 0, 'assignment', 1, ...
        'lap_type', 'normalized', 'sampling', 'VD', 'n_factor', 0.15*N/(k*log(k)));

    G.k = k;
    t1 = tic;
    my_CSC(G, param_CSC);
    tim = toc(t1);
    
    data.csc = tim;
    
    %CFEARS
    clear params_filt;
    params_filt.order = 50; % order of the polynomial approx

    t1 = tic;
    [Bk, ~, out_params] = gsp_eigenspace_estimation(G, k, params_filt);
    
    param_CSC = struct('poly_order', params_filt.order, 'features', 0, 'assignment', 1, ...
        'lap_type', 'normalized', 'sampling', 'uniform', 'n_factor', 0.15*N/(k*log(k)), ...
        'filt_sig', Bk, 'lk_est', out_params.lk));
    G.k = k;
    my_CSC(G, params_CSC);
    tim = toc(t1);
    
    data.our = tim;
end

