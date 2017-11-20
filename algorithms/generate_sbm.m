function [ G ] = generate_sbm( N, k, params )
%   GENERATE_SBM generates SBM graphs either directly from q2 and q1 or from 
%   the average degree and the ratio epsilon = q2 / q1.
%
%   generate_sbm(N, k, params) outputs a graph object based on the
%   parameters that are inputed. Params must either contain q1 and q2 or
%   deg and eps_factor or deg and eps_value. With eps_factor the constant
%   is multiplied by epsilon_max.
%
%   Optional parameter:
%    * 'until_spec_dec' tries to create a SBM until one can be
%    eigendecomposed.
%   * 'time_spec_dec' defines if spectral decomposition should be timed.
%   * 'lap_type' (default 'combinatorial') the Laplacian type.
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Created by Lionel Martin

if nargin < 3, error('Params is mandatory. Read help.'); end
if ~isfield(params, 'until_spec_dec'), params.until_spec_dec = false; end
if ~isfield(params, 'time_spec_dec'), params.time_spec_dec = false; end

methods = struct('with_q', isfield(params, 'q1') && isfield(params, 'q2'), ...
        'with_eps_factor', isfield(params, 'deg') && isfield(params, 'eps_factor'), ...
        'with_eps_value', isfield(params, 'deg') && isfield(params, 'eps_value'));


if ~(methods.with_q || methods.with_eps_factor || methods.with_eps_value)
    error('Params is incomplete. You must either provide q1 and q2, deg and eps_factor or deg and eps_value.');
end

if methods.with_q
    if methods.with_eps_factor || methods.with_eps_value
        warning('Construction run with q1 and q2 although other parameters where also inputed.');
    end
    gparams = struct('p', params.q1, 'q', params.q2);
elseif methods.with_eps_factor
    ec = (params.deg - sqrt(params.deg)) / (params.deg + sqrt(params.deg) * (k-1));
    eps_value = ec * params.eps_factor;
    q1 = params.deg * k / (N * (k-1) * eps_value + N - k);
    q2 = eps_value * q1;
    gparams = struct('p', q1, 'q', q2);
else
    q1 = params.deg * k / (N * (k-1) * params.eps_value + N - k);
    q2 = params.eps_value * q1;
    gparams = struct('p', q1, 'q', q2);
end

if isfield(params, 'lap_type'), gparams.lap_type = params.lap_type; end
G = gsp_stochastic_block_graph(N, k, gparams);

if params.until_spec_dec
    eigen_params = struct('timing', params.time_spec_dec);

    while ~isfield(G, 'Uk')
        try
            [G.Uk, G.ek, G.time_spec_dec] = compute_eigen(G, k, eigen_params);
        catch ME
            warning('Catched %s', ME.identifier);
            G = gsp_stochastic_block_graph(N, k, gparams);
        end
    end
end