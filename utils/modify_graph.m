function [new_graph] = modify_graph(G, param_modif)
%	MODIFY_GRAPH transforms an SBM into another one by modifying the connectivity
%   param_modif contains the parameters for the modification
%
%   prop_node_mod describes the proportion of nodes to randomly map to
%   another class and reconnect it following the SBM properties
%
%   prop_edge_mod describes the proportion of edges to redraw at random
%   following the belonging of the nodes
%
%   until_spec_dec is a boolean asking for computing until one can compute
%   the spectral decomposition for the new graph (stored in attributes U
%   and e of the structure).
%
%   p and q are the SBM constructive probabilities
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Created by Lionel Martin

assert(isfield(G, 'W'));
assert(isfield(G, 'k'));
assert(isfield(G, 'info') && isfield(G.info, 'node_com'));

assert(isfield(param_modif, 'p'));
assert(isfield(param_modif, 'q'));

if nargin < 2, param_modif = struct; end

if ~isfield(param_modif, 'prop_node_mod'), param_modif.prop_node_mod = 0; end
if ~isfield(param_modif, 'prop_edge_mod'), param_modif.prop_edge_mod = 0; end
if ~isfield(param_modif, 'until_spc_dec'), param_modif.until_spc_dec = false; end

infeasible = true;
while infeasible
    W_init = G.W;
    N = G.N;
    K = G.k;
    coms_init = G.info.node_com;

    %% NODE MOVEMENT
    nb_moved = round(param_modif.prop_node_mod * N);
    mov_idx = randperm(N, nb_moved);

    % Deconnecting the selected nodes from their old neighbors
    W_new = W_init;
    W_new(mov_idx, :) = 0;
    W_new(:, mov_idx) = 0;

    % Creating the new community list and putting zeros in the old one at nodes moved
    coms_new = coms_init;
    cand_new_com = repmat(1:K, nb_moved, 1)';
    cand_new_com(coms_new(mov_idx) + [0:K:K*(nb_moved-1)]) = [];
    coms_new(mov_idx) = cand_new_com(randi(K-1, 1, nb_moved) + [0:K-1:(K-1)*(nb_moved-1)]);
    coms_init(mov_idx) = 0;

    for ii=1:K
        % Reconnecting moved nodes within their new community (symmetrically)
        rect_sparse = gsp_sparse_mat(numel(find(coms_new(mov_idx)==ii)), numel(find(coms_init==ii)), param_modif.p);
        W_new(mov_idx(coms_new(mov_idx)==ii), coms_init==ii) = rect_sparse;
        W_new(coms_init==ii, mov_idx(coms_new(mov_idx)==ii)) = rect_sparse';

        for jj=1:K
            % Reconnecting moved nodes with other moved nodes
            if ii < jj
                W_new(mov_idx(coms_new(mov_idx)==ii), mov_idx(coms_new(mov_idx)==jj)) = gsp_sparse_mat(numel(find(coms_new(mov_idx)==ii)), numel(find(coms_new(mov_idx)==jj)), param_modif.q);
            elseif ii == jj
                W_new(mov_idx(coms_new(mov_idx)==ii), mov_idx(coms_new(mov_idx)==ii)) = gsp_sparse_mat_sym(numel(find(coms_new(mov_idx)==ii)), param_modif.p);
            else %(use symmetry)
                W_new(mov_idx(coms_new(mov_idx)==jj), mov_idx(coms_new(mov_idx)==ii)) = W_new(mov_idx(coms_new(mov_idx)==ii), mov_idx(coms_new(mov_idx)==jj))';
            end
            % Reconnecting moved nodes with the non moved nodes of other clusters
            rect_sparse = gsp_sparse_mat(numel(find(coms_init==ii)), numel(find(coms_new(mov_idx)==jj)), param_modif.q);
            W_new(coms_init==ii, mov_idx(coms_new(mov_idx)==jj)) = rect_sparse;
            W_new(mov_idx(coms_new(mov_idx)==jj), coms_init==ii) = rect_sparse';
        end
    end

    W_new(1:N+1:end) = 0;

    %% EDGE MOVEMENT
    edges = find(triu(W_new));
    M = numel(edges);
    not_edges = find(triu(ones(N), 1));
    not_edges = setdiff(not_edges, edges);
    nb_moved = round(param_modif.prop_edge_mod * M);

    %Deletion of edges (symmetrical)
    del_idx = edges(randperm(numel(edges), nb_moved));
    [I_rem, J_rem] = ind2sub([N, N], del_idx);
    W_new(sub2ind([N, N], [I_rem, J_rem], [J_rem, I_rem])) = 0;

    %Addition of edges
    weights = ones(size(not_edges)) * param_modif.q;
    class_diff = repmat(coms_new, N, 1) - repmat(coms_new', 1, N);   
    class_diff = class_diff(not_edges);
    weights(class_diff == 0) = param_modif.p;

    new_edges = datasample(not_edges, nb_moved, 'Replace', false, 'Weights', weights);
    [Is, Js] = ind2sub([N, N], new_edges);
    W_new = W_new + sparse([Is, Js], [Js, Is], 1, N, N);

    new_graph = gsp_graph(W_new);
    new_graph.info.node_com = coms_new;
    new_graph = gsp_create_laplacian(new_graph, 'normalized');
    if param_modif.until_spec_dec
        try
            tSC = tic;
            [U, e] = eigs(new_graph.L, G.k, 'sm');
            [new_graph.ek, feat_idx] = sort(diag(e));
            new_graph.Uk = U(:, feat_idx);
            infeasible = 0;
        catch ME
            infeasible = 1;
            warning('Catched %s', ME.identifier);
        end

        new_graph.time_spec_dec = toc(tSC);

    else
        infeasible = false;
    end
end

