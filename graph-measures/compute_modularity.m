function [ Q ] = compute_modularity( G, labels )
%   COMPUTE_MODULARITY computes the Louvain's modularity for a given graph
%   and partition. 
%
%   compute_modularity(G, part) returns the value of the modularity using
%   G.W for the adjacency matrix, G.d for the nodes' degree and labels for
%   the classes. There must be one label per node of the graph.
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Created by Lionel Martin
    assert(numel(labels)==G.N);

    M = sum(G.d);

    [~, ~, labels] = unique(labels);
    K = max(labels);

    Q = 0;
    for k=1:K
        s = find(labels == k);
        Q = Q + sum(vec(G.W(s, s) - (G.d(s) * G.d(s)' / M)));
    end

    Q = Q / M;
end
