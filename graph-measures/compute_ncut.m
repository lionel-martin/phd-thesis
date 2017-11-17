function [ nc ] = compute_ncut( G, labels )
%   COMPUTE_NCUT Compute the k-way normalized cut. The labels must be of
%   size G.N and assign every node to a class.
%
%   compute_ncut(G, labels) returns the value of the normalized cut using
%   G.W for the adjacency matrix and labels for the classes.
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Created by Lionel Martin
    assert(numel(labels) == G.N);

    cats = unique(labels);
    K = length(cats);

    nc = 0;

    for k = 1:K
        s = find(labels == cats(k));
        sc = find(labels ~= cats(k));

        nc = nc + compute_cut(G, s, sc) * (1/compute_volume(G, s) + 1/compute_volume(G, sc));
    end
end
