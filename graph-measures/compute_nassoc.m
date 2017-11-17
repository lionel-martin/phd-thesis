function [ na ] = compute_nassoc( G, labels )
%   COMPUTE_NASSOC Compute the NAssoc cost. The labels must be of size G.N
%   and assign every node to a class.
%
%   compute_nassoc(G, labels) returns the value of the normalized
%   associative cut using G.W for the adjacency matrix and labels for the classes.
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Created by Lionel Martin
    assert(numel(labels) == G.N);

    cats = unique(labels);
    K = length(cats);

    na = 0;

    for k = 1:K
        s = find(labels == cats(k));
        sc = find(labels ~= cats(k));

        na = na + (compute_cut(G, s, s) / compute_volume(G, s)) + (compute_cut(G, sc, sc) / compute_volume(G, sc));
    end
end
