function [ ch ] = compute_cheeger( G, labels )
%   COMPUTE_CHEEGER computes the Cheeger cut.
%   Labels must be of size G.N and assign every node to a class.
%
%   compute_cheeger(G, labels) returns the value of the Cheeger cut using
%   G.W for the adjacency matrix and labels for the classes.
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Created by Lionel Martin
    cats = unique(labels);
    K = length(cats);

    ch = 0;

    for k = 1:K
        s = find(labels == cats(k));
        sc = find(labels ~= cats(k));

        ch = ch + compute_cut(G, s, sc) / min(compute_volume(G, s), compute_volume(G, sc));
    end
end
