function [ c ] = compute_cut( G, idx_a, idx_b )
%   COMPUTE_CUT computes the cut between two subsets of the graph's nodes.
%   A cut corresponds to the sum of weights of the edges going from set A
%   to set B.
%
%   compute_cut(G, idx_a, idx_b) returns the value of the cut using
%   G.W for the adjacency matrix idx_a and idx_b contain the indicies of
%   the nodes for each class.
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Created by Lionel Martin
    assert(numel(intersect(idx_a, idx_b)) == 0);
    assert(numel(idx_a) + numel(idx_b) == G.N);

    c = sum(vec(G.W(idx_a, idx_b)));
end
