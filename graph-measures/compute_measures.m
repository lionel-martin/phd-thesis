function [ cuts, measures ] = compute_measures( G, labels, measures )
%   COMPUTE_CUTS computes different cuts based on the arguments.
%
%   [cuts, measures] = compute_measures(G, labels, measures) returns the value
%   of the selected measures. There must be one label per node of the graph.
%   Allowed measures are "cheeger", "ncut", "nassoc".
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Created by Lionel Martin
    assert(numel(labels) == G.N);

    cats = unique(labels);
    K = length(cats);

    measures = unique(measures);
    measures = measures(union(union(strmatch('cheeger', measures), strmatch('ncut', measures)), strmatch('nassoc', measures)));
    M = length(measures);

    cuts = zeros(1,M);

    for k = 1:K
        s = find(labels == cats(k));
        sc = find(labels ~= cats(k));

        for m = 1:M
            cuts(strmatch('cheeger', measures)) = cuts(strmatch('cheeger', measures)) + compute_cut(G, s, sc) / min(compute_volume(G, s), compute_volume(G, sc));
            cuts(strmatch('nassoc', measures)) = cuts(strmatch('nassoc', measures)) + (compute_cut(G, s, s) / compute_volume(G, s)) + (compute_cut(G, sc, sc) / compute_volume(G, sc));
            cuts(strmatch('ncut', measures)) = cuts(strmatch('ncut', measures)) + compute_cut(G, s, sc) * (1/compute_volume(G, s) + 1/compute_volume(G, sc));
    end
end
