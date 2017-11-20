function [ assignment ] = ebd_light_clust( G, k, filter )
%   EBD_LIGHT_CLUST computes the clustering assignment using the
%   "energy-based determination" of centroids and diffusion their labels to
%   the rest of the graph.
%
%   assignment = ebd_light_clust(G, k, filter) outputs the assignment
%   into k classes based on the inputed graph G. The diffusion is performed
%   with the filter described with an anonymous function filter.
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Created by Lionel Martin
    accumulated_energy = zeros(G.N, 1);
    filt_sigs = zeros(G.N, k);

    for ii = 1:k
        [min_tig, ~] = min(accumulated_energy);
        cent_i = datasample(find(accumulated_energy==min_tig), 1); %randomly break ties
        accumulated_energy(cent_i) = NaN;

        delta = zeros(G.N, 1);
        delta(cent_i) = 1;
        filt_sigs(:, ii) = filter(delta);

        accumulated_energy = accumulated_energy + (filt_sigs(:, ii)).^2;
    end

    [~, assignment] = max(filt_sigs, [], 2);
end