function [ assignment ] = as_light_clust( G, k, filter, d, params )
%   AS_LIGHT_CLUST computes the clustering assignment using the
%   "agglomerative selection" of centroids and diffusion their labels to
%   the rest of the graph.
%
%   assignment = as_light_clust(G, k, filter) outputs the assignment
%   into k classes based on the inputed graph G. The diffusion is performed
%   with the filter described with an anonymous function filter.
%
%   assignment = as_light_clust(G, k, filter, d) set the number of
%   candidates to d (default: 2klog(k)).
%
%   assignment = as_light_clust(G, k, filter, params) allows to compute
%   the all-nodes norm Tig with the approximated method of Perraudin.
%   Set params.tig_approx to true and define the number of random Bernoulli
%   signals with params.nb_sigs (default: 2klog(k)).
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Created by Lionel Martin

    if nargin < 4, d = round(2 * k * log(k)); end
    if nargin < 5, params = struct(); end
    if ~isfield(params, 'tig_approx'), params.tig_approx = G.N >= 1000; end
    if ~isfield(params, 'nb_sigs'), params.nb_sigs = round(2*k*log(k)); end
    
    if ~params.tig_approx && G.N > 5000
        warning('[AS_LIGHT_CLUST] All Tig norms will be computed for a large graph. Pass tig_approx=true to estimate the norms instead.');
    end

    if params.tig_approx
        R = 2 * randi(2, G.N, params.nb_sigs) - 3;
        filt_tig_norm = sum((filter(R)).^2, 2) / params.nb_sigs;
    else
        filt_sigs = filter(eye(G.N));
        filt_sigs_res = reshape(filt_sigs, G.N, [], G.N);
        filt_tig_norm = squeeze(sum(filt_sigs_res.^2, 1));
    end

    centers = datasample(1:G.N, d, 'Replace', false, 'Weights', 1./filt_tig_norm);
    deltas = sparse(centers, 1:d, 1, G.N, d);
    data = filter(deltas);

    center2center = data(centers, :);
    center2center(1:d+1:end) = 0;
    final_centers = centers;

    while numel(final_centers) > k
        [~, min_idx] = min(vec(center2center));
        
        % row and col are the two closest centers
        [row, col] = ind2sub(size(center2center), min_idx);
        center2center(:, row) = []; % remove the diffusion of row
        center2center(row, :) = []; % remove row as a potential center
        final_centers(row) = []; % remove rows from the centers list
        data(:, row) = []; % remove the diffusion from row
    end

    [~, assignment] = max(data, [], 2);
end