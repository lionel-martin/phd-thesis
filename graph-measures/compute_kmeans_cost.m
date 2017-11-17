function [ kcost ] = compute_kmeans_cost( features, assignment )
%   COMPUTE_KMEANS_COST computes the cost of the kmeans objective function
%   for a given assignment.
%
%   compute_kmeans_cost(features, assignment) returns the value of the
%   objective function of k-means for a set of features. The number of
%   datapoints in the features matrix and the assignment should be the same.
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Created by Lionel Martin
    assert(size(features, 1) == numel(assignment));
    N = numel(assignment);

    [~, ~, assignment] = unique(assignment);
    K = max(assignment);

    X = sparse(1:N, assignment, 1, N, K);
    X = X .* repmat(1./sqrt(sum(X, 1)), N, 1);

    kcost = norm((eye(N) - X * X') * features, 'fro')^2;
end
