function [ acc, costM ] = compute_accuracy( p1, p2, masking)
%   COMPUTE_ACCURACY compares two partitioning and compute the accuracy.
%   Accuracy is simply the number of elements that belong to the same 
%   class up to permutation of the classes. Partitions must be of same size.
%
%   compute_accuracy(p1, p2) returns the proportion of elements
%   identically classified between the two partitioning up to labels
%   permutation.
%
%   compute_accuracy(p1, p2, masking) only compare the elements where
%   masking is larger than 0.
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Created by Lionel Martin
    assert(numel(p1) == numel(p2));

    if nargin < 3
        masking = ones(numel(p1), 1);
    end
    assert(numel(p1) == numel(masking));
    p1 = p1(masking>0);
    p2 = p2(masking>0);

    N = numel(p1);

    [~, ~, p1] = unique(p1);
    N1 = max(p1);

    [~, ~, p2] = unique(p2);
    N2 = max(p2);

    costMat = sparse(p1, p2, 1, N1, N2);
    [~, cost] = munkres(-costMat); % we want the minimum

    acc = -cost / N;
    
    if nargout > 1
        costM = costMat;
    end
end
