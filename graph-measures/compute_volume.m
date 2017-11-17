function [ v ] = compute_volume( G, idx )
%   COMPUTE_VOLUME computes the volume of a subset of the graph's nodes.
%   The volume is the sum of the degrees of the nodes.
%
%   compute_volume(G, idx) returns the value of the volume using
%   G.d for the nodes' degree and idx for the indices of the nodes to
%   consider.
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Created by Lionel Martin
    v = sum(G.d(idx));
end
