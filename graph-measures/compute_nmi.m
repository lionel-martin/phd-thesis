function nmi = compute_nmi(p1, p2)
%   COMPUTE_NMI Compute normalized mutual information between two partitions p1 and p2.
%	I(p1,p2)/sqrt(H(p1)*H(p2)) of two discrete variables p1 and p2.
%	Input:
%	p1, p2: two integer vectors of the same length containing the assignments
%	Ouput:
%	nmi: normalized mutual information nmi=I(p1,p2)/sqrt(H(p1)*H(p2))
%	Written by Mo Chen (sth4nth@gmail.com), modified by Lionel Martin.
assert(numel(p1) == numel(p2));

N = numel(p1);
[~, ~, p1] = unique(p1);
N1 = max(p1);

[~, ~, p2] = unique(p2);
N2 = max(p2);

idx = 1:N;
M1 = sparse(idx, p1, 1, N, N1);
M2 = sparse(idx, p2, 1, N, N2);

% Joint distribution of p1 and p2
P12 = nonzeros(M1' * M2 / N);
H12 = -P12' * log2(P12);

P1 = mean(M1,1);
P2 = mean(M2,1);

% entropy of p1 and p2
H1 = -P1 * log2(P1)';
H2 = -P2 * log2(P2)';

% mutual information
MI = H1 + H2 - H12;

% normalized mutual information
nmi = sqrt((MI / H1) * (MI / H2));
