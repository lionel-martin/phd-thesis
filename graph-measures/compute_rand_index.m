function [ ri ] = compute_rand_index( p1, p2, varargin )
%   COMPUTE_RAND_INDEX Compute the rand index between two partitions.
%
%   compute_rand_index(p1, p2) computes the rand index between partitions p1 and p2.
%
%   compute_rand_index(p1, p2, 'adjusted') computes the adjusted rand index
%   between partitions p1 and p2. The adjustment accounts for chance
%   correlation.

    assert(numel(p1) == numel(p2));

    % Parse the input and throw errors
    adj = false;

    if nargin > 3
        error('Too many input arguments');
    end

    if nargin == 3
        if strcmp(varargin{1}, 'adjusted')
            adj = true;
        else
            error('%s is an unrecognized argument.', varargin{1});
        end
    end

    % Preliminary computations and cleansing of the partitions
    N = length(p1);
    [~, ~, p1] = unique(p1);
    N1 = max(p1);
    [~, ~, p2] = unique(p2);
    N2 = max(p2);

    % Create the matching matrix
    M1 = sparse(1:N, p1, 1, N, N1);
    M2 = sparse(1:N, p2, 1, N, N2);
    n = M1' * M2;
    
    % Calculate the chosen rand index
    switch adj
        case 0
            ss = sum(vec(n.^2));
            ss1 = sum(sum(n, 1).^2);
            ss2 = sum(sum(n, 2).^2);
            if N < 2
                ri = 0;
            else
                ri = 1 + (2*ss - ss1 - ss2)/(N * (N - 1));
            end

        case 1
            % sum of n_ij choose 2
            ssm = sum(vec(n) .* vec(n - 1)) / 2;

            % sum of n_i choose 2
            ni = sum(n, 2);
            sm1 = sum(vec(ni) .* vec(ni - 1)) / 2;
            
            % sum of n_j choose 2
            nj = sum(n, 1);
            sm2 = sum(vec(nj) .* vec(nj - 1)) / 2;

            sm = N * (N - 1) / 2;

            esm = (sm1 * sm2 / sm);
            num = ssm - esm;
            den = (sm1 + sm2)/2 - esm;
            ri = num / den;
    end 
end
