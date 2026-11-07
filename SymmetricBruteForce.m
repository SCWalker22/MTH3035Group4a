

A = [21  1  1  1  1  1  1  1;
    1 12  0  0  0  0  0  0;
    1  0 13  0  0  0  0  0;
    1  0  0 14  0  0  0  0;
    1  0  0  0 15  0  0  0;
    1  0  0  0  0 16  0  0;
    1  0  0  0  0  0 17  0;
    1  0  0  0  0  0  0 18];
% This matrix is interesting because the Graph representation is like a
% star centered on vertex 1. (the only Edges are "1 to j" for j=2,...,8).

% In the first step, eliminating Vertex 1 first would create a lot of
% fill-in, but eliminating any other Vertex would create no fill-in.

% So fill-in depends on when Vertex 1 is eliminated. There's less overall
% fill-in if Vertex 1 is eliminated later.

% Minimum degree algorithm would prevent this fill-in easily, since it will
% pick each Vertex 2,...,8 (which have degree=1)
% to eliminate before Vertex 1 (which has degree=7 initially)

n = size(A,1);

p = perms(1:n);

k = factorial(n);

result = zeros(k,n+2); % preallocate result matrix

for i = 1:k
    q = p(i,:);
    A_perm = A(q,q); % only symmetric permutations of A

    [L,flag] = chol(A_perm,'lower'); % double output form (no partial pivoting) of 'chol' 
    m = nnz(L);
    fill_in = m - nnz(tril(A_perm));

    result(i,:) = [i, fill_in, q]; 
    % matrix of "results" 
    % col 1 = index
    % col 2 = fill in, 
    % col 3+ = permutation of n = order of vertex elimination
end

sorted_result = sortrows(result, [2,1])

best = sorted_result(1,:)
worst = sorted_result(k,:)





flag % if this is non-zero, then the original matrix isnt spd