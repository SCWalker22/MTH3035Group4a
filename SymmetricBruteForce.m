
% A = [11, 1, 1, 1;
%     1, 12, 1, 0;
%     1, 1, 13, 0;
%     1, 0, 0, 14];
% Random 4x4 example matrix

% A = [11, 0, 0, 0, 1, 1, 1;
%     0, 12, 0, 0, 1, 0, 1;
%     0, 0, 13, 0, 1, 1, 1;
%     0, 0, 0, 14, 0, 0, 0;
%     1, 1, 1, 0, 15, 0, 0;
%     1, 0, 1, 0, 0, 16, 0;
%     1, 1, 1, 0, 0, 0, 17];
% Random 7x7 example matrix, where fill-in is unavoidable.


A = [21  1  1  1  1  1  1  1;
   1 12  0  0  0  0  0  0;
   1  0 13  0  0  0  0  0;
   1  0  0 14  0  0  0  0;
   1  0  0  0 15  0  0  0;
   1  0  0  0  0 16  0  0;
   1  0  0  0  0  0 17  0;
   1  0  0  0  0  0  0 18];
% This 8x8 matrix is interesting because the Graph representation is like a
% star centered on vertex 1. (the only Edges are "1 to j" for j=2,...,8).

% In the first step, eliminating Vertex 1 first would create a lot of
% fill-in, but eliminating any other Vertex would create no fill-in.

% So fill-in depends on when Vertex 1 is eliminated. There's less overall
% fill-in if Vertex 1 is eliminated later.
% This can be seen in the histogram result, as it goes down in sequence
% 21,15,10,6,3,1,0 ("n choose 2").

% Minimum degree algorithm would prevent this fill-in easily, since it will
% pick each Vertex 2,...,8 (which have degree=1)
% to eliminate before Vertex 1 (which has degree=7 initially)




% check A spd
[A_L, flag] = chol(A,'lower');
if flag+~(issymmetric(A)) >0 
    error('A is not spd.') 
end



n = size(A,1);

p = perms(1:n);

k = factorial(n);

result = zeros(k,n+2); % preallocate result matrix

for i = 1:k
    q = p(i,:);
    A_perm = A(q,q); % only symmetric permutations of A

    [L] = chol(A_perm,'lower'); % single output form (no partial pivoting) of 'chol' 
    m = nnz(L);
    fill_in = m - nnz(tril(A_perm));

    result(i,:) = [i, fill_in, q]; 
end
% matrix of "results" 
% col 1 = index
% col 2 = fill in, 
% col 3+ = permutation of n = order of vertex elimination

sorted_result = sortrows(result, [2,1]);

best = sorted_result(1,:);
worst = sorted_result(k,:);





% super basic histogram of "number of permutations that give x amount of fill-in"
fill_in_vector = sorted_result(:,2);
histogram(fill_in_vector);
ylabel('Number of permutations');
xlabel('Fill-in produced');






% sparsity comparison plot for A, AMD A, RCM A.

amd_perm = amd(sparse(A));
amd_L = chol(A(amd_perm,amd_perm),'lower');

rcm_perm = symrcm(sparse(A));
rcm_L = chol(A(rcm_perm,rcm_perm),'lower');


figure

subplot(2,3,1)
spy(A)
title('Original Matrix A')

subplot(2,3,2) 
spy(A(amd_perm,amd_perm))
title('AMD ordered A')

subplot(2,3,3)
spy(A(rcm_perm,rcm_perm))
title('RCM ordered A')

subplot(2,3,4) 
spy(A_L)
title('Cholesky factor of A')

subplot(2,3,5) 
spy(amd_L)
title('Cholesky factor of AMD ordered A')

subplot(2,3,6) 
spy(rcm_L)
title('Cholesky factor of RCM ordered A')



