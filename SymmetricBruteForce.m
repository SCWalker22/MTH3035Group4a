
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

A = [11, 0, 0, 0, 0, 1, 1;
    0, 12, 1, 1, 1, 0, 1;
    0, 1, 13, 0, 0, 0, 1;
    0, 1, 0, 14, 0, 0, 1;
    0, 1, 0, 0, 15, 0, 0;
    1, 0, 0, 0, 0, 16, 0;
    1, 1, 1, 1, 0, 0, 17];
% Random 7x7 example matrix, where AMD<RCM (fill-in).
% *Its hard to get an example with RCM<AMB.

% A = [21  1  1  1  1  1  1  1;
%    1 12  0  0  0  0  0  0;
%    1  0 13  0  0  0  0  0;
%    1  0  0 14  0  0  0  0;
%    1  0  0  0 15  0  0  0;
%    1  0  0  0  0 16  0  0;
%    1  0  0  0  0  0 17  0;
%    1  0  0  0  0  0  0 18];
% % This 8x8 matrix is interesting because the graph representation is like a
% star centered on Vertex 1. (the only Edges are "1 to j" for j=2,...,8).

% In the first step, eliminating Vertex 1 first would create a lot of
% fill-in, but eliminating any other Vertex would create no fill-in.

% So fill-in depends on when Vertex 1 is eliminated. There's less overall
% fill-in if Vertex 1 is eliminated later.
% This can be seen in the histogram result, as it goes down in sequence
% 21,15,10,6,3,1,0 ("n choose 2").

% Minimum degree algorithm would prevent this fill-in easily, since it will
% pick each Vertex 2,...,8 (which have degree=1)
% to eliminate before Vertex 1 (which has degree=7 initially)


% A = [21,  1,  1,  1,  1,  1,  1,  1,  1,  1;
%     1, 12,  0,  0,  0,  0,  0,  0,  0,  1;
%     1,  0, 13,  0,  0,  0,  0,  0,  0,  1;
%     1,  0,  0, 14,  0,  0,  0,  0,  1,  0;
%     1,  0,  0,  0, 15,  0,  0,  1,  0,  0;
%     1,  0,  0,  0,  0, 16,  1,  0,  0,  0;
%     1,  0,  0,  0,  0,  1, 17,  1,  0,  1;
%     1,  0,  0,  0,  1,  0,  1, 18,  1,  0;
%     1,  0,  0,  1,  0,  0,  0,  1, 19,  0;
%     1,  1,  1,  0,  0,  0,  1,  0,  0, 20];
% Random 10x10 example matrix, where optimal<AMD<RCM (fill-in).


% check A spd
[A_L, flag] = chol(A,'lower');
if flag+~(issymmetric(A)) >0 
    error('A is not spd.') 
end


% Finding best perm
n = size(A,1);
p = perms(1:n);
k = factorial(n);

result = zeros(k,n+2); 

for i = 1:k
    q = p(i,:);
    A_perm = A(q,q); 

    [L] = chol(A_perm,'lower'); 

    fill_in = nnz(L) - nnz(tril(A_perm));

    result(i,:) = [i, fill_in, q]; 
end
% matrix of "results" 
% col 1 = index
% col 2 = fill in, 
% col 3+ = permutation of n = order of vertex elimination

sorted_result = sortrows(result, [2,1]);

best_perm = sorted_result(1,3:end);
worst_perm = sorted_result(k,3:end);

min_fill_in = sorted_result(1,2); % = best_perm_fill_in 



% super basic histogram of "number of permutations that give x amount of fill-in"
fill_in_vector = sorted_result(:,2);
histogram(fill_in_vector);
ylabel('Number of permutations');
xlabel('Fill-in produced');






% Sparsity comparison plot for: original A, optimal A, AMD A, RCM A.
% Also compares their fill-in after chol fac.

[original_L] = chol(A,'lower');
original_fill_in = nnz(original_L) - nnz(tril(A));

[best_perm_L] = chol(A(best_perm, best_perm),'lower');
best_perm_fill_in = nnz(best_perm_L) - nnz(tril(A(best_perm, best_perm)));


amd_perm = amd(sparse(A));
[amd_L] = chol(A(amd_perm,amd_perm),'lower');
amd_fill_in = nnz(amd_L) - nnz(tril(A(amd_perm,amd_perm)));


rcm_perm = symrcm(sparse(A));
[rcm_L] = chol(A(rcm_perm,rcm_perm),'lower');
rcm_fill_in = nnz(rcm_L) - nnz(tril(A(rcm_perm,rcm_perm)));




figure

subplot(2,4,1)
spy(A)
title('A')

subplot(2,4,2)
spy(A(best_perm,best_perm))
title('optimal A')

subplot(2,4,3) 
spy(A(amd_perm,amd_perm))
title('AMD A')

subplot(2,4,4)
spy(A(rcm_perm,rcm_perm))
title('RCM A')

subplot(2,4,5) 
spy(A_L)
title('chol A')
xlabel(sprintf('Fill-in: %d', original_fill_in))

subplot(2,4,6) 
spy(best_perm_L)
title('chol optimal A')
xlabel(sprintf('Fill-in: %d', best_perm_fill_in))

subplot(2,4,7) 
spy(amd_L)
title('chol AMD A')
xlabel(sprintf('Fill-in: %d', amd_fill_in))


subplot(2,4,8) 
spy(rcm_L)
title('chol RCM A')
xlabel(sprintf('Fill-in: %d', rcm_fill_in))


%todo
% add markowitz pivoting and nested dissection
% add labels for where amd, rcm, mpiv, nd. perms appear on the histogram

dissect(A)
