% Produces sparsity plots for matrices and their lower triangular matrix
% after LU decomposition.

% Used to compare how well ordering methods perform against each other, the 
% original matrix and the optimal matrix (found through brute forcing every
% permutation).

%rng(1) % set seed for reproducibility

n = 6; % square matrix dimensions
density = 1 / 3; % density of generated matrix

% generate random n x n matrix with density
A = sprand(n, n, density) + speye(n);

p = perms(1:n); % all permutations

k = factorial(n); % length of p

fill_ins = zeros(k, k); % fill-in for each permutation

% brute force
for i = 1:k
    % reorder A
    p_i = p(i, :);
    for j = 1:k
        p_j = p(j, :);

        A_perm = A(p_i, p_j);

        % fill-in
        % second input stops lu() from pivoting matrix when performing LU
        [L, U] = lu(A_perm, 0);
        fill_in = nnz(L) + nnz(U) - n - nnz(A_perm);

        fill_ins(i, j) = fill_in;
    end
end

% optimal A, choosing first permutation with min fill-in
min_ids = find(fill_ins == min(fill_ins(:)));
x = p(mod(min_ids(1), k), :); % permutation of first min's x component
y = p(ceil(min_ids(1) / k), :); % permutation of first min's y component
A_opt = A(x, y);

disp("Percentage of permutations with minimum fill-in (= " + ...
    fill_ins(min_ids(1)) + "): " + length(min_ids) / k ^ 2)

% RCM A
RCM_perm = symrcm(A);
A_RCM = A(RCM_perm, RCM_perm);

% AMD A
AMD_perm = amd(A);
A_AMD = A(AMD_perm, AMD_perm);

% Markowitz A
[P, Q] = MarkowitzPivoting(full(A));
A_Markowitz = A(P, Q);

% nested dissection A
dissection_perm = dissect(A);
A_dissection = A(dissection_perm, dissection_perm);

% cell array of matrices
matrices = {A, A_opt, A_RCM, A_AMD, A_Markowitz, A_dissection};
matrices_names = ["Original", "Optimal", "RCM", "AMD", ...
    "Markowitz", "Dissection"];
m = length(matrices);

% histogram for number of permutations that give each fill-in value
figure
histogram(fill_ins);
title("Permutation count for each fill-in value")
ylabel("Number of permutations");
xlabel("Fill-in");

% sparsity plots
figure
for i = 1:m
    matrix_i = matrices{i};

    % fill-in
    % second input stops lu() from pivoting matrix when performing LU
    [L, U] = lu(matrix_i, 0);
    matrix_i_fill = nnz(L) + nnz(U) - n - nnz(matrix_i);

    % plotting
    % sparsity pattern
    subplot(2, m, i)
    spy(matrix_i)
    title(matrices_names(i) + " A")

    % L sparsity pattern
    subplot(2, m, i + m)
    spy(L)
    title("Lower " + matrices_names(i) + " A")
    xlabel("Fill-in: " + matrix_i_fill)
end
