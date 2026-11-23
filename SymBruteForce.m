% Produces sparsity plots for matrices and their lower triangular matrix
% after Cholesky decomposition.

% Used to compare how well ordering methods perform against each other, the 
% original matrix and the optimal matrix (found through brute forcing every
% permutation).

%rng(1) % set seed for reproducibility

n = 10; % square matrix dimensions
density = 1 / 3; % density of generated matrix

% generate random n x n spd matrix with density
A = sprandsym(n, density) + n * eye(n);

p = perms(1:n); % all permutations

k = factorial(n); % length of p

fill_ins = zeros(k, 1); % fill-in for each permutation

% brute force
for i = 1:k
    % reorder A
    perm = p(i, :);
    A_perm = A(perm, perm);

    % fill-in
    A_chol = chol(A_perm, "lower");
    fill_in = nnz(A_chol) - nnz(tril(A_perm));

    fill_ins(i) = fill_in;
end

% optimal A
min_ids = find(fill_ins == min(fill_ins)); % row indexes with min fill-in
opt_perm = p(min_ids(1), :); % choosing first permutation with min fill-in
A_opt = A(opt_perm, opt_perm);

disp("Percentage of permutations with minimum fill-in (= " + ...
    fill_ins(min_ids(1)) + "): " + length(min_ids) / k)

% RCM A
RCM_perm = symrcm(A);
A_RCM = A(RCM_perm, RCM_perm);

% AMD A
AMD_perm = symamd(A);
A_AMD = A(AMD_perm, AMD_perm);

% Markowitz A
Markowitz_perm = MarkowitzPivoting(A);
A_Markowitz = A(Markowitz_perm, Markowitz_perm);

% nested dissection A
dissection_perm = dissect(A);
A_dissection = A(dissection_perm, dissection_perm);

% cell array of matrices
matrices = {A, A_opt, A_RCM, A_AMD, A_Markowitz, A_dissection};
matrices_names = ["Original", "Optimal", "RCM", "AMD", "Markowitz", "Dissection"];
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
    matrix_i_chol = chol(matrix_i, "lower");
    matrix_i_fill = nnz(matrix_i_chol) - nnz(tril(matrix_i)); 

    % plotting
    % sparsity pattern
    subplot(2, m, i)
    spy(matrix_i)
    title(matrices_names(i) + " A")

    % chol sparsity pattern
    subplot(2, m, i + m)
    spy(matrix_i_chol)
    title("Chol " + matrices_names(i) + " A")
    xlabel("Fill-in: " + matrix_i_fill)
end
