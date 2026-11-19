% Calculates total fill-in of matrix

% plots sparsity of matrix and lower matrix
function plot_sparsity(matrix, name, m, index)
    [L, U, P] = lu(matrix);

    % fill-in
    matrix_fill = nnz(L) + nnz(U) - nnz(P * matrix);

    % plotting
    % sparsity pattern
    subplot(2, m, index)
    spy(matrix)
    title(name + " A")

    % lower sparsity pattern
    subplot(2, m, index + m)
    spy(L)
    title("Lower " + name + " A")
    xlabel("Total fill-in: " + matrix_fill)
end

%rng(1) % set seed

n = 10; % square matrix dimensions
density = 1 / 3; % density of generated matrix

% generate random n x n matrix with density
A = full(spones(sprand(n, n, density)));

p = perms(1:n); % all permutations

k = factorial(n); % length of p

fill_ins = zeros(k, 1); % fill-in for each permutation

% brute force
for i = 1:k
    % reorder A
    perm = p(i, :);
    A_perm = A(perm, perm);

    % fill-in
    [L, U, P] = lu(A_perm);
    fill_in = nnz(L) + nnz(U) - nnz(P * A_perm);

    fill_ins(i) = fill_in;
end

% optimal A
[~, I] = min(fill_ins); % row index with minimum fill-in
opt_perm = p(I, :);
A_opt = A(opt_perm, opt_perm);

% RCM A
RCM_perm = symrcm(A);
A_RCM = A(RCM_perm, RCM_perm);

% AMD A
AMD_perm = amd(A);
A_AMD = A(AMD_perm, AMD_perm);

% Markowitz A
Markowitz_perm = MarkowitzPivoting(A);
A_Markowitz = A(Markowitz_perm, Markowitz_perm);

% nested dissection A
dissection_perm = dissect(A);
A_dissection = A(dissection_perm, dissection_perm);

% cell array of matrices
matrices = {A, A_opt, A_RCM, A_AMD, A_Markowitz, A_dissection};
m = length(matrices);

% names of matrices
matrices_names = ["Original", "Optimal", "RCM", "AMD", "Markowitz", "Dissection"];

% histogram for number of permutations that give each fill-in value
figure
histogram(fill_ins);
title("Permutation count for each fill-in value")
ylabel("Number of permutations");
xlabel("Fill-in");

% sparsity plots
figure
for i = 1:m
    plot_sparsity(matrices{i}, matrices_names(i), m, i)
end
