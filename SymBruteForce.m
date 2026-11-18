% Notes: 
%   - different version of brute force, not necessarily better than 
%     other, used code from other
%   - histogram not added
%   - A's of interest not added

% plots sparsity of matrix and cholesky of matrix
function plot_sparsity(matrix, name, m, index)
    matrix_chol = chol(matrix, "lower");
    
    % fill-in
    matrix_fill = nnz(matrix_chol) - nnz(tril(matrix)); 

    % plotting
    % sparsity pattern
    subplot(2, m, index)
    spy(matrix)
    title(name + " A")

    % chol sparsity pattern
    subplot(2, m, index + m)
    spy(matrix_chol)
    title("Chol " + name + " A")
    xlabel("Fill-in: " + matrix_fill)
end

%rng(1) % set seed

n = 10; % square matrix dimensions
density = 1 / 3; % density of generated matrix

% generate random n x n spd matrix with density
A = spones(sprandsym(n, density)) + n * eye(n);

p = perms(1:n); % all permutations

k = factorial(n); % length of p

% storing results, rows: [fill-in, permutation]
results = zeros(k, n + 1);

% brute force
for i = 1:k
    % reorder A
    q = p(i, :);
    A_perm = A(q, q);

    % fill-in
    L = chol(A_perm, "lower");
    fill_in = nnz(L) - nnz(tril(A_perm));

    results(i, :) = [fill_in, q];
end

% optimal A
[~, I] = min(results(:, 1)); % row with minimum fill-in
min_fill = results(I, :);
opt_perm = min_fill(2:end);
A_opt = A(opt_perm, opt_perm);

% RCM A
RCM_perm = symrcm(A);
A_RCM = A(RCM_perm, RCM_perm);

% AMD A
AMD_perm = amd(A);
A_AMD = A(AMD_perm, AMD_perm);

% Markowitz A
Markowitz_perm = MarkowitzPivoting(A);
Markowitz_A = A(Markowitz_perm, Markowitz_perm);

% cell array of matrices
matrices = {A, A_opt, A_RCM, A_AMD, Markowitz_A};
m = length(matrices);

% names of matrices
matrices_names = ["Original", "Optimal", "RCM", "AMD", "Markowitz"];

figure
for i = 1:m
    plot_sparsity(matrices{i}, matrices_names(i), m, i)
end
