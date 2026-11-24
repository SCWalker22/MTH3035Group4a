include_original = true;

n_iterations = 100; %Number of random matrices to generate & optimise
n_bins = 20; % Number of bins for histogram

% Set desired size of matrices
n_rows = 100;
n_cols = n_rows;

methods = ["Original", "Approximate Minimum Degree", "Reverse Cuthill-Mckee", "Markowitz Pivoting"];

% Arrays to store results
fill_ins = zeros(n_iterations, 3 + include_original);

% Loop through n iterations
for i = 1:n_iterations
    S = sprand(n_rows, n_cols, 0.1) + eye(n_rows);
    A_temp = full(S);
    A = transpose(A_temp)*A_temp; % Ensure A is SPD by taking A^T A

    [original_L] = chol(A,'lower');
    original_fill_in = nnz(original_L) - nnz(tril(A));

    % Compute permutation using Approx Min Degree
    amd_perm = amd(sparse(A));
    [amd_L] = chol(A(amd_perm,amd_perm),'lower');
    amd_fill_in = nnz(amd_L) - nnz(tril(A(amd_perm,amd_perm)));

    %Compute permutation using Reverse Cuthill-McKee
    rcm_perm = symrcm(sparse(A));
    [rcm_L] = chol(A(rcm_perm,rcm_perm),'lower');
    rcm_fill_in = nnz(rcm_L) - nnz(tril(A(rcm_perm,rcm_perm)));

    % Compute permutation using Markowitz Pivoting
    [pivot_row_perm, pivot_col_perm] = MarkowitzPivotingOld(A);
    pivot_L = chol(A(pivot_row_perm, pivot_col_perm), 'lower');
    pivot_fill_in = nnz(pivot_L) - nnz(tril(A(pivot_row_perm, pivot_col_perm)));


    % Add num non-zeros to arrays
    if include_original
        fill_ins(i, 1) = original_fill_in;
    end
    fill_ins(i, 1+include_original) = amd_fill_in;
    fill_ins(i, 2+include_original) = rcm_fill_in;
    fill_ins(i, 3+include_original) = pivot_fill_in;
end

% Find min and max of x, so all graphs line up
x_left = min(min(fill_ins));
x_right = max(max(fill_ins));

% Start graphing
figure()

[m, n] = size(fill_ins);

for i = 1:n
    subplot(n, 1, i)
    histogram(fill_ins(:, i), n_bins)
    title(methods(i + 1 - include_original))
    xlim([x_left, x_right])
end