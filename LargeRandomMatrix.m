include_original = true;

n_iterations = 100; %Number of random matrices to generate & optimise
n_bins = 20; % Number of bins for histogram

% Set desired size of matrices
n_rows = 100;
n_cols = n_rows;

% Arrays to store results
original_fill_ins = zeros(n_iterations, 1);
amd_fill_ins = zeros(n_iterations, 1);
rcm_fill_ins = zeros(n_iterations, 1);
pivot_fill_ins = zeros(n_iterations, 1);
partial_pivot_fill_ins = zeros(n_iterations, 1);

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

    % Compute permutation using Partial Pivoting
    [partial_pivot_row_perm, partial_pivot_col_perm] = MarkowitzPivoting(A);
    partial_pivot_L = chol(A(partial_pivot_row_perm, partial_pivot_col_perm), 'lower');
    partial_pivot_fill_in = nnz(pivot_L) - nnz(tril(A(partial_pivot_row_perm, partial_pivot_col_perm)));
    
    % Add num non-zeros to arrays
    original_fill_ins(i) = original_fill_in;
    amd_fill_ins(i) = amd_fill_in;
    rcm_fill_ins(i) = rcm_fill_in;
    pivot_fill_ins(i) = pivot_fill_in;
    partial_pivot_fill_ins(i) = partial_pivot_fill_in;
end

% Find min and max of x, so all graphs line up
x_left = min(min(min(amd_fill_ins), min(rcm_fill_ins)), min(pivot_fill_ins));
x_right = max(max(max(amd_fill_ins), max(rcm_fill_ins)), max(pivot_fill_ins));

if include_original
    x_left = min(min(original_fill_ins), x_left);
    x_right = max(max(original_fill_ins), x_right);
end

figure()
if include_original
    subplot(4, 1, 1)
    histogram(original_fill_ins, n_bins)
    title('Original Fill-in')
    xlim([x_left, x_right])
end

subplot(3 + include_original, 1, 1 + include_original)
histogram(amd_fill_ins, n_bins)
title('Approximate Minimum Degree')
xlim([x_left, x_right])

subplot(3 + include_original, 1, 2 + include_original)
histogram(rcm_fill_ins, n_bins)
title('Reverse Cuthill-McKee')
xlim([x_left, x_right])

subplot(3 + include_original, 1, 3 + include_original)
histogram(pivot_fill_ins, n_bins)
title('Markowitz Pivoting')
xlim([x_left, x_right])

% subplot(4, 1, 4)
% histogram(partial_pivot_fill_ins, n_bins)
% title('Partial Pivoting')
% xlim([x_left, x_right])