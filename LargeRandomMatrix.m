n_iterations = 100;
n_bins = 20;
amd_fill_ins = zeros(n_iterations, 1);
rcm_fill_ins = zeros(n_iterations, 1);
pivot_fill_ins = zeros(n_iterations, 1);
partial_pivot_fill_ins = zeros(n_iterations, 1);
for i = 1:n_iterations
    n_rows = 100;
    n_cols = n_rows;
    S = sprand(n_rows, n_cols, 0.1) + eye(n_rows);
    A_temp = full(S);
    A = transpose(A_temp)*A_temp;

    % check A spd
    [A_L, flag] = chol(A,'lower');
    if flag+~(issymmetric(A)) >0 
        A = A + transpose(A); % Make A SPD by adding transpose
        % error('A is not spd.') 
    end

    amd_perm = amd(sparse(A));
    [amd_L] = chol(A(amd_perm,amd_perm),'lower');
    amd_fill_in = nnz(amd_L) - nnz(tril(A(amd_perm,amd_perm)));


    rcm_perm = symrcm(sparse(A));
    [rcm_L] = chol(A(rcm_perm,rcm_perm),'lower');
    rcm_fill_in = nnz(rcm_L) - nnz(tril(A(rcm_perm,rcm_perm)));

    [pivot_row_perm, pivot_col_perm] = MarkowitzPivotingOld(A);
    pivot_L = chol(A(pivot_row_perm, pivot_col_perm), 'lower');
    pivot_fill_in = nnz(pivot_L) - nnz(tril(A(pivot_row_perm, pivot_col_perm)));

    [partial_pivot_row_perm, partial_pivot_col_perm] = MarkowitzPivoting(A);
    partial_pivot_L = chol(A(partial_pivot_row_perm, partial_pivot_col_perm), 'lower');
    partial_pivot_fill_in = nnz(pivot_L) - nnz(tril(A(partial_pivot_row_perm, partial_pivot_col_perm)));
    
    amd_fill_ins(i) = amd_fill_in;
    rcm_fill_ins(i) = rcm_fill_in;
    pivot_fill_ins(i) = pivot_fill_in;
    partial_pivot_fill_ins(i) = partial_pivot_fill_in;
end

x_left = min(min(min(amd_fill_ins), min(rcm_fill_ins)), min(pivot_fill_ins));
x_right = max(max(max(amd_fill_ins), max(rcm_fill_ins)), max(pivot_fill_ins));

figure()

subplot(3, 1, 1)
histogram(amd_fill_ins, n_bins)
title('Approximate Minimum Degree')
xlim([x_left, x_right])

subplot(3, 1, 2)
histogram(rcm_fill_ins, n_bins)
title('Reverse Cuthill-McKee')
xlim([x_left, x_right])

subplot(3, 1, 3)
histogram(pivot_fill_ins, n_bins)
title('Markowitz Pivoting')
xlim([x_left, x_right])

% subplot(4, 1, 4)
% histogram(partial_pivot_fill_ins, n_bins)
% title('Partial Pivoting')
% xlim([x_left, x_right])