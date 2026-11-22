n_rows = 100;
n_cols = n_rows;
S = sprand(n_rows, n_cols, 0.1) + eye(n_rows); %Add to trasnpose
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

% figure

% subplot(2,4,1)
% spy(A)
% title('A')

% subplot(2,4,2)
% spy(A(pivot_row_perm,pivot_col_perm))
% title('Markowitz A')

% subplot(2,4,3) 
% spy(A(amd_perm,amd_perm))
% title('AMD A')

% subplot(2,4,4)
% spy(A(rcm_perm, rcm_perm))
% title('RCM A')

% subplot(2,4,5) 
% spy(A_L)
% title('chol A')
% xlabel(sprintf('Fill-in: %d', original_fill_in))

% subplot(2,4,6) 
% spy(pivot_L)
% title('chol Pivot A')
% xlabel(sprintf('Fill-in: %d', pivot_fill_in))

% subplot(2,4,7) 
% spy(amd_L)
% title('chol AMD A')
% xlabel(sprintf('Fill-in: %d', amd_fill_in))


% subplot(2,4,8) 
% spy(rcm_L)
% title('chol RCM A')
% xlabel(sprintf('Fill-in: %d', rcm_fill_in))

% %subplot(2, 4, 9)
% %spy(pivot_L)
% %title('chol Pivot A')
% %xlabel(sprintf('Fill-in: %d', pivot_fill_in))

% %todo
% % add nested dissection
% % add labels for where amd, rcm, mpiv, nd. perms appear on the histogram

% dissect(A)
