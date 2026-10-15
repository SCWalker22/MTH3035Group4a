function [cost] = count_cost(A)
% Count the cost of doing Gaussian elimination on the sparse
% matrix A if the sparsity pattern is known in advance
% Brute force calculation stepping through the elimination process
% and counting floating point operations.

% The entries of A are
% 0 if guaranteed empty;
% 1 if guaranteed non-zero (possible pivot)
% 0.5 if not known - need to eliminate but not guaranteed pivot

[m,n] = size(A);

% Asave = A;
% 
% if m ~= n
%     disp('A must be square')
%     pause
% end
% 

% Flag to indicate invalid cost;
% sg = -1 (=> negative cost) means diagonal element was not
% a guaranteed pivot at some point.
sg = 1;

% Keep track of fill in
nfill = 0;

% Downward sweep
costdown = 0;
for c = 1:n
    if A(c,c) < 1
        disp(['A(',num2str(c),',',num2str(c),') is not a guaranteed pivot'])
        sg = -1;
        pause
    end
    % Cost of rescaling row c = number of nonzeros in this row
    nnz_rowc = sum(sign(A(c,:)));
    costdown = costdown + nnz_rowc;
    % Which rows need to be eliminated in column c? Loop over rows.
    for r = c+1:n
        if A(r,c) > 0
            % We need to eliminate in this row
            costdown = costdown + nnz_rowc;
            % Count any fill-in
            nfill_rowr = sum(A(c,:) > 0 & (A(r,:) == 0));
            nfill = nfill + nfill_rowr;
            % Eliminate A(r,c)
            A(r,c) = 0;
            % and update this row with fill in where needed
            A(r,c+1:n) = max(A(r,c+1:n), A(c,c+1:n));
%             Uncomment these lines to visualise the elimination
%             step-by-step, for checking
%             subplot(1,2,1)
%             spy(A)
%             The second plot will differ from the first if any entries are
%             0.5 rather than 1
%             subplot(1,2,2)
%             spy(A==1)
%             pause
        end
    end
end

% Upward sweep
% The upward sweep produces no fill-in, so just count the number of
% nonzero entries, minus the diagonal.
costup = sum(sign(A),'all') - n;

disp('costdown  costup  fill-in')
[costdown,costup,nfill]

% Total cost
cost = costdown + costup;
% Flag invalid result
cost = sg*cost;

end