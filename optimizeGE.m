function [Aopt,rperm,cperm] = optimizeGE(A)

% Attempt to reorder the rows and columns of matrix A
% so as to minimize the cost of a subsequent Gaussian
% Elimination

% The entries of A are
% 0 if guaranteed empty;
% 1 if guaranteed non-zero (possible pivot)
% 0.5 if not known - need to eliminate but not guaranteed pivot

% Aopt is the reordered matrix and rperm and cperm index
% the reordered rows and columns

% Markersize for spy
ms = 32;

[m,n] = size(A);

% Initialize permutation arrays
rperm = [1:n];
cperm = [1:n];

% Make a working copy
B = A(rperm,cperm);

% And check the GE cost for this matrix
spy(B,ms);
set(gca,'fontsize',16)
cost = count_cost(B)

% Loop over pivot points
for  p = 1:n

    % Assume p-1 rows and columns are already done
    % Loop over what remains to find the next best pivot
    bestcost = 1e8;
    bestr = 0;
    bestc = 0;
    % Loop over rows and columns of the remaining lower right square matrix
    for r = p:n
        for c = p:n
            % Must be 1 to be candidate pivot
            if B(r,c) == 1
                % Work out the cost and benefit of eliminating column c by pivoting on (r,c)
                nnz_rowr = sum(sign(B(r,p:n)));
                n2elim = sum(sign(B(p:n,c)));
                pivot_cost = nnz_rowr*n2elim;
                benefit = n2elim - 1;
                fill_in = 0;
                for r1 = p:n
                    if r1 ~= r
                        % Check if we're eliminating
                        if B(r1,c) > 0
                            % Find fill-in
                            nfill_rowr = sum(B(r,:) > 0 & (B(r1,:) == 0));
                            fill_in = fill_in + nfill_rowr;
                        end
                    end
                end
                % Net cost
                cost = pivot_cost - benefit + fill_in;
                disp([' Pivot ',num2str(r),' ',num2str(c),...
                      ' cost ',num2str(pivot_cost),'-',num2str(benefit),'+',num2str(fill_in),...
                      '=',num2str(cost)])
                % Is it best?
                if cost < bestcost
                    bestcost = cost;
                    bestr = r;
                    bestc = c;
                end
            end
        end
    end
    % Implement the best pivot
    disp(['p = ',num2str(p),' Best pivot is ',num2str(bestr),' ',num2str(bestc),' cost ',num2str(bestcost)])
    if bestr == 0
        disp('No pivot found')
        B
        pause
    end
    % Update the permutation vectors rperm and cperm
    temp = rperm(p);
    rperm(p) = rperm(bestr);
    rperm(bestr) = temp;
    temp = cperm(p);
    cperm(p) = cperm(bestc);
    cperm(bestc) = temp;
    % Update the working copy matrix B 
    B = A(rperm,cperm);

    % And check the full GE cost for this matrix
    spy(B,ms);
    set(gca,'fontsize',16)
    title(['Done pivot ',num2str(bestr),' ',num2str(bestc)])
    cost = count_cost(B)
    pause

    % Now implement the eliminations of the first p pivots
    % so that we know where we stand, ready to find next pivot
    for c = 1:p
        % B(c,c) should now be a guaranteed pivot
        % Which rows need to be eliminated?
        for r = c+1:n
            if B(r,c) == 1
                % Need to eliminate in this row
                B(r,c) = 0;
                B(r,c+1:n) = max(B(r,c+1:n), B(c,c+1:n));
            elseif B(r,c) == 0.5
                % Might need to eliminate in this row but can't be sure
                B(r,c) = 0;
                B(r,c+1:n) = max(B(r,c+1:n), 0.5);
            end
        end
    end

end

% And save final matrix
Aopt = A(rperm,cperm);