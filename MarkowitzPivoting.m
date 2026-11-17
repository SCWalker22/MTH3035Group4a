function [P, Q] = MarkowitzPivoting(A)
    [m, n] = size(A); % Find the size of the matrix
    rperm = 1:n;
    cperm = 1:n;

    %Create duplicate of array A
    B = A(rperm, cperm);
    best_row = 1;
    best_col = 1;

    %Create initial permutations for rows and cols
    best_rperm = 1:n;
    best_cperm = 1:n;
    
    for p = 1:n %Loop over pivots - we will have to use n pivots
        best_markowitz_cost = inf; %Set initial cost infinite, so we always find a pivot
        A = B(best_rperm, best_cperm); %Rearrange A as defined so far
        for r = p:n %Loop through Rows
            for c = p:n % Loop through cols
                if A(r, c) ~= 0 % If this entry is non-zero
                    c_i = sum(sign(A(r, p:n))); %Num zeros in col
                    r_i = sum(sign(A(p:n, c))); % Num zeros in row
                    markowitz_cost = (c_i - 1)*(r_i - 1); %Calculate Markowitz cost
                    if markowitz_cost < best_markowitz_cost %Update best pivot if needed
                        best_row = r;
                        best_col = c;
                        best_markowitz_cost = markowitz_cost;
                    end
                end
            end
        end
        % Update the row, col permutations to the optimum
        previous_r = best_rperm(p);
        previous_c = best_cperm(p);
        best_rperm(p) = best_rperm(best_row);
        best_cperm(p) = best_cperm(best_col);
        best_rperm(best_row) = previous_r;
        best_cperm(best_col) = previous_c;
    end
    P = best_rperm; % Return optimal permutations
    Q = best_cperm;
end


%A = [11, 0, 0, 0, 0, 1, 1;
 %   0, 12, 1, 1, 1, 0, 1;
  %  0, 1, 13, 0, 0, 0, 1;
   % 0, 1, 0, 14, 0, 0, 1;
    %0, 1, 0, 0, 15, 0, 0;
    %1, 0, 0, 0, 0, 16, 0;
    %1, 1, 1, 1, 0, 0, 17];

%[P, Q] = MarkowitzPivoting(A)


