function [P, Q] = MarkowitzPivoting(A)
    % Markowitz tolerance (range: [1e-4, 0.999], set to 0 to ignore)    
    u = 0.0078125;

    [m, n] = size(A); % Matrix dimensions
    
    % Duplicate of array A
    B = A;

    best_row = -1;
    best_col = -1;
    
    % Initial permutations for rows and cols
    best_rperm = 1:m;
    best_cperm = 1:n;

    % Loop over min(m, n) pivots
    for p = 1:min(m, n)
        % Set initial cost infinite, so we always find a pivot
        best_markowitz_cost = Inf;
        
        % Lowest cost is 0, if a 0 cost pivot is found, break from loops
        zero_cost = false;

        % Rearrange A
        A = B(best_rperm, best_cperm); 
        
        % Loop through remaining rows
        for r = p:m
            % Loop through remaining cols
            for c = p:n
                val = A(r, c); % value at current row and column
                % Ignore non-zero pivots 
                if val ~= 0
                    col = A(p:m, c); % current column

                    % Markowitz tolerance
                    if abs(val) >= u * max(abs(col))
                        nnz_r = nnz(A(r, p:n)); % Non-zeros in row
                        nnz_c = nnz(col); % Non-zeros in col
                        
                        % Calculate Markowitz cost
                        markowitz_cost = (nnz_r - 1) * (nnz_c - 1);

                        % Update best pivot
                        if markowitz_cost < best_markowitz_cost
                            best_row = r;
                            best_col = c;

                            % Break if zero cost is found
                            if markowitz_cost == 0
                                zero_cost = true;
                                break
                            end
                            
                            best_markowitz_cost = markowitz_cost;
                        end
                    end
                end
            end

            % Break if zero cost is found
            if zero_cost
                break
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

    % Return optimal permutations
    P = best_rperm;
    Q = best_cperm;
end
