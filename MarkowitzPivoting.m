function [P, Q] = markowitz_pivot(A)
    [m, n] = size(A);
    rperm = 1:n;
    cperm = 1:n;

    %Create duplicate of array A
    B = A(rperm, cperm);
    best_row = 1;
    best_col = 1;

    best_rperm = zeros(1,n);
    best_cperm = zeros(1,n);
    
    for p = 1:n %Loop over pivots
        best_markowitz_cost = inf;
        for r = p:n
            for c = p:n
                if B(r, c) == 1
                    c_i = sum(sign(B(r, p:n)));
                    r_i = sum(sign(B(p:n, c)));
                    markowitz_cost = (c_i - 1)*(r_i - 1);
                    if markowitz_cost < best_markowitz_cost
                        best_row = r;
                        best_col = c;
                        best_markowitz_cost = markowitz_cost;
                    end
                end
            end
        end
        best_rperm(p) = best_row;
        best_cperm(p) = best_col;
    end
    P = best_rperm;
    Q = best_cperm;
end


A = [11, 0, 0, 0, 0, 1, 1;
    0, 12, 1, 1, 1, 0, 1;
    0, 1, 13, 0, 0, 0, 1;
    0, 1, 0, 14, 0, 0, 1;
    0, 1, 0, 0, 15, 0, 0;
    1, 0, 0, 0, 0, 16, 0;
    1, 1, 1, 1, 0, 0, 17];

markowitz_pivot(A)


