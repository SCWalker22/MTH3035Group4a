% Produces a heatmap for frequencies of entries for each method (if method
% is ever 'optimal') when the method produces a low enough fill-in.

% Used to see if there's a pattern for when each method is 'optimal' and
% which methods consistently produce 'low' fill-ins.

% Compares methods to each other to find the 'optimal' method. If one
% method consistently produces the 'optimal' ordering, can't see patterns
% very well. Remove dominant method from list if needed.

% Use plot tolerance 'p_tol' to ignore the 'n x n - (p_tol + n)' lowest 
% entries (filter highest values in frequency matrix to see shape of 
% 'optimal' matrices).
%   e.g. if you want to plot the 100 most frequent non-diagonal values, set 
%        p_tol = 100.

n = 100; % square matrix dimensions
density =  5 / n; % density of generated matrix
t = 1000; % test number

% cell array of ordering method functions
methods = {@symrcm, @symamd, @MarkowitzPivoting, @dissect}; 
method_names = ["RCM", "AMD", "Markowitz", "Dissection"];
m = length(methods);

% frequency of matrix for each method when method produces 'optimal' perm
method_A_freq = repmat({zeros(n, n)}, m, 1);

fill_in = zeros(m, 1); % fill-in for each method at each step
count = zeros(m, 1); % count for times each method produces 'optimal' perm

% plot the 'p_tol' highest non-diagonal values of each frequency matrix
p_tol = 250;

% testing
for i = 1:t
    % generate random n x n spd matrix with density
    A = spones(sprandsym(n, density)) + n * eye(n);

    % fill-in for each method
    for j = 1:m
        % reorder A
        perm = methods{j}(A);
        A_perm = A(perm, perm);

        % fill-in
        A_perm_chol = chol(A_perm, "lower");
        fill_in(j) = nnz(A_perm_chol) - nnz(tril(A_perm));
    end

    % minimum fill-in indexes
    min_ids = find(fill_in == min(fill_in));

    % adding A to optimal frequency matrices
    for j = min_ids
        method_A_freq{j} = method_A_freq{j} + A;
        count(j) = count(j) + 1;
    end
end

% plotting shapes of overall matrices as heatmap
for i = 1:m
    c = count(i);

    % ignore empty frequency matrices
    if c ~= 0
        A_freq = method_A_freq{i};

        % diagonal values are much larger than any other point (due to
        % sparsity), set diagonals to maximum non-diagonal entry
        max_value = max(1, max(A_freq(A_freq < c)));
        A_freq(c <= A_freq) = max_value;

        max_values = maxk(A_freq(:), p_tol + n); % accounting for diagonals
        min_mv = min(max_values); % lowest values that pass threshold

        % ignore 'n x n - (tol + n)' lower value points
        A_freq(A_freq < min(max_values)) = 0;

        % final matrix's heatmap plot
        figure
        heatmap(A_freq, "GridVisible", "off", "Colormap", copper, ...
            "CellLabelColor", "none", ...
            "ColorLimits", [min_mv max_value])
        title(method_names(i) + " frequency plot")
    end

    disp(method_names(i) + " optimal method count: " + c)
end
