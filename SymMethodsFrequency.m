% Produces a heatmap for frequencies of entries for each method (if method
% is ever 'optimal') when the method produces a low enough fill-in.

% Used to see if there's a pattern for when each method is 'optimal' and
% which methods consistently produce 'low' fill-ins.

% Can't find the actual optimal permutation and fill-in using brute force 
% for high n and high tests, instead use a fill-in threshold for methods - 
% when a method produces fill-in lower than the threshold, consider it to 
% have produced an 'optimal' permutation.

% Use plot tolerance 'p_tol' to ignore the 'n x n - (2 * p_tol + n)' lowest 
% entries (filters highest values in frequency matrix to see shape of 
% 'optimal' matrices).
%   - e.g. if you want to plot the 100 most frequent non-diagonal values, 
%     set p_tol = 50.
% Note that the actual number of points plotted will have a minimum of 2 *
% p_tol but can have a high ceiling - if the lowest frequency number 
% that passes the threshold is very common in the frequency matrix.
%   - a good f_tol should produce a low-medium optimal method count a 
%     method to see an accurate frequency plot.

n = 100; % square matrix dimensions
density =  2 / n; % density of generated matrix
t = 2500; % test number

% cell array of ordering method functions
methods = {@symrcm, @symamd, @MarkowitzPivoting, @dissect}; 
method_names = ["RCM", "AMD", "Markowitz", "Dissection"];
m = length(methods);

% if a permutation is found with fill-in <= f_tol, found an 'optimal' order
% use 'SymMethodsComparison.m' to get a good f_tol value
f_tol = 50;

% frequency of matrix for each method when method produces 'optimal' perm
method_A_freq = repmat({zeros(n, n)}, m, 1);

count = zeros(m, 1); % count for times each method produces 'optimal' perm

% plot the '2 * p_tol' highest non-diagonal values of each frequency matrix
p_tol = 100;

% testing
for i = 1:t
    %rng(i); % set seed for reproducibility

    % generate random n x n spd matrix with density
    A = sprandsym(n, density) + n * eye(n);

    % fill-in for each method
    for j = 1:m
        % reorder A
        perm = methods{j}(A);
        A_perm = A(perm, perm);

        % fill-in
        A_perm_chol = chol(A_perm, "lower");
        fill_in = nnz(A_perm_chol) - nnz(tril(A_perm));

        % adding A to optimal frequency matrices
        if fill_in <= f_tol
            method_A_freq{j} = method_A_freq{j} + spones(A);
            count(j) = count(j) + 1;
        end
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
        A_freq(A_freq >= c) = max_value;

        % lowest values that pass threshold
        min_mv = min(maxk(A_freq(:), p_tol + n));

        % ignore points that don't meet the threshold
        A_freq(A_freq < min_mv) = 0;

        % final matrix's heatmap plot
        figure
        heatmap(A_freq, "GridVisible", "off", "Colormap", copper, ...
            "CellLabelColor", "none", ...
            "ColorLimits", [max(min_mv - 1, 0) max_value])
        title(method_names(i) + " frequency plot")
    end

    disp(method_names(i) + " optimal method count: " + c)
end
