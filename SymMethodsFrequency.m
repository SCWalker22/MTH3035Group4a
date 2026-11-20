% Produces a heatmap for frequencies of entries for each method (if method
% is ever optimal) when the method produces the lowest fill-in.
% Can be used to see if there's a pattern for when each method is optimal
% Use a tolerance to ignore entries with frequency lower than
% tolerance (tolerance is a percentage) for visualisation.

% Notes:
%   - unsure on variable values
%       - increased density, not sure what's best for a large sparse matrix
%       - tolerance (use higher tol for amd - almost always produces lowest 
%         fill-in)
%   - changing variables might give a better pattern for each method??
%       - can't see a pattern when there's one extremely dominant method
%           - changing the methods comparison list to remove the dominant
%             method(s) can give a good comparison

n = 100; % square matrix dimensions
density =  5 / n; % density of generated matrix
t = 1000; % test number

% cell array of ordering method functions
methods = {@symrcm, @symamd, @MarkowitzPivoting, @dissect}; 
method_names = ["RCM", "AMD", "Markowitz", "Dissection"];
m = length(methods);

% frequency of matrix for each method when method produces lowest fill-in
method_A_freq = repmat({zeros(n, n)}, m, 1);

fill_in = zeros(m, 1); % fill-in for each method at each step
count = zeros(m, 1); % count for times each method produces lowest fill-in

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
    min_ids = find(fill_in == min(fill_in))';

    % adding A to optimal frequency matrices
    for j = min_ids
        method_A_freq{j} = method_A_freq{j} + A;
        count(j) = count(j) + 1;
    end
end

% plotting shapes of overall matrices as heatmap
for i = 1:m
    c = count(i);

    % ignore empty matrices
    if c ~= 0
        % tolerance
        tol = 0.3;

        A_freq = method_A_freq{i};

        % diagonal values are much larger than any other point (due to
        % sparsity), set diagonals to maximum non-diagonal entry
        max_value = max(A_freq(A_freq < c));
        A_freq(c <= A_freq) = max_value;

        % ignore points lower than tolerance
        A_freq(A_freq < max_value * tol) = 0;
        
        % final matrix's heatmap plot
        figure
        heatmap(A_freq, "GridVisible", "off", "Colormap", parula, ...
            "CellLabelColor", "none", ...
            "ColorLimits", [floor(max_value * tol) max_value])
        title(method_names(i) + " frequency plot")
    end

    disp(method_names(i) + " optimal method count: " + c)
end
