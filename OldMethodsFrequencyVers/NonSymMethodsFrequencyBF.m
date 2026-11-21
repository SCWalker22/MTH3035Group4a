% Produces a heatmap for frequencies of entries for each method (if method
% is ever 'optimal') when the method produces 'lowest' fill-in.

% Used to see if there's a pattern for when each method is 'optimal' and
% which methods consistently produce 'low' fill-ins.

% Can't find the optimal permutation and fill-in using brute force, instead
% pick 'p_n' random permutations and find the minimum fill-in for current 
% A. Use this minimum fill-in as a benchmark for ordering methods:
%   - if a method produces a fill-in above the BF fill-in, not optimal
%   - if a method produces a fill-in equal to the BF fill-in, optimal
%   - if any method produces a fill-in less than the BF fill-in, BF optimal 
%     permutation isn't valid, move to next test

% For the non-symmetric case, finding an optimal ordering using BF with a
% fill-in less than the ordering methods requires a high 'p_n'. Using this
% method to find the patterns of optimal orderings isn't great.

% Use plot tolerance 'p_tol' to ignore the 'n x n - (p_tol + n)' lowest 
% entries (filter highest values in frequency matrix to see shape of 
% 'optimal' matrices).
%   e.g. if you want to plot the 100 most frequent non-diagonal values, set 
%        p_tol = 100.

n = 50; % square matrix dimensions
density = 1 / n; % density of generated matrix
t = 500; % test number

% cell array of ordering method functions
methods = {@symrcm, @symamd, @MarkowitzPivoting, @dissect}; 
method_names = ["RCM", "AMD", "Markowitz", "Dissection"];
m = length(methods);

% if a permutation in brute force is found <= tol, use that as 'optimal'
f_tol = 20;

% frequency of matrix for each method when method produces 'optimal' perm
method_A_freq = repmat({zeros(n, n)}, m, 1);

count = zeros(m, 1); % count for times each method produces 'optimal' perm

% plot the 'p_tol' highest non-diagonal values of each frequency matrix
p_tol = 250;

p_n = 1e4; % number of random permutations to test in brute force

% random permutations (using same random permutations for each test)
p = zeros(p_n, n);
for i = 1:p_n
    p(i, :) = randperm(n);
end

% testing
for i = 1:t
    % generate random n x n matrix with density
    A = full(spones(sprand(n, n, density)));

    % fill-in for non-BF methods
    fill_ins = zeros(1, m);
    for j = 1:m
        % reorder A
        % MarkowitzPivoting returns row and column permutations
        if isequal(methods{j}, @MarkowitzPivoting)
            [P, Q] = MarkowitzPivoting(A);
            A_perm = A(P, Q);
        else
            perm = methods{j}(A);
            A_perm = A(perm, perm);
        end

        % fill-in
        [L, U, P] = lu(A_perm);
        fill_ins(j) = nnz(L) - nnz(tril(P * A_perm));
    end
    methods_min_fill = min(fill_ins);
    
    %disp("Min fill-in found by methods: " + methods_min_fill)
    
    % brute force
    min_id = -1; % minimum fill-in index
    BF_min_fill = Inf; % minimum fill-in value
    for j = 1:p_n
        % reorder A
        perm = p(j, :);
        A_perm = A(perm, perm);
    
        % fill-in
        [L, U, P] = lu(A_perm);
        fill_in = nnz(L) - nnz(tril(P * A_perm));
    
        % if optimal permutation found (or one better than ordering
        % methods), break
        if fill_in <= f_tol || fill_in < methods_min_fill
            min_id = j;
            BF_min_fill = fill_in;
            break
        elseif fill_in < BF_min_fill
            min_id = j;
            BF_min_fill = fill_in;
        end
    end
    
    %disp("Min fill-in found by brute force: " + BF_min_fill)

    % adding A to frequency matrices only if optimal matrix found by BF
    if BF_min_fill <= methods_min_fill
        for j = 1:m
            fill_in = fill_ins(j);

            if fill_in == BF_min_fill
                method_A_freq{j} = method_A_freq{j} + A;
                count(j) = count(j) + 1;
            end
        end
    end

    % use to see how long your tests are taking (if BF is taking too long)
    if mod(i, t / 10) == 0
        disp(i + " tests done.")
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
