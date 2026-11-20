% Notes:
%   - Use '''h = histfit(fill_ins(:, i), bins, "gamma"); 
%          h(2).Color = [1 0 0];''' to plot gamma fit on top
%   - MarkowitzPivoting slowing down higher n runs
%   - random bins calc function used, there's probably a simpler one.
%   - random quantities chosen as examples (for analysis)

n = 100; % square matrix dimensions
density =  2 / n; % density of generated matrix
t = 1000; % test number

% cell array of ordering method functions
methods = {@symrcm, @amd, @MarkowitzPivoting, @dissect}; 
method_names = ["RCM", "AMD", "Markowitz", "Dissection"];
m = length(methods);

% fill-ins for each method
fill_ins = zeros(t, m);

% testing
for i = 1:t
    % generate random n x n matrix with density
    A = full(spones(sprand(n, n, density)));

    % fill-in for each method
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
        fill_ins(i, j) = nnz(L) + nnz(U) - nnz(P * A_perm);
    end
end

% quantities of each method
disp("Fill-in quantities for each method" + newline)
for i = 1:m
    fill_in = fill_ins(:, i);
    [prcnt, edges] = histcounts(fill_in, 'Normalization', 'probability');
    disp(method_names(i) + newline + newline + ...
        "Bin percentage: " + num2str(prcnt) + newline + ...
        "Bin edges: " + num2str(edges) + newline + ...
        "Mean: " + mean(fill_in) + newline +  ...
        "Standard deviation: " + std(fill_in) + newline + ...
        "Minimum: " + min(fill_in) + newline + ...
        "Maximum: " + max(fill_in) + newline + ...
        "----------------------------------------------------" + newline)
end

% plotting histogram
figure
sgtitle("Permutation count for each fill-in value")

bins = floor(15 * log10(n + 1) - 10); % histogram bins
[min_fill_in, max_fill_in] = bounds(fill_ins(:)); % min and max fill-in overall

for i = 1:m
    subplot(ceil(m / 2), 2, i)
    histogram(fill_ins(:, i), bins)
    subtitle(method_names(i))
    ylabel("Number of permutations");
    xlabel("Fill-in");
    %xlim([min_fill_in, max_fill_in]); % useful to compare distributions
    xlim([min_fill_in, Inf]) % useful to see individual distributions
end
