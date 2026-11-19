% Note: gamma seems to be better fit than normal (large n)

n = 100; % square matrix dimensions
density =  2 / n; % density of generated matrix
t = 1000; % test number

% cell array of ordering method functions
methods = {@symrcm, @amd, @MarkowitzPivoting, @dissect}; 
method_names = ["RCM", "AMD", "Markowitz", "Dissection"];
m = length(methods);

% fill-ins for each method
fill_ins = zeros(t, m);

for i = 1:t
    % generate random n x n matrix with density
    A = full(spones(sprand(n, n, density)));

    % fill-in for each method
    for j = 1:m
        % reorder A
        perm = methods{j}(A);
        A_perm = A(perm, perm);

        % fill-in
        [L, U, P] = lu(A_perm);
        fill_ins(i, j) = nnz(L) + nnz(U) - nnz(P * A_perm);
    end
end

% plotting histogram with fitted normal distribution
figure
sgtitle("Permutation count for each fill-in value")

bins = floor(15 * log10(n + 1) - 10); % histogram bins
[min_fill_in, max_fill_in] = bounds(fill_ins(:)); % min and max fill-in overall

for i = 1:m
    subplot(2, ceil(m / 2), i)
    histfit(fill_ins(:, i), bins, "normal");
    subtitle(method_names(i))
    ylabel("Number of permutations");
    xlabel("Fill-in");
    xlim([min_fill_in, max_fill_in]);
end
