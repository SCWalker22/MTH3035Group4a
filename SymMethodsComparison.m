% Notes:
%   - normal dist good fit (for large n) (gamma is also okay)
%   - MarkowitzPivoting slowing down higher n runs
%   - random bins calc function used, probably a simpler one.

n = 100; % square matrix dimensions
density =  2 / n; % density of generated matrix
t = 1000; % test number

% cell array of ordering method functions
methods = {@symrcm, @symamd, @MarkowitzPivoting, @dissect}; 
method_names = ["RCM", "AMD", "Markowitz", "Dissection"];
m = length(methods);

% fill-ins for each method
fill_ins = zeros(t, m);

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
        fill_ins(i, j) = nnz(A_perm_chol) - nnz(tril(A_perm));
    end
end

% plotting histogram with fitted normal distribution
figure
sgtitle("Permutation count for each fill-in value")
bins = floor(15 * log10(n + 1) - 10); % histogram bins
for i = 1:m
    subplot(2, ceil(m / 2), i)
    histfit(fill_ins(:, i), bins);
    subtitle(method_names(i))
    ylabel("Number of permutations");
    xlabel("Fill-in");
end
