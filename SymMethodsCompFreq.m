% Produces a heatmap for frequencies of entries for when the first method
% produces a lower fill-in than the second method (only compares two
% methods).

% Used to see if there's a pattern for when one method is better than the
% other in terms of fill-in.

% Use plot tolerance 'p_tol' to ignore the 'n x n - 2(p_tol + n)' lowest 
% entries (filters highest values in frequency matrix to see shape of 
% matrices).
%   - e.g. if you want to plot the 100 most frequent non-diagonal values, 
%     set p_tol = 50.
% Note that the actual number of points plotted will have a minimum of 2 *
% p_tol but can have a high ceiling - if the lowest frequency number 
% that passes the threshold is very common in the frequency matrix.

n = 100; % square matrix dimensions
density = 0.1; % density of generated matrix
t = 50000; % test number

% only interested in when the first ordering method in the list performs 
% better than the second
methods = {@dissect, @symamd}; 
method_names = ["Dissection", "AMD"];

% frequency of matrix for first method when method performs better than the
% second
A_freq = zeros(n);

% count for times first method performs better than second
count = 0;

fill_ins = zeros(t, 2); % fill-ins for each method 

% plot the '2 * p_tol' highest non-diagonal values of frequency matrix
p_tol = 200;

% testing
for i = 1:t
    rng(i); % set seed for reproducibility

    % generate random n x n spd matrix with density
    A = sprandsym(n, density) + n * eye(n);

    % fill-in for each method
    fill_in_1 = Inf; % fill_in for first method
    for j = 1:2
        % reorder A
        perm = methods{j}(A);
        A_perm = A(perm, perm);

        % fill-in
        A_perm_chol = chol(A_perm, "lower");
        fill_in = nnz(A_perm_chol) - nnz(tril(A_perm));

        fill_ins(i, j) = fill_in;

        % comparing fill-ins
        if fill_in_1 < fill_in
            A_freq = A_freq + spones(A);
            count = count + 1;
        end
        fill_in_1 = fill_in;
    end
end

% plotting shapes of overall matrix as heatmap

% diagonal values are much larger than any other point (due to
% sparsity), set diagonals to maximum non-diagonal entry
max_value = max(1, max(A_freq(A_freq < count)));
A_freq(A_freq >= count) = max_value;

% lowest values that pass threshold
min_mv = min(maxk(A_freq(:), p_tol + n));

% ignore points that don't meet the threshold
A_freq(A_freq < min_mv) = 0;

% final matrix's heatmap plot
figure
heatmap(A_freq, "GridVisible", "off", "Colormap", copper, ...
    "CellLabelColor", "none", ...
    "ColorLimits", [max(min_mv - 1, 0) max_value])
title(method_names(1) + " frequency plot")

disp(method_names(1) + " was better than " + method_names(2) + " " + ...
    count + " times.")

% overlapping histogram plot of fill-ins for both methods
%figure
%bins = max(3, floor(15 * log10(n + 1) - 10)); % histogram bins
%histogram(fill_ins(:, 1), bins, "DisplayName", method_names(1));
%hold on
%histogram(fill_ins(:, 2), bins, "DisplayName", method_names(2));
%ylabel("No. permutations");
%xlabel("Fill-in");
%legend
