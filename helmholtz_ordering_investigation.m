function results = helmholtz_ordering_investigation()
% HELMHOLTZ_ORDERING_INVESTIGATION
% Compare sparsity orderings for a 2-D finite-difference Helmholtz problem.
%
% The domain is the unit square with a 10-by-10 grid of interior unknowns.
% The resulting coefficient matrix is 100-by-100.  Methods in this analysis:
%   1. Natural ordering
%   2. Reverse Cuthill-McKee (symrcm)
%   3. Exact minimum degree (implemented below)
%   4. Symmetric approximate minimum degree (symamd)
%   5. Pattern-only Markowitz pivoting (implemented below)
%
% It reports factor-pattern nonzeros, fill, bandwidth and scalar updates.
% It does not report wall-clock timing or numerical pivot stability.

clc;

%% 1. Build the finite-difference Helmholtz matrix
m = 10;                         % 10 x 10 interior grid = 100 unknowns
N = m^2;
h = 1/(m + 1);                 % grid spacing on the unit square
k = 4;                          % wavenumber; chosen below first discrete eigenfrequency

one = ones(m,1);
T = spdiags([-one, 2*one, -one], [-1, 0, 1], m, m);
I = speye(m);

% Discretisation of -Delta(u) - k^2*u = f with homogeneous Dirichlet BCs:
% A = (I kron T + T kron I)/h^2 - k^2 I.
A = (kron(I,T) + kron(T,I))/(h^2) - (k^2)*speye(N);
A = sparse(A);
S = spones(A);                  % structural pattern: 1 = nonzero, 0 = zero

fprintf('\n2-D Helmholtz finite-difference model\n');
fprintf('Grid: %d x %d interior points\n', m, m);
fprintf('Unknowns: %d\n', N);
fprintf('Matrix size: %d x %d\n', size(A,1), size(A,2));
fprintf('Structural nonzeros: %d of %d (%.2f%% density)\n', ...
    nnz(S), N^2, 100*nnz(S)/(N^2));
fprintf('Grid spacing h = %.6f, wavenumber k = %.3f\n', h, k);
fprintf('Symmetry check ||A-A''||_inf = %.3e\n', norm(A-A','inf'));

%% 2. Obtain the four non-natural orderings
pNatural = 1:N;
pRCM = symrcm(S);
pMD = exact_minimum_degree(S);
pAMD = symamd(S);
[rMark, cMark, markowitzTrace] = markowitz_pattern_order(S);

names = {'Natural'; 'Reverse Cuthill-McKee'; 'Exact minimum degree'; ...
         'Approximate minimum degree'; 'Markowitz pivoting'};
rPerms = {pNatural, pRCM, pMD, pAMD, rMark};
cPerms = {pNatural, pRCM, pMD, pAMD, cMark};

%% 3. Symbolically eliminate each reordered matrix
nMethods = numel(names);
stats = cell(nMethods,1);
Lpatterns = cell(nMethods,1);
Upatterns = cell(nMethods,1);

for q = 1:nMethods
    [stats{q}, Lpatterns{q}, Upatterns{q}] = ...
        symbolic_lu_stats(S, rPerms{q}, cPerms{q});
end

% Put the results in a MATLAB table.
Bandwidth = zeros(nMethods,1);
FillPositions = zeros(nMethods,1);
LNonzeros = zeros(nMethods,1);
UNonzeros = zeros(nMethods,1);
FactorNonzeros = zeros(nMethods,1);
ScalarUpdates = zeros(nMethods,1);
FlopEstimate = zeros(nMethods,1);

for q = 1:nMethods
    Bandwidth(q) = stats{q}.Bandwidth;
    FillPositions(q) = stats{q}.FillPositions;
    LNonzeros(q) = stats{q}.LNonzeros;
    UNonzeros(q) = stats{q}.UNonzeros;
    FactorNonzeros(q) = stats{q}.FactorNonzeros;
    ScalarUpdates(q) = stats{q}.ScalarUpdates;
    FlopEstimate(q) = stats{q}.FlopEstimate;
end

results = table(names, Bandwidth, FillPositions, LNonzeros, UNonzeros, ...
    FactorNonzeros, ScalarUpdates, FlopEstimate);

fprintf('\nSymbolic Gaussian-elimination comparison\n');
disp(results);

fprintf('Initial matrix nonzeros: %d\n', nnz(S));
fprintf(['FactorNonzeros = nnz(L) + nnz(U) - N. FillPositions counts ', ...
         'new structural entries introduced during elimination.\n']);
fprintf(['FlopEstimate = 2*(number of scalar Schur-complement updates) ', ...
         '+ number of multipliers); it is an estimate, not a timed flop count.\n']);

%% 4. Display the reordered sparsity patterns
figure('Name','Helmholtz matrix after reordering','Color','w');
for q = 1:nMethods
    subplot(2,3,q);
    spy(S(rPerms{q},cPerms{q}));
    title(names{q}, 'Interpreter','none');
    xlabel(sprintf('nnz = %d', nnz(S)));
end
sgtitle('100-by-100 Helmholtz matrix: sparsity patterns after reordering');

%% 5. Display the symbolic L and U patterns
figure('Name','Symbolic LU factor patterns','Color','w');
for q = 1:nMethods
    subplot(2,nMethods,q);
    spy(Lpatterns{q});
    title([names{q} ' : L'], 'Interpreter','none');
    subplot(2,nMethods,q+nMethods);
    spy(Upatterns{q});
    title([names{q} ' : U'], 'Interpreter','none');
end
sgtitle('Symbolic LU factors after each ordering');

%% 6. Show the pivot choices made by the Markowitz heuristic
MarkowitzPivotHistory = array2table(markowitzTrace, ...
    'VariableNames', {'OriginalRow','OriginalColumn','MarkowitzScore','NewFill'});
fprintf('\nFirst 15 Markowitz pivots (indices refer to original A):\n');
disp(MarkowitzPivotHistory(1:min(15,height(MarkowitzPivotHistory)),:));

% Save a copy of the comparison table for the project report.
writetable(results, 'helmholtz_ordering_results.csv');
fprintf('Wrote helmholtz_ordering_results.csv to the current folder.\n');

end


function p = exact_minimum_degree(S)

n = size(S,1);
G = logical(full(spones(S) + spones(S')));
G(1:n+1:end) = false;
active = true(n,1);
p = zeros(1,n);

for step = 1:n
    degree = sum(G,2);
    degree(~active) = inf;
    [~,v] = min(degree);       % exact degree in the current filled graph
    p(step) = v;

    neighbours = find(G(v,:) & active.');
    if numel(neighbours) > 1
        G(neighbours,neighbours) = true;  % eliminate v: connect its neighbours
        G(1:n+1:end) = false;
    end

    G(v,:) = false;
    G(:,v) = false;
    active(v) = false;
end
end


function [rperm,cperm,trace] = markowitz_pattern_order(S)
% Pattern-only Markowitz pivoting.
% For each nonzero candidate a(i,j) in the active submatrix, use
% (row_nnz(i)-1)*(column_nnz(j)-1) as the Markowitz fill estimate.
% The selected row and column are moved to the next pivot position.

B = logical(full(S));
n = size(B,1);
rperm = 1:n;
cperm = 1:n;
trace = zeros(n,4);

for step = 1:n
    C = B(step:n,step:n);
    rowCount = sum(C,2);
    colCount = sum(C,1);
    [ii,jj] = find(C);

    if isempty(ii)
        error('Markowitz ordering failed: no structural pivot in the active submatrix.');
    end

    score = (rowCount(ii)-1) .* (colCount(jj).'-1);
    [bestScore,which] = min(score);
    pivotRow = step + ii(which) - 1;
    pivotCol = step + jj(which) - 1;

    % Apply the pivot's independent row and column permutations.
    if pivotRow ~= step
        B([step pivotRow],:) = B([pivotRow step],:);
        rperm([step pivotRow]) = rperm([pivotRow step]);
    end
    if pivotCol ~= step
        B(:,[step pivotCol]) = B(:,[pivotCol step]);
        cperm([step pivotCol]) = cperm([pivotCol step]);
    end

    % Symbolic Schur-complement update: every row with a nonzero in the
    % pivot column is connected to every column with a nonzero in the pivot row.
    rows = find(B(step+1:n,step)) + step;
    cols = find(B(step,step+1:n)) + step;
    newFill = 0;
    if ~isempty(rows) && ~isempty(cols)
        block = B(rows,cols);
        newFill = nnz(~block);
        B(rows,cols) = true;
    end
    B(step+1:n,step) = false;

    trace(step,:) = [rperm(step), cperm(step), bestScore, newFill];
end
end


function [stats,Lpat,Upat] = symbolic_lu_stats(S,rperm,cperm)
% Count structural fill and symbolic LU factor entries for a fixed ordering.
% The diagonal is counted once in FactorNonzeros = nnz(L)+nnz(U)-n.

n = size(S,1);
B = logical(full(S(rperm,cperm)));
Lpat = logical(eye(n));
Upat = logical(eye(n));
fillPositions = 0;
scalarUpdates = 0;
flopEstimate = 0;

for step = 1:n
    rows = find(B(step+1:n,step)) + step;
    cols = find(B(step,step+1:n)) + step;

    Lpat(rows,step) = true;
    Upat(step,cols) = true;

    nRows = numel(rows);
    nCols = numel(cols);
    scalarUpdates = scalarUpdates + nRows*nCols;
    flopEstimate = flopEstimate + 2*nRows*nCols + nRows;

    if nRows > 0 && nCols > 0
        block = B(rows,cols);
        fillPositions = fillPositions + nnz(~block);
        B(rows,cols) = true;
    end

    % Entries below the pivot are now stored in L, not in the active matrix.
    B(step+1:n,step) = false;
end

[ii,jj] = find(S(rperm,cperm));
if isempty(ii)
    bandwidth = 0;
else
    bandwidth = max(abs(ii-jj));
end

stats.Bandwidth = bandwidth;
stats.FillPositions = fillPositions;
stats.LNonzeros = nnz(Lpat);
stats.UNonzeros = nnz(Upat);
stats.FactorNonzeros = stats.LNonzeros + stats.UNonzeros - n;
stats.ScalarUpdates = scalarUpdates;
stats.FlopEstimate = flopEstimate;
end
