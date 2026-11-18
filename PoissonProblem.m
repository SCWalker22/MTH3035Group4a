% fully vectorised version of 'Poisson_problem.m'

m = 10; % grid size
del_x = 1 / m; % grid spacing

alpha = 2; % coefficient alpha
N = m * m; % number of matrix elements

ir = (mod(1:N, m) + floor((0:N - 1) / m) * m + 1)'; % periodic wrap, ir is i right
il = (mod(-1:N - 2, m) + floor((0:N - 1) / m) * m + 1)'; % periodic wrap, il is i left
ju = mod(m:N - 1 + m, N)' + 1; % periodic wrap, ju is j up
jd = mod(-m:N - 1 - m, N)' + 1; % periodic wrap, jd is j down


% f_{i+1} + f_{i-1} + f_{i+n} +f_{i-n} - (4 + \alpha^2 * \del{x^2})*f_{i} = \del{x^2}*r_{i}

rows = repelem(1:N, 5)';

cols_ = [(1:N)' ir il ju jd]';
cols = cols_(:);

val = -(4 + alpha ^ 2 * del_x ^ 2);
vals_ = [val 1 1 1 1]';
vals = repmat(vals_, N, 1);

A = sparse(rows, cols, vals, N, N);
spy(A)
%disp(full(A))

r = -del_x ^ 2 * randi([0 1], N, 1);  %ones(N,1);
f = A \ r;
%disp(f)
