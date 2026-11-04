% Solving the poisson problem

m = 3; % define grid size
del_x = 1/m; % define grid spacing

alpha = 2; % define coefficient alpha
N = m*m; %number of matrix elements


r = -del_x^2 * randi([0 1], N, 1)   %ones(N,1); 



rows = zeros(5*N,1); cols = zeros(5*N,1); vals = zeros(5*N,1);
k = 0;
for I= 1:N
    ir = mod(I,m) + floor((I-1)/m)*m + 1; % i+1, periodic wrap, ir is i right
    il = mod(I-2,m) + floor((I-1)/m)*m + 1; % i-1, periodic wrap, il is i left
    ju = mod(I-1+m,N) + 1; % j+1, periodic wrap, ju is j up
    jd = mod(I-1-m,N) + 1; % j-1, periodic wrap, jd is j down
       
    % f_{I+1} + f_{I-1} + f_{I+M} +f_{I-M} - (4 + \alpha^2 * \del{x^2})*f_{I} = \del{x^2}*r_{I}
    % forms the equation above by adding add to the kth row the coefficient
    % of f as vals(k) at the column which is the Ith index of f 
    k = k+1; rows(k) = I; cols(k) = I;  vals(k) = -(4 + alpha^2*del_x^2);
    k = k+1; rows(k) = I; cols(k) = ir; vals(k) = 1;
    k = k+1; rows(k) = I; cols(k) = il; vals(k) = 1;
    k = k+1; rows(k) = I; cols(k) = ju; vals(k) = 1;
    k = k+1; rows(k) = I; cols(k) = jd; vals(k) = 1;
end



A = sparse(rows, cols, vals, N, N);

spy(A)
fprintf('A is:\n');
disp(full(A))


f = A \ r;

fprintf('f is:\n');
disp(full(f))
