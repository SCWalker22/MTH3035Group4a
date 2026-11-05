n = 4;

A = [1 2 4 0;...
    4 0 0 3;...
    0 2 0 8;...
    1 0 13 2];

p = perms(1:n); %Create all permutations

p(1, :)

[k, l] = size(p);
for i=1:k
    q = p(i, :);
    for j=1:k
        o = p(j, :);
        
        A_t = A(q, o);
        S = sparse(A_t);
        nnz(lu(S))
    end
end