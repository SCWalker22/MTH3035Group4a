n = 4;

A = eye(n);

p = perms(1:n)

p(1, :)

[k, l] = size(p)
for i=1:k
    q = p(i, :)
    for j=1:k
        o = p(j, :)
        
        A_t = A(q, o);
        S = sparse(A_t);
        lu(S)
    end
end