function x = tridiagSparseSolve(a,b,c,f)
%TRIDIAGSPARSESOLVE solves the tridiagonal matrix represented by the
%vectors a, b and c. Nicked (with minor modifications) from
%https://people.sc.fsu.edu/~jburkardt/classes/math2071_2020/tridiagonal/tridiagonal.pdf

n = size(a,1);

for j = 1:(n-1)
    s = a(j+1) / b(j);
    b(j+1) = b(j+1) - s*c(j);
    f(j+1) = f(j+1) - s*f(j);
end

x = zeros(n,1);
x(n) = f(n)/b(n);
for j = (n-1):-1:1
    x(j) = (f(j) - c(j) * x(j+1)) / b(j);
end
