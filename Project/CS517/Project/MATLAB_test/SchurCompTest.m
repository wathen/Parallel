n = 10;
m = 5;

A = randn(10);
B = randn(n,m);
C = randn(m,n);

K = [A B;
     C sparse(m,m)];


L = [speye(n) sparse(m,n)';
     C/A, speye(m)];
U = [speye(n) A\B;
     sparse(m,n), speye(m)];

D = [A sparse(m,n)';
         sparse(m,n) -C*(A\B)];
     
norm(full(K)-full(L*D*U))

eig(full(K)),eig(full(D))
