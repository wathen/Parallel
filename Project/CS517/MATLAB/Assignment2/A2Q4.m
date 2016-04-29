n = 10;

A = randn(n);

a1 = A(1,:);
a2 = A(2,:);

[V,D] = eig(A);
D = diag(D)
x1 = V(:,1);

A1 = A-x1*a1;
