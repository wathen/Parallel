m = 512;

e = ones(m,1);

A = spdiags([-e 2*e -e],[-1 0 1 ],m,m);
A = kron(speye(m),A)+kron(A,speye(m));

p = symrcm(A);
b = randn(m^2,1);

L = ichol(A,struct('droptol',1e-1));
pcg(A,b,1e-6,1000,L,L');
L = ichol(A(p,p),struct('droptol',1e-1));
pcg(A(p,p),b,1e-6,1000,L,L');
