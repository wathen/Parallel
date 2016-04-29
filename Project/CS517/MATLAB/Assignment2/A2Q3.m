close all
A = [2 0 0 0 8 0 2 4 0 0;
     0 2 1 0 0 0 0 0 5 0;
     0 1 2 0 0 0 0 0 0 9;
     0 0 0 2 0 3 0 0 7 0;
     8 0 0 0 2 0 0 7 0 1;
     0 0 0 3 0 2 0 0 5 0;
     2 0 0 0 0 0 2 0 0 0;
     4 0 0 0 7 0 0 2 10 6;
     0 5 0 7 0 5 0 10 2 3;
     0 0 9 0 1 0 0 6 3 2;];

A= A+14*eye(10);
figure; spy(A)
pp = symrcm(A);
figure; spy(A(pp,pp))
r = chol(A(pp,pp));
figure;spy(r+r')


p = symamd(A);
figure; spy(A(p,p))
r = chol(A(p,p));
figure;spy(r+r')