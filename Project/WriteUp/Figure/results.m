A = [4.82e+00 ,  3.94e+01 ,  3.04e+02 ,  2.46e+03; 
3.08e+00 ,  2.09e+01 ,  1.58e+02 ,  1.30e+03 ; 
1.50e+00 ,  1.06e+01 ,  8.49e+01 ,  6.60e+02 ; 
7.88e-01 ,  5.86e+00 ,  4.65e+01 ,  3.58e+02 ; 
 6.72e-01 ,  3.33e+00 ,  2.40e+01 ,  1.88e+02];
P = [1 2 4 8 16];
for i = 1:4
    loglog(P,A(:,i));hold on
end