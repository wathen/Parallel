from numpy import *
from scipy import *
from scipy.sparse import *
from scipy.linalg import norm
import time
import ipdb


def ConvectionDiffusion(beta,gamma,n):

    # Creating h
    h = 1.0/(n+1)
    e = ones((1,n))

    # Creating sparse diagonal matrices
    I = spdiags(e,0,n,n)
    I1 =spdiags(e,1,n,n)
    I2 = spdiags(e,-1,n,n)

    # Creating 1D Convection-Diffusion matricies
    Abeta = 2*I +(beta-1)*I1 - (beta+1)*I2
    Agamma = 2*I +(gamma-1)*I1 - (gamma+1)*I2

    # Creating 2D Convection-Diffusion matrix
    A = kron(Abeta,I)+kron(I,Agamma)

    return A



def MatVec(A,x):
    col,row = A.nonzero()
    nn,mm = x.shape
    A = A.tocsr()
    y = zeros((nn,1))
    for i in xrange(0,len(col)):
        y[row[i],0] = y[row[i],0]+A[row[i],col[i]]*x[row[i]]
    return y






beta = 1
gamma = 1

n = 40
start_time = time.time()
A = ConvectionDiffusion(beta,gamma,n)
print time.time()-start_time
x = random.rand(n**2,1)
# print A.todense()
start_time = time.time()
y = MatVec(A,x)
print time.time()-start_time

start_time = time.time()
y1 = A*x
print time.time()-start_time


print linalg.norm(y-y1)