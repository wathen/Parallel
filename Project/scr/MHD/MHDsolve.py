import petsc4py
import sys

petsc4py.init(sys.argv)

from petsc4py import PETSc
from scipy.sparse import csr_matrix
from scipy.io import loadmat
import numpy as np

def load_sparse_csr(filename):
    loader = np.load(filename)
    return csr_matrix((  loader['data'], loader['indices'], loader['indptr']),
                         shape = loader['shape'])
def Scipy2PETSc(A):
    A = A.tocsr()
    return PETSc.Mat().createAIJ(size=A.shape, csr=(A.indptr, A.indices, A.data))

def Multigrid(A):
    ksp = PETSc.KSP()
    ksp.create(comm=PETSc.COMM_WORLD)
    pc = ksp.getPC()
    ksp.setType('preonly')
    pc.setType('gamg')
    ksp.setOperators(A,A)
    return ksp

level = str(3)
HX = [Scipy2PETSc(loadmat('Mat/'+level+'A0.mat')['Mat/'+level+'A0']),
    [Scipy2PETSc(loadmat('Mat/'+level+'A10.mat')['Mat/'+level+'A10']),
    Scipy2PETSc(loadmat('Mat/'+level+'A11.mat')['Mat/'+level+'A11'])],
    Multigrid(Scipy2PETSc(loadmat('Mat/'+level+'A2.mat')['Mat/'+level+'A2'])),
    Scipy2PETSc(loadmat('Mat/'+level+'A3.mat')['Mat/'+level+'A3']),
    Scipy2PETSc(loadmat('Mat/'+level+'A4.mat')['Mat/'+level+'A4']),1,
    Scipy2PETSc(loadmat('Mat/'+level+'A6.mat')['Mat/'+level+'A6'])]

print HX