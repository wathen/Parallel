import petsc4py
import sys

petsc4py.init(sys.argv)

from petsc4py import PETSc
from scipy.sparse import csr_matrix
from scipy.io import loadmat
import numpy as np

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

def MultigridCG(A):
    ksp = PETSc.KSP()
    ksp.create(comm=PETSc.COMM_WORLD)
    pc = ksp.getPC()
    ksp.setType('cg')
    pc.setType('gamg')
    ksp.setOperators(A,A)
    return ksp

def JacobiCG(A):
    ksp = PETSc.KSP()
    ksp.create(comm=PETSc.COMM_WORLD)
    pc = ksp.getPC()
    ksp.setType('cg')
    pc.setType('jacobi')
    ksp.setOperators(A,A)
    return ksp

level = str(3)
HX = [Scipy2PETSc(loadmat('Mat/'+level+'A0.mat')['Mat/'+level+'A0']),
    [Scipy2PETSc(loadmat('Mat/'+level+'A10.mat')['Mat/'+level+'A10']),
    Scipy2PETSc(loadmat('Mat/'+level+'A11.mat')['Mat/'+level+'A11'])],
    Multigrid(Scipy2PETSc(loadmat('Mat/'+level+'A2.mat')['Mat/'+level+'A2'])),
    Multigrid(Scipy2PETSc(loadmat('Mat/'+level+'A3.mat')['Mat/'+level+'A3'])),
    MultigridCG(Scipy2PETSc(loadmat('Mat/'+level+'A4.mat')['Mat/'+level+'A4'])),
    1,
    Scipy2PETSc(loadmat('Mat/'+level+'A6.mat')['Mat/'+level+'A6'])]

Fluid = [Multigrid(Scipy2PETSc(loadmat('Mat/'+level+'F0.mat')['Mat/'+level+'F0'])), JacobiCG(Scipy2PETSc(loadmat('Mat/'+level+'F1.mat')['Mat/'+level+'F1']))]

Fp = Scipy2PETSc(loadmat('Mat/'+level+'Fp.mat')['Mat/'+level+'Fp'])
F = Multigrid(Scipy2PETSc(loadmat('Mat/'+level+'F.mat')['Mat/'+level+'F']))
K = Scipy2PETSc(loadmat('Mat/'+level+'K.mat')['Mat/'+level+'K'])





