import petsc4py
import sys

petsc4py.init(sys.argv)

from petsc4py import PETSc
from scipy.sparse import csr_matrix
from scipy.io import loadmat
import numpy as np
import MHDprec
from dolfin import tic, toc

def Scipy2PETSc(A):
    A = A.tocsr()
    return PETSc.Mat().createAIJ(size=A.shape, csr=(A.indptr, A.indices, A.data))
def arrayToVec(vecArray):
    vec = PETSc.Vec().create(comm=PETSc.COMM_WORLD)
    vec.setSizes(len(vecArray))
    vec.setUp()
    (Istart,Iend) = vec.getOwnershipRange()
    return vec.createWithArray(vecArray[Istart:Iend],
            comm=PETSc.COMM_WORLD)
    vec.destroy()


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
kspF = Multigrid(Scipy2PETSc(loadmat('Mat/'+level+'F.mat')['Mat/'+level+'F']))
K = Scipy2PETSc(loadmat('Mat/'+level+'K.mat')['Mat/'+level+'K'])

dim = np.load('Mat/'+level+'dim.mat.npy')
b = arrayToVec(np.load('Mat/'+level+'b.mat.npy'))
x = arrayToVec(np.load('Mat/'+level+'x.mat.npy'))


ksp = PETSc.KSP()
ksp.create(comm=PETSc.COMM_WORLD)
pc = ksp.getPC()
ksp.setType('fgmres')
pc.setType('python')
pc.setType(PETSc.PC.Type.PYTHON)
OptDB = PETSc.Options()
OptDB['ksp_gmres_restart'] = 200
# FSpace = [Velocity,Magnetic,Pressure,Lagrange]
reshist = {}
def monitor(ksp, its, fgnorm):
    reshist[its] = fgnorm
    print its,"    OUTER:", fgnorm
# ksp.setMonitor(monitor)
ksp.max_it = 500
# W = Fspace
# FFSS = [W.sub(0),W.sub(1),W.sub(2),W.sub(3)]
pc.setPythonContext(MHDprec.Parallel(dim,kspF, Fluid[0], Fluid[1],Fp, HX[3], HX[4], HX[2], HX[0], HX[1], HX[6],1e-6,1))
#OptDB = PETSc.Options()

# OptDB['pc_factor_mat_solver_package']  = "mumps"
# OptDB['pc_factor_mat_ordering_type']  = "rcm"
# ksp.setFromOptions()
scale = b.norm()
b = b/scale
u = b.duplicate()
ksp.setOperators(K,K)
del K
tic()
ksp.solve(b,u)
print toc()
# Mits +=dodim
u = u*scale
print ksp.its
# MO.PrintStr("Number iterations = "+str(ksp.its),60,"+","\n\n","\n\n")


