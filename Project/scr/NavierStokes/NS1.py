

import petsc4py
import sys
petsc4py.init(sys.argv)
from petsc4py import PETSc

from dolfin import *
# from MatrixOperations import *
import numpy as np
import os
import scipy.io
#from PyTrilinos import Epetra, EpetraExt, AztecOO, ML, Amesos
#from scipy2Trilinos import scipy_csr_matrix2CrsMatrix
import PETScIO as IO
import time
import common
import CheckPetsc4py as CP
import NSpreconditioner
from scipy.sparse import  spdiags
import MatrixOperations as MO
import ExactSol
import NSprecondSetup
import matplotlib.pylab as plt
# parameters["form_compiler"]["optimize"]     = True
# parameters["form_compiler"]["cpp_optimize"] = True
#MO.SwapBackend('epetra')
#os.system("echo $PATH")

m = 8
errL2u =np.zeros((m-1,1))
errH1u =np.zeros((m-1,1))
errL2p =np.zeros((m-1,1))

l2uorder =  np.zeros((m-1,1))
H1uorder =np.zeros((m-1,1))
l2porder =  np.zeros((m-1,1))
NN = np.zeros((m-1,1))
DoF = np.zeros((m-1,1))
Vdim = np.zeros((m-1,1))
Qdim = np.zeros((m-1,1))
Wdim = np.zeros((m-1,1))
l2uorder = np.zeros((m-1,1))
l2porder = np.zeros((m-1,1))
nonlinear = np.zeros((m-1,1))
SolTime = np.zeros((m-1,1))
AvIt = np.zeros((m-1,1))
NonLinearIts = np.zeros((m-1,1))
nn = 2

dim = 2
Solver = 'PCD'
Saving = 'no'
case = 4
# parameters['linear_algebra_backend'] = 'uBLAS'
# parameters = CP.ParameterSetup()
def LOG(arg):
    if INFO:
        print(arg)




for xx in xrange(1,m):
    print xx
    NN[xx-1] = xx+1
    nn = 2**(NN[xx-1])
    nn = int(nn)
    mesh = UnitSquareMesh(nn,nn)
    # tic()
    parameters["form_compiler"]["quadrature_degree"] = -1

    parameters['reorder_dofs_serial'] = False
    V = VectorFunctionSpace(mesh, "CG", 2)
    Q = FunctionSpace(mesh, "CG", 1)
    mpi_comm = mpi_comm_world()
    my_rank = MPI.rank(mpi_comm)
    num_threads = dolfin.parameters["num_threads"]
    parameters["num_threads"] = 0
    # if my_rank == 0:
    ones = Function(Q)
    ones.vector()[:]=(1)
    print ones.shape
    parameters["num_threads"] = num_threads
    # QQ = VectorFunctionSpace(mesh,"B",3)
    # V = V+QQ
    parameters['reorder_dofs_serial'] = False
    # print 'time to create function spaces', toc(),'\n\n'
    W = V*Q
    Vdim[xx-1] = V.dim()
    Qdim[xx-1] = Q.dim()
    Wdim[xx-1] = W.dim()
    print "\n\nV:  ",Vdim[xx-1],"Q:  ",Qdim[xx-1],"W:  ",Wdim[xx-1],"\n\n"

    def boundary(x, on_boundary):
        return on_boundary



    R = 10.0
    # MU = Constant(0.01)
    # MU = 1000.0
    MU = 1.0

    (u, p) = TrialFunctions(W)
    (v, q) = TestFunctions(W)




    class r0(Expression):
        def __init__(self):
            self.M = 1
        def eval_cell(self, values, x, ufc_cell):
            values[0] = 1.
            values[1] = 1.
        def value_shape(self):
            return (2,)
    f = r0()

    n = FacetNormal(mesh)
    h = CellSize(mesh)
    h_avg =avg(h)
    d = 0
    u_k,p_k = common.Stokes(V,Q,f,f,[1,1,MU], FS = "CG",InitialTol = 1e-6)
    # p_k.vector()[:] = p_k.vector().array()
    pConst = - assemble(p_k*dx)/assemble(ones*dx)
    p_k.vector()[:] += pConst
    # u_k = Function(V)
    # p_k = Function(Q)
    # plot(u_k)
    # plot(p_k)
    uOld = np.concatenate((u_k.vector().array(),p_k.vector().array()), axis=0)
    r = IO.arrayToVec(uOld)

    a11 = MU*inner(grad(v), grad(u))*dx(mesh) + inner((grad(u)*u_k),v)*dx(mesh) + (1./2)*div(u_k)*inner(u,v)*dx(mesh) - (1./2)*inner(u_k,n)*inner(u,v)*ds(mesh)
    a12 = div(v)*p*dx
    a21 = div(u)*q*dx
    # L1  = inner(v-0.1*h*h*grad(q), f)*dx
    print inner(v, f).shape()
    L1  = inner(v, f)*dx
    a = a11-a12-a21


    r11 = MU*inner(grad(v), grad(u_k))*dx(mesh) + inner((grad(u_k)*u_k),v)*dx(mesh) + (1./2)*div(u_k)*inner(u_k,v)*dx(mesh) - (1./2)*inner(u_k,n)*inner(u_k,v)*ds(mesh)
    r12 = div(v)*p_k*dx
    r21 = div(u_k)*q*dx
    RHSform = r11-r12-r21
    # -r22
    # RHSform = 0

    p11 = inner(u,v)*dx
    # p12 = div(v)*p*dx
    # p21 = div(u)*q*dx
    p22 = inner(p,q)*dx
    prec = p11 +p22
    bc = DirichletBC(W.sub(0),Expression(("0.0","0.0")), boundary)
    bcs = [bc]

    eps = 1.0           # error measure ||u-u_k||
    tol = 1.0E-4      # tolerance
    iter = 0            # iteration counter
    maxiter = 20        # max no of iterations allowed
    # parameters = CP.ParameterSetup()
    outerit = 0

    N = FacetNormal(mesh)
    h = CellSize(mesh)
    h_avg =avg(h)
    alpha = 10.0
    gamma =10.0

    (pQ) = TrialFunction(Q)
    (qQ) = TestFunction(Q)
    print MU
    Mass = assemble(inner(pQ,qQ)*dx)
    L = MU*(inner(grad(qQ), grad(pQ))*dx(mesh))

    O = inner(inner(grad(pQ),u_k),qQ)*dx(mesh)\
            - (1./2)*inner(u_k,n)*inner(qQ,pQ)*ds(mesh) \
            -(1/2)*(dot(u_k('+'),N('+'))+dot(u_k('-'),N('-')))*avg(inner(qQ,pQ))*ds(mesh) \
            -dot(avg(qQ),dot(outer(pQ('+'),N('+'))+outer(pQ('-'),N('-')),avg(u_k)))*dS(mesh)


    Laplacian = assemble(L)
    Laplacian = CP.Assemble(Laplacian)
    Mass = CP.Assemble(Mass)
    kspA, kspQ = NSprecondSetup.PCDKSPlinear(Mass,Laplacian)

    u_is = PETSc.IS().createGeneral(range(W.sub(0).dim()))
    p_is = PETSc.IS().createGeneral(range(W.sub(0).dim(),W.sub(0).dim()+W.sub(1).dim()))
    # print L
    SolutionTime = 0
    while eps > tol and iter < maxiter:
        iter += 1
        x = Function(W)

        uu = Function(W)
        tic()
        AA, bb = assemble_system(a, L1-RHSform, bcs)
        A,b = CP.Assemble(AA,bb)
        print toc()
        # b = b.getSubVector(t_is)
        F = A.getSubMatrix(u_is,u_is)

        kspF = NSprecondSetup.LSCKSPnonlinear(F)

        x = b.duplicate()
        ksp = PETSc.KSP()
        ksp.create(comm=PETSc.COMM_WORLD)
        # ksp.setTolerances(1e-5)
        ksp.setType('gmres')
        pc = ksp.getPC()


        pc.setType(PETSc.PC.Type.PYTHON)
        if Solver == "PCD":
            Fp = assemble(L+O)
            Fp = CP.Assemble(Fp)
            pc.setPythonContext(NSpreconditioner.NSPCD(W, kspF, kspA, kspQ,Fp))
        else:
            pc.setPythonContext(NSpreconditioner.NSLSC(W, kspF, kspBQB, dBt) )
        # elif Solver == "PCD":
            # F = assemble(fp)
            # F = CP.Assemble(F)
        #     pc.setPythonContext(NSprecond.PCD(W, A, Mass, F, L))

        ksp.setOperators(A)
        OptDB = PETSc.Options()
        ksp.max_it = 1000
        # OptDB['pc_factor_shift_amount'] = 1
        # OptDB['pc_factor_mat_ordering_type'] = 'rcm'
        # OptDB['pc_factor_mat_solver_package']  = 'mumps'
        ksp.setFromOptions()


        x = r

        toc()
        ksp.solve(b, x)

        time = toc()
        print time
        MO.StrTimePrint("Solve time",time)
        SolutionTime = SolutionTime +time
        outerit += ksp.its
        print "==============", ksp.its
        u1 = u.getSubVector(is0)
        p1 = u.getSubVector(is1)
        print u1.norm()
        print p1.norm()
        u2 = Function(V)
        u2.vector()[:] = u1.vector().array() + u_k.vector().array()
        p2 = Function(Q)
        p2.vector()[:] = p1.vector().array() + p_k.vector().array()
        p2.vector()[:] += - assemble(p2*dx)/assemble(ones*dx)
        u_k.assign(u2)
        p_k.assign(p2)

        # plot(p_k)

        uOld = np.concatenate((u_k.vector().array(),p_k.vector().array()), axis=0)
        r = IO.arrayToVec(uOld)
        # plot(p_k)



    # del  solver




import pandas as pd
# tableTitles = ["Total DoF","V DoF","Q DoF","AvIt","V-L2","V-order","P-L2","P-order"]
# tableValues = np.concatenate((Wdim,Vdim,Qdim,AvIt,errL2u,l2uorder,errL2p,l2porder),axis=1)
# df = pd.DataFrame(tableValues, columns = tableTitles)
# pd.set_option('precision',3)
# print df
# print df.to_latex()

# print "\n\n   Velocity convergence"
# VelocityTitles = ["Total DoF","V DoF","Soln Time","AvIt","V-L2","L2-order","V-H1","H1-order"]
# VelocityValues = np.concatenate((Wdim,Vdim,SolTime,AvIt,errL2u,l2uorder,errH1u,H1uorder),axis=1)
# VelocityTable= pd.DataFrame(VelocityValues, columns = VelocityTitles)
# pd.set_option('precision',3)
# VelocityTable = MO.PandasFormat(VelocityTable,"V-L2","%2.4e")
# VelocityTable = MO.PandasFormat(VelocityTable,'V-H1',"%2.4e")
# VelocityTable = MO.PandasFormat(VelocityTable,"H1-order","%1.2f")
# VelocityTable = MO.PandasFormat(VelocityTable,'L2-order',"%1.2f")
# print VelocityTable

# print "\n\n   Pressure convergence"
# PressureTitles = ["Total DoF","P DoF","Soln Time","AvIt","P-L2","L2-order"]
# PressureValues = np.concatenate((Wdim,Qdim,SolTime,AvIt,errL2p,l2porder),axis=1)
# PressureTable= pd.DataFrame(PressureValues, columns = PressureTitles)
# pd.set_option('precision',3)
# PressureTable = MO.PandasFormat(PressureTable,"P-L2","%2.4e")
# PressureTable = MO.PandasFormat(PressureTable,'L2-order',"%1.2f")
# print PressureTable

# print "\n\n"

# LatexTitles = ["l","DoFu","Dofp","V-L2","L2-order","V-H1","H1-order","P-L2","PL2-order"]
# LatexValues = np.concatenate((NN,Vdim,Qdim,errL2u,l2uorder,errH1u,H1uorder,errL2p,l2porder), axis=1)
# LatexTable = pd.DataFrame(LatexValues, columns = LatexTitles)
# pd.set_option('precision',3)
# LatexTable = MO.PandasFormat(LatexTable,"V-L2","%2.4e")
# LatexTable = MO.PandasFormat(LatexTable,'V-H1',"%2.4e")
# LatexTable = MO.PandasFormat(LatexTable,"H1-order","%1.2f")
# LatexTable = MO.PandasFormat(LatexTable,'L2-order',"%1.2f")
# LatexTable = MO.PandasFormat(LatexTable,"P-L2","%2.4e")
# LatexTable = MO.PandasFormat(LatexTable,'PL2-order',"%1.2f")
# print LatexTable.to_latex()


# print "\n\n\n\n"

# LatexTitles = ["l","DoFu","Dofp","Soln Time","AvIt","Non-Lin its"]
# LatexValues = np.concatenate((NN,Vdim,Qdim, SolTime,AvIt, NonLinearIts), axis=1)
# LatexTable = pd.DataFrame(LatexValues, columns = LatexTitles)
# pd.set_option('precision',3)
# LatexTable = MO.PandasFormat(LatexTable,'AvIt',"%3.1f")
# print LatexTable.to_latex()


# plot(ua)
# plot(interpolate(ue,V))

# plot(pp)
# plot(interpolate(p0,Q))

# interactive()

plt.show()

