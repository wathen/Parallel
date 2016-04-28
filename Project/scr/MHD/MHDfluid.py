#!/usr/bin/python

# interpolate scalar gradient onto nedelec space
import petsc4py
import sys

petsc4py.init(sys.argv)

from petsc4py import PETSc
from dolfin import *
Print = PETSc.Sys.Print
# from MatrixOperations import *
import numpy as np
import PETScIO as IO
import common
import scipy
import scipy.io
import time
import scipy.sparse as sp

import BiLinear as forms
import IterOperations as Iter
import MatrixOperations as MO
import CheckPetsc4py as CP
import ExactSol
import Solver as S
import MHDmatrixPrecondSetup as PrecondSetup
import NSprecondSetup
import MHDprec as MHDpreconditioner
import memory_profiler
import gc
import MHDmulti
import MHDmatrixSetup as MHDsetup
#@profile
def Scipy2PETSc(A):
    A = A.tocsr()
    return PETSc.Mat().createAIJ(size=A.shape, csr=(A.indptr, A.indices, A.data))

def PETSc2Scipy(A):
    row, col, value = A.getValuesCSR()
    return sp.csr_matrix((value, col, row), shape=A.size)

def StoreMatrix(A,name):
    sA = A
    test ="".join([name,".mat"])
    scipy.io.savemat( test, {name: sA},oned_as='row')

def SaveHX(HX,level):
    StoreMatrix(PETSc2Scipy(HX[0]),'Mat/'+ str(level) +'A0')
    StoreMatrix(PETSc2Scipy(HX[1][0]),'Mat/'+ str(level) +'A10')
    StoreMatrix(PETSc2Scipy(HX[1][1]),'Mat/'+ str(level) +'A11')
    StoreMatrix(PETSc2Scipy(HX[2].getOperators()[0]),'Mat/'+ str(level) +'A2')
    StoreMatrix(PETSc2Scipy(HX[3].getOperators()[0]),'Mat/'+ str(level) +'A3')
    StoreMatrix(PETSc2Scipy(HX[4].getOperators()[0]),'Mat/'+ str(level) +'A4')
    # np.save(HX[5].array,'Mat/'+ str(level) +'A5')
    StoreMatrix(PETSc2Scipy(HX[6]),'Mat/'+ str(level) +'A6')

def SaveF(F,level):
    StoreMatrix(PETSc2Scipy(F[0].getOperators()[0]),'Mat/'+ str(level) +'F0')
    StoreMatrix(PETSc2Scipy(F[1].getOperators()[0]),'Mat/'+ str(level) +'F1')

def arrayToVec(vecArray):
    vec = PETSc.Vec().create(comm=PETSc.COMM_WORLD)
    vec.setSizes(len(vecArray))
    vec.setUp()
    (Istart,Iend) = vec.getOwnershipRange()
    return vec.createWithArray(vecArray[Istart:Iend],
            comm=PETSc.COMM_WORLD)
    vec.destroy()


def foo():
    m = 2


    errL2u =np.zeros((m-1,1))
    errH1u =np.zeros((m-1,1))
    errL2p =np.zeros((m-1,1))
    errL2b =np.zeros((m-1,1))
    errCurlb =np.zeros((m-1,1))
    errL2r =np.zeros((m-1,1))
    errH1r =np.zeros((m-1,1))



    l2uorder =  np.zeros((m-1,1))
    H1uorder =np.zeros((m-1,1))
    l2porder =  np.zeros((m-1,1))
    l2border =  np.zeros((m-1,1))
    Curlborder =np.zeros((m-1,1))
    l2rorder =  np.zeros((m-1,1))
    H1rorder = np.zeros((m-1,1))

    NN = np.zeros((m-1,1))
    DoF = np.zeros((m-1,1))
    Velocitydim = np.zeros((m-1,1))
    Magneticdim = np.zeros((m-1,1))
    Pressuredim = np.zeros((m-1,1))
    Lagrangedim = np.zeros((m-1,1))
    Wdim = np.zeros((m-1,1))
    iterations = np.zeros((m-1,1))
    SolTime = np.zeros((m-1,1))
    udiv = np.zeros((m-1,1))
    MU = np.zeros((m-1,1))
    level = np.zeros((m-1,1))
    NSave = np.zeros((m-1,1))
    Mave = np.zeros((m-1,1))
    TotalTime = np.zeros((m-1,1))

    nn = 2

    dim = 2
    ShowResultPlots = 'yes'
    split = 'Linear'

    MU[0]= 1e0
    for xx in xrange(1,m):
        # print xx
        level[xx-1] = xx + 6
        nn = 2**(level[xx-1])



        # Create mesh and define function space
        #nn = int(nn)
        #NN[xx-1] = nn/2
        # parameters["form_compiler"]["quadrature_degree"] = 6
        # parameters = CP.ParameterSetup()
        #mesh = UnitSquareMesh(nn,nn)
        
        N = int(sys.argv[1])
        n = int(2**N)
        mesh = UnitCubeMesh(n,n,n)
        order = 2
        parameters['reorder_dofs_serial'] = False
        Velocity = VectorFunctionSpace(mesh, "CG", order)
        Pressure = FunctionSpace(mesh, "CG", order-1)
        Magnetic = FunctionSpace(mesh, "N1curl", order-1)
        Lagrange = FunctionSpace(mesh, "CG", order-1)
        W = MixedFunctionSpace([Velocity, Pressure, Magnetic,Lagrange])
        # W = Velocity*Pressure*Magnetic*Lagrange
        Velocitydim[xx-1] = Velocity.dim()
        Pressuredim[xx-1] = Pressure.dim()
        Magneticdim[xx-1] = Magnetic.dim()
        Lagrangedim[xx-1] = Lagrange.dim()
        Wdim[xx-1] = W.dim()
        print "\n\nW:  ",Wdim[xx-1],"Velocity:  ",Velocitydim[xx-1],"Pressure:  ",Pressuredim[xx-1],"Magnetic:  ",Magneticdim[xx-1],"Lagrange:  ",Lagrangedim[xx-1],"\n\n"
        dim = [Velocity.dim(), Pressure.dim(), Magnetic.dim(), Lagrange.dim()]


        def boundary(x, on_boundary):
            return on_boundary




        # bc = [u0,p0,b0,r0]
        FSpaces = [Velocity,Pressure,Magnetic,Lagrange]


        (u, b, p, r) = TrialFunctions(W)
        (v, c, q, s) = TestFunctions(W)
        kappa = 10.0
        Mu_m =10.0
        MU = 1.0/1
        IterType = 'Full'
        Split = "No"
        Saddle = "No"
        Stokes = "No"
        SetupType = 'python-class'
        F_NS = Expression(("0.0", "0.0","0.0"))
        if kappa == 0:
            F_M = Expression(("0.0", "0.0","0.0"))
        else:
            F_M = Expression(("0.0", "0.0","0.0"))
        params = [kappa,Mu_m,MU]

        # MO.PrintStr("Seting up initial guess matricies",2,"=","\n\n","\n")
        # BCtime = time.time()
        # BC = MHDsetup.BoundaryIndices(mesh)
        # MO.StrTimePrint("BC index function, time: ", time.time()-BCtime)
        # Hiptmairtol = 1e-2
        # HiptmairMatrices = PrecondSetup.MagneticSetup(Magnetic, Lagrange, b0, r0, Hiptmairtol, params)


        # MO.PrintStr("Setting up MHD initial guess",5,"+","\n\n","\n\n")
        # u_k,p_k,b_k,r_k = common.InitialGuess(FSpaces,[u0,p0,b0,r0],[F_NS,F_M],params,HiptmairMatrices,1e-10,Neumann=Expression(("0","0")),options ="New")
        # b_t = TrialFunction(Velocity)
        # c_t = TestFunction(Velocity)

        # ones = Function(Pressure)
        # ones.vector()[:]=(0*ones.vector().array()+1)
        # # pConst = - assemble(p_k*dx)/assemble(ones*dx)
        # p_k.vector()[:] += - assemble(p_k*dx)/assemble(ones*dx)
        u_k = Function(Velocity)
        p_k = Function(Pressure)
        b_k = Function(Magnetic)
        r_k = Function(Lagrange)

        x = Iter.u_prev(u_k,p_k,b_k,r_k)
        # w = Function(W)
        # x =
        # KSPlinearfluids, MatrixLinearFluids = PrecondSetup.FluidLinearSetup(Pressure, MU)
        # kspFp, Fp = PrecondSetup.FluidNonLinearSetup(Pressure, MU, u_k)
        # #plot(b_k)

        ns,maxwell,CoupleTerm,Lmaxwell,Lns = forms.MHD2D(mesh, W,F_M,F_NS, u_k,b_k,params,IterType,"CG",Saddle,Stokes)
        RHSform = forms.PicardRHS(mesh, W, u_k, p_k, b_k, r_k, params,"CG",Saddle,Stokes)

        bcu = DirichletBC(W.sub(0),Expression(("0.0","0.0","0.0")), boundary)
        bcb = DirichletBC(W.sub(2),Expression(("0.0","0.0","0.0")), boundary)
        bcr = DirichletBC(W.sub(3),Expression(("0.0")), boundary)
        bcs = [bcu,bcb,bcr]

        # parameters['linear_algebra_backend'] = 'uBLAS'

        eps = 1.0           # error measure ||u-u_k||
        tol = 1.0E-4     # tolerance
        iter = 0            # iteration counter
        maxiter = 1       # max no of iterations allowed
        SolutionTime = 0
        outer = 0
        # parameters['linear_algebra_backend'] = 'uBLAS'

        # FSpaces = [Velocity,Magnetic,Pressure,Lagrange]

        if IterType == "CD":
            MO.PrintStr("Setting up PETSc "+SetupType,2,"=","\n","\n")
            Alin = MHDsetup.Assemble(W,ns,maxwell,CoupleTerm,Lns,Lmaxwell,RHSform,bcs+BC, "Linear",IterType)
            Fnlin,b = MHDsetup.Assemble(W,ns,maxwell,CoupleTerm,Lns,Lmaxwell,RHSform,bcs+BC, "NonLinear",IterType)
            A = Fnlin+Alin
            A,b = MHDsetup.SystemAssemble(FSpaces,A,b,SetupType,IterType)
            u = b.duplicate()


        u_is = PETSc.IS().createGeneral(range(Velocity.dim()))
        NS_is = PETSc.IS().createGeneral(range(Velocity.dim()+Pressure.dim()))
        M_is = PETSc.IS().createGeneral(range(Velocity.dim()+Pressure.dim(),W.dim()))
        OuterTol = 1e-5
        InnerTol = 1e-5
        NSits =0
        Mits =0
        TotalStart =time.time()
        SolutionTime = 0
        while eps > tol  and iter < maxiter:
            iter += 1
            # MO.PrintStr("Iter "+str(iter),7,"=","\n\n","\n\n")
            # AssembleTime = time.time()

            t0 = Timer("assemble_system")
            AA, bb = assemble_system(maxwell+ns+CoupleTerm, (Lmaxwell + Lns) - RHSform,  bcs)
            del(t0)
#            if IterType == "CD":
#                # MO.StrTimePrint("MHD CD RHS assemble, time: ", time.time()-AssembleTime)
#                b = MHDsetup.Assemble(W,ns,maxwell,CoupleTerm,Lns,Lmaxwell,RHSform,bcs+BC, "CD",IterType)
#            else:
#                # MO.PrintStr("Setting up PETSc "+SetupType,2,"=","\n","\n")
#                if  Split == "Yes":
#                    if iter == 1:
#                        Alin = MHDsetup.Assemble(W,ns,maxwell,CoupleTerm,Lns,Lmaxwell,RHSform,bcs+BC, "Linear",IterType)
#                        Fnlin,b = MHDsetup.Assemble(W,ns,maxwell,CoupleTerm,Lns,Lmaxwell,RHSform,bcs+BC, "NonLinear",IterType)
#                        A = Fnlin+Alin
#                        A,b = MHDsetup.SystemAssemble(FSpaces,A,b,SetupType,IterType)
#                        u = b.duplicate()
#                    else:
#                        Fnline,b = MHDsetup.Assemble(W,ns,maxwell,CoupleTerm,Lns,Lmaxwell,RHSform,bcs+BC, "NonLinear",IterType)
#                        A = Fnlin+Alin
#                        A,b = MHDsetup.SystemAssemble(FSpaces,A,b,SetupType,IterType)
#                else:
#                    AA, bb = assemble_system(maxwell+ns+CoupleTerm, (Lmaxwell + Lns) - RHSform,  bcs)
                    # A,b = CP.Assemble(AA,bb)
            # assemble(maxwell+ns+CoupleTerm)
            # if iter == 1:
            # MO.StrTimePrint("MHD total assemble, time: ", time.time()-AssembleTime)

            # u = b.duplicate()
            # kspFp, Fp = PrecondSetup.FluidNonLinearSetup(Pressure, MU, u_k)
            # print "Inititial guess norm: ",  u.norm(PETSc.NormType.NORM_INFINITY)

            # if IterType == 'Full':

            #     n = FacetNormal(mesh)
            #     mat =  as_matrix([[b_k[1]*b_k[1],-b_k[1]*b_k[0]],[-b_k[1]*b_k[0],b_k[0]*b_k[0]]])
            #     a = params[2]*inner(grad(b_t), grad(c_t))*dx(W.mesh()) + inner((grad(b_t)*u_k),c_t)*dx(W.mesh()) +(1./2)*div(u_k)*inner(c_t,b_t)*dx(W.mesh()) - (1./2)*inner(u_k,n)*inner(c_t,b_t)*ds(W.mesh())+kappa/Mu_m*inner(mat*b_t,c_t)*dx(W.mesh())
            #     ShiftedMass = assemble(a)
            #     bcu.apply(ShiftedMass)
            #     ShiftedMass = CP.Assemble(ShiftedMass)
            #     kspF = NSprecondSetup.LSCKSPnonlinear(ShiftedMass)
            # else:
            #     F = A.getSubMatrix(u_is,u_is)
            #     kspF = NSprecondSetup.LSCKSPnonlinear(F)

            # Options = 'p4'
            # for items in KSPlinearfluids:
            #     print(items)

            # SaveHX(HiptmairMatrices, int(level[xx-1][0]))
            # SaveF(KSPlinearfluids, int(level[xx-1][0]))
            # StoreMatrix(PETSc2Scipy(A),'Mat/' + str(int(level[xx-1][0])) +'K')
            # StoreMatrix(PETSc2Scipy(Fp),'Mat/' + str(int(level[xx-1][0])) + 'Fp')
            # StoreMatrix(PETSc2Scipy(kspF.getOperators()[0]),
            #     'Mat/' + str(int(level[xx-1][0])) +'F')
            # # b.array.tofile('Mat/'+ str(int(level[xx-1][0])) +'b.mat')
            # # x.array.tofile('Mat/'+ str(int(level[xx-1][0])) +'x.mat')
            # # np.array(dim).
            # np.save('Mat/'+ str(int(level[xx-1][0])) +'dim.mat',np.array(dim))
            # np.save('Mat/'+ str(int(level[xx-1][0])) +'b.mat',b.array)
            # np.save('Mat/'+ str(int(level[xx-1][0])) +'x.mat',x.array)
            # # np.save(b.array,'Mat/'+ str(level) +'/

    
    # interactive()
foo()

list_timings()
