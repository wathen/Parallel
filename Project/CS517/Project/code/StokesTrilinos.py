#!/usr/bin/python
from dolfin import *
import numpy as np
import matplotlib.pylab as plt
import os
import scipy.io
from PyTrilinos import Epetra, EpetraExt, AztecOO, ML, Amesos
from scipy2Trilinos import scipy_csr_matrix2CrsMatrix
def GradientCorrectSize(A,name,xdim,ydim):
     from PyTrilinos import EpetraExt
     from numpy import array,loadtxt
     import scipy.sparse as sps
     import scipy.io
     test ="".join([name,".txt"])
     EpetraExt.RowMatrixToMatlabFile(test,A)
     data = loadtxt(test)
     col,row,values = data[:,0]-1,data[:,1]-1,data[:,2]
     Asparse = sps.csr_matrix((values, (row, col)))
     As = Asparse[0:xdim,0:ydim]
     comm = Epetra.PyComm()
     Ap = scipy_csr_matrix2CrsMatrix(As, comm)

     return Ap

m = 5
errL2u = np.zeros((m-1,1))
errL2p = np.zeros((m-1,1))
l2uorder = np.zeros((m-1,1))
l2porder = np.zeros((m-1,1))
NN = np.zeros((m-1,1))
DoF = np.zeros((m-1,1))
Vdim = np.zeros((m-1,1))
Qdim = np.zeros((m-1,1))
Wdim = np.zeros((m-1,1))
iterations = np.zeros((m-1,1))
SolTime = np.zeros((m-1,1))
udiv = np.zeros((m-1,1))
nn = 2

dim = 2
Solving = 'Iterative'
ShowResultPlots = 'no'
ShowErrorPlots = 'no'
EigenProblem = 'no'
SavePrecond = 'no'
case = 1
parameters['linear_algebra_backend'] = 'Epetra'


for xx in xrange(1,m):
    print xx
    nn = 2**(xx)
    # Create mesh and define function space
    nn = int(nn)
    NN[xx-1] = nn

    if dim == 3:
        mesh = UnitCubeMesh(nn,nn,nn)
    else:
        mesh = UnitSquareMesh(nn,nn)

    #parameters['reorder_dofs_serial'] = False

    V =FunctionSpace(mesh, "BDM", 2 )
    Q =FunctionSpace(mesh, "DG",  1)
    Vdim[xx-1] = V.dim()
    print "\n\n DoF ", V.dim()
    # creating trial and test function s
    v = TestFunction(V)
    u = TrialFunction(V)
    p = TestFunction(Q)
    q = TrialFunction(Q)
    W = V*Q
    (u,p) = TestFunctions(W)
    (v,q) = TrialFunctions(W)
    def boundary(x, on_boundary):
        return on_boundary

    if dim == 3:
        u0 = Expression(("0","0","0"))
    else:
        u0 = Expression(("20*x[0]*pow(x[1],3)","5*pow(x[0],4)-5*pow(x[1],4)"))
        p0 = Expression("60*pow(x[0],2)*x[1]-20*pow(x[1],3)")

    bc = DirichletBC(W.sub(0),u0, boundary)

    if dim == 3:
        f = Expression(('- 2*(x[1]*x[1]-x[1])*(x[0]*x[0]-x[0])-2*(x[0]*x[0]-x[0])*(x[2]*x[2]-x[2])-2*(x[1]*x[1]-x[1])*(x[2]*x[2]-x[2])', \
        '- 2*(x[1]*x[1]-x[1])*(x[0]*x[0]-x[0])-2*(x[0]*x[0]-x[0])*(x[2]*x[2]-x[2])-2*(x[1]*x[1]-x[1])*(x[2]*x[2]-x[2])', \
        '- 2*(x[1]*x[1]-x[1])*(x[0]*x[0]-x[0])-2*(x[0]*x[0]-x[0])*(x[2]*x[2]-x[2])-2*(x[1]*x[1]-x[1])*(x[2]*x[2]-x[2])'))
    else:
        f = Expression(("0","0"))

    N = FacetNormal(mesh)
    n = N
    h = CellSize(mesh)
    h_avg =avg(h)
    alpha = 10.0
    gamma =10.0
    n = FacetNormal(mesh)
    h = CellSize(mesh)
    h_avg =avg(h)
    d = 0
    a11 = inner(grad(v), grad(u))*dx \
        - inner(avg(grad(v)), outer(u('+'),N('+'))+outer(u('-'),N('-')))*dS \
        - inner(outer(v('+'),N('+'))+outer(v('-'),N('-')), avg(grad(u)))*dS \
        + alpha/h_avg*inner(outer(v('+'),N('+'))+outer(v('-'),N('-')),outer(u('+'),N('+'))+outer(u('-'),N('-')))*dS \
        - inner(outer(v,N), grad(u))*ds \
        - inner(grad(v), outer(u,N))*ds \
        + gamma/h*inner(v,u)*ds
    a12 = div(v)*p*dx
    a21 = div(u)*q*dx
    L1  =  inner(v,f)*dx + gamma/h*inner(u0,v)*ds - inner(grad(v),outer(u0,N))*ds
    a = a11-a12-a21
    i = p*q*dx
    L1  =  inner(v,f)*dx + gamma/h*inner(u0,v)*ds - inner(grad(v),outer(u0,N))*ds

    Gr = inner(u,grad(q))*dx
    Laplacian = dot(grad(q), grad(p))*dx \
       - dot(avg(grad(q)), jump(p, N))*dS \
       - dot(jump(q, N), avg(grad(p)))*dS \
       + alpha/h_avg*dot(jump(q, N), jump(p, N))*dS \
       - dot(q*N, grad(p))*ds \
       - dot(grad(q), p*N)*ds \
       + gamma/h*q*p*ds

    L = assemble(Laplacian)
    gradient = assemble(Gr)
    tic()
    AA, bb = assemble_system(a, L1, bc)
    A_epetra = AA.down_cast().mat()
    print toc()
    b_epetra = bb.down_cast().vec()
    zeros = 0*bb.array()
    del bb
    x_epetra = Epetra.Vector(zeros)
    PP, Pb = assemble_system(a11+i,L1,bc)
    P_epetra = PP.down_cast().mat()
    if (Solving == 'Iterative'):

        MLList = {

            "default values" : "maxwell",
            "max levels" : 10,
            "output" : 10,
            "PDE equations" : 1,
            "increasing or decreasing" : "decreasing",
            "aggregation: type" : "Uncoupled-MIS",
            "aggregation: damping factor" : 1.3333,
            "coarse: max size" : 75,
            "aggregation: threshold" : 0.0,
            "smoother: sweeps" : 2,
            "smoother: damping factor" : 0.67,
            "smoother: type" : "MLS",
            "smoother: MLS polynomial order" : 4,
            "smoother: pre or post" : "both",
            "coarse: type" : "Amesos-KLU",
            "prec type" : "MGV",
            "print unused" : -2
        }

        MLList = {
            "default values" : "maxwell",
            "max levels"                                     : 10,
            "prec type"                                        : "MGV",
            "increasing or decreasing"               : "decreasing",
            "aggregation: type"                          : "Uncoupled-MIS",
            "aggregation: damping factor"         : 4.0/3.0,
            "eigen-analysis: type"                      : "cg",
            "eigen-analysis: iterations"              : 10,
            "smoother: sweeps"                          : 1,
            "smoother: damping factor"              : 1.0,
            "smoother: pre or post"                     : "both",
            "smoother: type"                               : "Hiptmair",
            "subsmoother: type"                         : "Chebyshev",
            "subsmoother: Chebyshev alpha"    : 27.0,
            "subsmoother: node sweeps"           : 4,
            "subsmoother: edge sweeps"           : 4,
            "PDE equations" : 1,
            "coarse: type"                                   : "Amesos-MUMPS",
            "coarse: max size"                           : 128

        }
        # MLList = {"max levels"        : 3,
        #   "output"            : 10,
        #   "smoother: type"    : "symmetric Gauss-Seidel",
        #   "aggregation: type" : "Uncoupled"
        #  }
        As = GradientCorrectSize(gradient.down_cast().mat(),"name",V.dim(),Q.dim())
        ML_Hiptmair = ML.MultiLevelPreconditioner(P_epetra,gradient.down_cast().mat(),L.down_cast().mat(),MLList)
        # ML_Hiptmair = ml.MultiLevelPreconditioner(A.down_cast().mat(),MLList)
        ML_Hiptmair.ComputePreconditioner()

        solver = AztecOO.AztecOO(A_epetra, x_epetra, b_epetra)
        solver.SetPrecOperator(ML_Hiptmair)
        solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_cg);
        solver.SetAztecOption(AztecOO.AZ_output, 50);
        err = solver.Iterate(155000, 1e-10)

        iterations[xx-1] = solver.NumIters()
        SolTime[xx-1] = solver.SolveTime()
        del solver, ML_Hiptmair

    if (Solving == 'Iterative' or Solving == 'Direct'):
        if dim == 3:
            ue = Expression(('x[0]*x[1]*x[2]*(x[1]-1)*(x[2]-1)*(x[0]-1)','x[0]*x[1]*x[2]*(x[1]-1)*(x[2]-1)*(x[0]-1)','x[0]*x[1]*x[2]*(x[1]-1)*(x[2]-1)*(x[0]-1)'))
        else:
            ue = Expression(("20*x[0]*pow(x[1],3)","5*pow(x[0],4)-5*pow(x[1],4)"))
            pe = Expression("60*pow(x[0],2)*x[1]-20*pow(x[1],3)")

        u = interpolate(ue,V)

        Nv  = u.vector().array().shape

        # X = IO.vecToArray(x)
        # x = X[0:Vdim[xx-1][0]]
        x = x_epetra[0:Nv[0]]
        ua = Function(V)
        ua.vector()[:] = x.array
        x = x_epetra[0:Nv[0]]
        ua = Function(V)
        ua.vector()[:] = x.array
        # udiv[xx-1] = assemble(div(ua)*dx)


        errL2u[xx-1] = errornorm(ue,ua,norm_type="L2", degree_rise=4,mesh=mesh)

        if xx == 1:
            l2uorder[xx-1] = 0
        else:
            l2uorder[xx-1] =  np.abs(np.log2(errL2u[xx-2]/errL2u[xx-1]))

        print errL2u[xx-1]


if (ShowErrorPlots == 'yes'):
    plt.loglog(NN,errL2u)
    plt.title('Error plot for CG2 elements - Velocity L2 convergence = %f' % np.log2(np.average((errL2u[0:m-2]/errL2u[1:m-1]))))
    plt.xlabel('N')
    plt.ylabel('L2 error')

    plt.show()


if (Solving == 'Iterative' or Solving == 'Direct'):
    print "\n\n"
    print "          ==============================="
    print "                  Results Table"
    print "          ===============================\n\n"
    import pandas as pd
    if (Solving == 'Iterative'):
        tableTitles = ["DoF","# iters","Soln Time","V-L2","V-order","||div u_h||"]
        tableValues = np.concatenate((Vdim,iterations,SolTime,errL2u,l2uorder,udiv),axis=1)
    elif (Solving == 'Direct'):
        tableTitles = ["DoF","Soln Time","V-L2","V-order","||div u_h||"]
        tableValues = np.concatenate((Vdim,SolTime,errL2u,l2uorder,udiv),axis=1)

    df = pd.DataFrame(tableValues, columns = tableTitles)
    pd.set_printoptions(precision=3)
    print df
    print "\n\n"
    print "Velocity Elements rate of convergence ", np.log2(np.average((errL2u[0:m-2]/errL2u[1:m-1])))
    print "Pressure Elements rate of convergence ", np.log2(np.average((errL2p[0:m-2]/errL2p[1:m-1])))
    print df.to_latex()


if (SavePrecond == 'yes'):
    scipy.io.savemat('eigenvalues/Wdim.mat', {'Wdim':Wdim-1},)


if (ShowResultPlots == 'yes'):
    plot(ua)
    plot(interpolate(ue,V))

    interactive()

