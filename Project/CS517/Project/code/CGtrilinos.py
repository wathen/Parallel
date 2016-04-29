#!/usr/bin/python
from dolfin import *
import numpy as np
import matplotlib.pylab as plt
import os
import scipy.io
from PyTrilinos import Epetra, EpetraExt, AztecOO, ML, Amesos


m = 6
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

dim = 3
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

    parameters['reorder_dofs_serial'] = False
    V =VectorFunctionSpace(mesh, "CG", 1 )
    Vdim[xx-1] = V.dim()
    print "\n\n DoF ", V.dim()
    # creating trial and test function s
    v = TestFunction(V)
    u = TrialFunction(V)
    def boundary(x, on_boundary):
        return on_boundary

    if dim == 3:
        u0 = Expression(("0","0","0"))
    else:
        u0 = Expression(("0","0"))


    bc = DirichletBC(V,u0, boundary)

    if dim == 3:
        f = Expression(('- 2*(x[1]*x[1]-x[1])*(x[0]*x[0]-x[0])-2*(x[0]*x[0]-x[0])*(x[2]*x[2]-x[2])-2*(x[1]*x[1]-x[1])*(x[2]*x[2]-x[2])', \
        '- 2*(x[1]*x[1]-x[1])*(x[0]*x[0]-x[0])-2*(x[0]*x[0]-x[0])*(x[2]*x[2]-x[2])-2*(x[1]*x[1]-x[1])*(x[2]*x[2]-x[2])', \
        '- 2*(x[1]*x[1]-x[1])*(x[0]*x[0]-x[0])-2*(x[0]*x[0]-x[0])*(x[2]*x[2]-x[2])-2*(x[1]*x[1]-x[1])*(x[2]*x[2]-x[2])'))
    else:
         f = Expression(('- 2*(x[1]*x[1]-x[1])-2*(x[0]*x[0]-x[0])','-2*(x[0]*x[0]-x[0]) - 2*(x[1]*x[1]-x[1])'))


    a = inner(grad(v), grad(u))*dx
    L1  =  inner(v,f)*dx
    tic()
    AA, bb = assemble_system(a, L1, bc)
    A_epetra = AA.down_cast().mat()
    print toc()
    b_epetra = bb.down_cast().vec()
    zeros = 0*bb.array()
    del bb
    x_epetra = Epetra.Vector(zeros)

    if (Solving == 'Iterative'):

        MLList = {"max levels"        : 10,
          "output"            : 10,
          "smoother: type"    : "symmetric Gauss-Seidel",
          "aggregation: type" : "Uncoupled"
         }
        # ML_Hiptmair = ML.MultiLevelPreconditioner(P_epetra,gradient.down_cast().mat(),L.down_cast().mat(),MLList)
        ML_Hiptmair = ML.MultiLevelPreconditioner(A_epetra,False)
        ML_Hiptmair.SetParameterList(MLList)
        ML_Hiptmair.ComputePreconditioner()

        solver = AztecOO.AztecOO(A_epetra, x_epetra, b_epetra)
        solver.SetPrecOperator(ML_Hiptmair)
        solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_cg);
        solver.SetAztecOption(AztecOO.AZ_output, 1);
        err = solver.Iterate(155000, 1e-10)

        iterations[xx-1] = solver.NumIters()
        SolTime[xx-1] = solver.SolveTime()
        print iterations[xx-1]
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
        tableTitles = ["DoF","# iters","Soln Time"]
        tableValues = np.concatenate((Vdim,iterations,SolTime),axis=1)
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

