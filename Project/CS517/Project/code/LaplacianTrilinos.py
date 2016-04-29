#!/usr/bin/python
from dolfin import *
import numpy as np
import matplotlib.pylab as plt
import os
import scipy.io
from PyTrilinos import Epetra, EpetraExt, AztecOO, ML, Amesos
from scipy2Trilinos import scipy_csr_matrix2CrsMatrix
from scipy.sparse import csr_matrix, dia_matrix
from scipy.sparse.linalg import spsolve

from pyamg.aggregation import smoothed_aggregation_solver
from pyamg.multilevel import multilevel_solver
from pyamg.relaxation.smoothing import change_smoothers
from pyamg.relaxation.relaxation import make_system
from pyamg.relaxation.relaxation import gauss_seidel
from pyamg.krylov._cg import cg


def hiptmair_smoother(A,x,b,D,iterations=4,sweep='symmetric'):
    A,x,b = make_system(A,x,b,formats=['csr','bsr'])
    gauss_seidel(A,x,b,iterations=4,sweep='forward')
    r = b-A*x
    x_G = np.zeros(D.shape[1])
    A_G,x_G,b_G = make_system(D.T*A*D,x_G,D.T*r,formats=['csr','bsr'])
    gauss_seidel(A_G,x_G,b_G,iterations=4,sweep='symmetric')
    x[:] += D*x_G
    gauss_seidel(A,x,b,iterations=4,sweep='backward')

def setup_hiptmair(lvl,iterations=4,sweep='symmetric'):
    D = lvl.D
    def smoother(A,x,b):
        hiptmair_smoother(A,x,b,D,iterations=iterations,sweep=sweep)
    return smoother

def edgeAMG(Anode,Acurl,D):
    nodalAMG = smoothed_aggregation_solver(Anode,max_coarse=10,keep=True)


    ##
    # construct multilevel structure
    levels = []
    levels.append( multilevel_solver.level() )
    levels[-1].A = Acurl
    levels[-1].D = D
    for i in range(1,len(nodalAMG.levels)):
        A = levels[-1].A
        Pnode = nodalAMG.levels[i-1].AggOp
        P = findPEdge(D, Pnode)
        R = P.T
        levels[-1].P = P
        levels[-1].R = R
        levels.append( multilevel_solver.level() )
        A = R*A*P
        D = csr_matrix(dia_matrix((1.0/((P.T*P).diagonal()),0),shape=(P.shape[1],P.shape[1]))*(P.T*D*Pnode))
        levels[-1].A = A
        levels[-1].D = D

    edgeML = multilevel_solver(levels)
    for i in range(0,len(edgeML.levels)):
        edgeML.levels[i].presmoother = setup_hiptmair(levels[i])
        edgeML.levels[i].postsmoother = setup_hiptmair(levels[i])
    return edgeML


def findPEdge ( D, PNode):
    ###
    # use D to find edges
    # each row has exactly two non zeros, a -1 marking the start node, and 1 marking the end node
    numEdges = D.shape[0]
    edges = np.zeros((numEdges,2))
    DRowInd = D.nonzero()[0]
    DColInd = D.nonzero()[1]
    for i in range(0,numEdges):
        if (1.0-0.0001 < np.abs(D[DRowInd[2*i],DColInd[2*i]]) and np.abs(D[DRowInd[2*i],DColInd[2*i]])  < 1.0+0.0001 ):  # first index is start, second is end
            edges[DRowInd[2*i],0] = DColInd[2*i]
            edges[DRowInd[2*i],1] = DColInd[2*i+1]
        else :  # first index is end, second is start
            edges[DRowInd[2*i],0] = DColInd[2*i+1]
            edges[DRowInd[2*i],1] = DColInd[2*i]

    ###
    # now that we have the edges, we need to find the nodal aggregates
    # the nodal aggregates are the columns


    aggs = PNode.nonzero()[1] # each row has 1 nonzero and that column is its aggregate
    numCoarseEdges = 0
    row = []
    col = []
    data = []
    coarseEdges = {}
    for i in range(0,edges.shape[0]):
        coarseV1 = aggs[edges[i,0]]
        coarseV2 = aggs[edges[i,1]]
        if ( coarseV1 != coarseV2 ): # this is a coarse edges
            #check if in dictionary
            if ( coarseEdges.has_key((coarseV1,coarseV2)) ):
                row.append(i)
                col.append(coarseEdges[(coarseV1,coarseV2)])
                data.append(1)
            elif ( coarseEdges.has_key((coarseV2,coarseV1))):
                row.append(i)
                col.append(coarseEdges[(coarseV2,coarseV1)])
                data.append(-1)
            else :
                coarseEdges[(coarseV1,coarseV2)] = numCoarseEdges
                numCoarseEdges = numCoarseEdges + 1
                row.append(i)
                col.append(coarseEdges[(coarseV1,coarseV2)])
                data.append(1)

    PEdge = csr_matrix( (data, (row,col) ),shape=(numEdges,numCoarseEdges) )
    return PEdge




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
     print xdim,ydim
     As = Asparse[xdim:,0:xdim]
     As.eliminate_zeros()
     As =  As.T
     comm = Epetra.PyComm()
     Ap = scipy_csr_matrix2CrsMatrix(As, comm)

     return Ap

def RemoveZeros(A,name  ):
     from PyTrilinos import EpetraExt
     from numpy import array,loadtxt
     import scipy.sparse as sps
     import scipy.io
     test ="".join([name,".txt"])
     EpetraExt.RowMatrixToMatlabFile(test,A)
     data = loadtxt(test)
     col,row,values = data[:,0]-1,data[:,1]-1,data[:,2]
     Asparse = sps.csr_matrix((values, (row, col)))
     Asparse.eliminate_zeros()
     comm = Epetra.PyComm()
     Ap = scipy_csr_matrix2CrsMatrix(Asparse, comm)

     return Ap

m = 10
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
parameters['linear_algebra_backend'] = 'uBLAS'


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

    V =FunctionSpace(mesh, "BDM",1 )
    Q =FunctionSpace(mesh, "DG",  1)
    Vdim[xx-1] = V.dim()
    print "\n\n DoF ", V.dim()
    # creating trial and test function s
    v = TestFunction(V)
    u = TrialFunction(V)
    p = TestFunction(Q)
    q = TrialFunction(Q)
    W = V*Q
    (u1,p1) = TestFunctions(W)
    (v1,q1) = TrialFunctions(W)
    def boundary(x, on_boundary):
        return on_boundary

    if dim == 3:
        u0 = Expression(("0","0","0"))
    else:
        u0 = Expression(('x[0]*x[1]','x[0]*x[1]'))

    bc = DirichletBC(V,u0, boundary)

    if dim == 3:
        f = Expression(('- 2*(x[1]*x[1]-x[1])*(x[0]*x[0]-x[0])-2*(x[0]*x[0]-x[0])*(x[2]*x[2]-x[2])-2*(x[1]*x[1]-x[1])*(x[2]*x[2]-x[2])', \
        '- 2*(x[1]*x[1]-x[1])*(x[0]*x[0]-x[0])-2*(x[0]*x[0]-x[0])*(x[2]*x[2]-x[2])-2*(x[1]*x[1]-x[1])*(x[2]*x[2]-x[2])', \
        '- 2*(x[1]*x[1]-x[1])*(x[0]*x[0]-x[0])-2*(x[0]*x[0]-x[0])*(x[2]*x[2]-x[2])-2*(x[1]*x[1]-x[1])*(x[2]*x[2]-x[2])'))
    else:
        f = Expression(('- 2*(x[1]*x[1]-x[1])-2*(x[0]*x[0]-x[0])','-2*(x[0]*x[0]-x[0]) - 2*(x[1]*x[1]-x[1])'))

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
    a = inner(grad(v), grad(u))*dx \
        - inner(avg(grad(v)), outer(u('+'),N('+'))+outer(u('-'),N('-')))*dS \
        - inner(outer(v('+'),N('+'))+outer(v('-'),N('-')), avg(grad(u)))*dS \
        + alpha/h_avg*inner(outer(v('+'),N('+'))+outer(v('-'),N('-')),outer(u('+'),N('+'))+outer(u('-'),N('-')))*dS \
        - inner(outer(v,N), grad(u))*ds \
        - inner(grad(v), outer(u,N))*ds \
        + gamma/h*inner(v,u)*ds

    L1  =  inner(v,f)*dx + gamma/h*inner(u0,v)*ds - inner(grad(v),outer(u0,N))*ds

    Gr = inner(u1,grad(q1))*dx
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
    A = AA.sparray()
    print toc()
    b = bb.array()
    x = 0*b
    # ml = smoothed_aggregation_solver(A,max_coarse=10)
    # residuals = []
    # x = ml.solve(b,tol=1e-10,accel='cg',residuals=residuals)
    Grad = gradient.sparray().T
    G = Grad[V.dim():,0:V.dim()]
    # G = gradient.sparray()[V.dim():,0:V.dim()]
    L = L.sparray()
    B = np.ones((V.dim(),1))


    from pyamg.gallery import poisson
    from pyamg.aggregation import smoothed_aggregation_solver,rootnode_solver,adaptive_sa_solver
    from pyamg import ruge_stuben_solver
    from pyamg.relaxation.smoothing import change_smoothers
    from pyamg.util.linalg import norm
    from scipy import rand, array, mean
    B = np.ones((A.shape[0],1))
    # ml = smoothed_aggregation_solver(A,B=B, max_coarse=100,strength=('symmetric', {'theta' : 0}),smooth=('energy', {'degree':6}),improve_candidates=('strength_based_schwarz',{'sweep':'symmetric','iterations': 4}),coarse_solver='splu')
    #
    ml = adaptive_sa_solver(A)
    #
    ## Set all levels to use gauss_seidel's defaults
    # smoothers = 'strength_based_strength_based_schwarz'
    # change_smoothers(ml, presmoother=smoothers, postsmoother=smoothers)
    # residuals=[]
    # x = ml.solve(b, tol=1e-8, residuals=residuals)
    # print len(residuals)
    # #
    # ## Set all levels to use three iterations of gauss_seidel's defaults
    # smoothers = ('strength_based_strength_based_schwarz', {'iterations' : 3})
    # change_smoothers(ml, presmoother=smoothers, postsmoother=None)
    # residuals=[]
    # x = ml.solve(b, tol=1e-8, residuals=residuals)
    # print len(residuals)

    #
    ## Set level 0 to use gauss_seidel's defaults, and all
    ## subsequent levels to use 5 iterations of cgnr
    # smoothers = ('kaczmarz_gauss_seidel',{'sweep':'symmetric','iterations': 4})
    # change_smoothers(ml, presmoother=smoothers, postsmoother=smoothers)
    residuals=[]
    # M = ml.aspreconditioner(cycle='V')
    print ml
    x, info = cg(A, b, tol=1e-8, maxiter=3000, M=M,residuals=residuals)
    print len(residuals)
    # ml = edgeAMG(L,A,G.T)
    # MLOp = ml.aspreconditioner()

    # r_edgeAMG = []
    # r_None = []
    # r_SA = []

    # ml_SA = smoothed_aggregation_solver(A)
    # ML_SAOP = ml_SA.aspreconditioner()
    # x_prec,info = cg(A,b,x,M=MLOp,tol=1e-8,residuals=r_edgeAMG)
    iterations[xx-1] = len(residuals)
    if (Solving == 'Iterative' or Solving == 'Direct'):
        if dim == 3:
            ue = Expression(('x[0]*x[1]*x[2]*(x[1]-1)*(x[2]-1)*(x[0]-1)','x[0]*x[1]*x[2]*(x[1]-1)*(x[2]-1)*(x[0]-1)','x[0]*x[1]*x[2]*(x[1]-1)*(x[2]-1)*(x[0]-1)'))
        else:
            ue = Expression(('(x[1]*x[1]-x[1])*(x[0]*x[0]-x[0])+x[0]*x[1]','(x[1]*x[1]-x[1])*(x[0]*x[0]-x[0])+x[0]*x[1]'))

        u = interpolate(ue,V)

        Nv  = u.vector().array().shape

        # X = IO.vecToArray(x)
        # x = X[0:Vdim[xx-1][0]]
        # x = x_epetra[0:Nv[0]]
        ua = Function(V)
        ua.vector()[:] = x
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

