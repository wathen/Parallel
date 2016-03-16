import petsc4py
import sys
petsc4py.init(sys.argv)
from petsc4py import PETSc
from dolfin import *
import MatrixOperations as MO
import numpy as np

parameters['linear_algebra_backend'] = 'PETSc'
parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "quadrature"
parameters['reorder_dofs_serial'] = False


class Approx(object):

    def __init__(self, Fp, Mp, Lp):
        self.Fp = Fp
        self.Mp = Mp
        self.Lp = Lp

    def create(self, pc):
        kspL = PETSc.KSP()
        kspL.create(comm=PETSc.COMM_WORLD)
        pcL = kspL.getPC()
        kspL.setType('preonly')
        pcL.setType('gamg')
        kspL.setFromOptions()
        kspL.setOperators(self.Lp, self.Lp)
        self.kspL = kspL

        kspM = PETSc.KSP()
        kspM.create(comm=PETSc.COMM_WORLD)
        pcM = kspM.getPC()
        kspM.setType('cg')
        pcM.setType('jacobi')
        kspM.setFromOptions()
        kspM.setOperators(self.Mp, self.Mp)
        self.kspM = kspM

    def apply(self, pc, x, y):
        y1 = x.duplicate()
        y2 = x.duplicate()
        y3 = x.duplicate()
        self.kspL.solve(-x, y1)
        self.Fp.mult(y1, y2)
        self.kspM.solve(y3, y)


n = int(2**6)
mesh = UnitSquareMesh(n, n)

V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)
W = V*Q

(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)

class r0(Expression):
    def __init__(self):
        self.M = 1
    def eval_cell(self, values, x, ufc_cell):
        if abs(1-x[1]) < 1e-6:
            values[0] = 1.0
            values[1] = 0.0
        else:
            values[0] = 0.0
            values[1] = 0.0
    def value_shape(self):
        return (2,)
f = r0()

u_k = Function(V)
p_k = Function(Q)
n = FacetNormal(mesh)
nu = 1.0
a11_laplacian = nu*inner(grad(v), grad(u))*dx(mesh)
a11_advection = inner(grad(u)*u_k, v)*dx(mesh) \
    + (1./2)*div(u_k)*inner(u, v)*dx(mesh) \
    - (1./2)*inner(u_k, n)*inner(u, v)*ds(mesh)
a12 = -div(v)*p*dx
a21 = -div(u)*q*dx
a = a11_laplacian+a11_advection+a12+a21

r11 = nu*inner(grad(v), grad(u_k))*dx(mesh) \
    + inner((grad(u_k)*u_k), v)*dx(mesh) \
    + (1./2)*div(u_k)*inner(u_k, v)*dx(mesh) \
    - (1./2)*inner(u_k, n)*inner(u_k, v)*ds(mesh)
r12 = -div(v)*p_k*dx
r21 = -div(u_k)*q*dx
RHSform = r11+r12+r21

l = inner(f, v)*dx
precNS = nu*inner(grad(v), grad(u))*dx(mesh) \
     + inner((grad(u)*u_k), v)*dx(mesh) \
     + (1./2)*div(u_k)*inner(u, v)*dx(mesh) \
     - (1./2)*inner(u_k, n)*inner(u, v)*ds(mesh) \
     + inner(p, q)*dx(mesh) \
     - div(v)*p*dx(mesh)
precS = nu*inner(grad(v), grad(u))*dx(mesh) \
      + (1./nu)*inner(p, q)*dx

def boundary(x, on_boundary):
    return on_boundary

eps = 1
tol = 1e-5
maxiter = 21
iter = 0
is0 = PETSc.IS().createGeneral(W.sub(0).dofmap().dofs())
is1 = PETSc.IS().createGeneral(W.sub(1).dofmap().dofs())

p = TrialFunction(Q)
q = TestFunction(Q)

fp = nu*inner(grad(q), grad(p))*dx(mesh) \
   + inner(inner(grad(p), u_k), q)*dx(mesh) \
   - (1./2)*inner(u_k, n)*inner(q, p)*ds(mesh) \
   - (1/2)*(dot(u_k('+'), n('+')) + dot(u_k('-'), n('-'))) \
   *avg(inner(q, p))*ds(mesh) \
   -dot(avg(q), dot(outer(p('+'), n('+')) + outer(p('-'),n('-')), \
    avg(u_k)))*dS(mesh)
Lp = as_backend_type(assemble(nu*inner(grad(q), grad(p))*dx(mesh))).mat()
Mp = as_backend_type(assemble((1./nu)*inner(q, p)*dx(mesh))).mat()

while eps > tol and iter < maxiter:
    A = PETScMatrix()
    P = PETScMatrix()
    b = PETScVector()
    b1 = PETScVector()
    if iter == 0:
        bc = DirichletBC(W.sub(0), f, boundary)
        assemble_system(a11_laplacian + a12 + a21, l, bc, A_tensor = A, b_tensor = b)
        assemble_system(precS, l, bc, A_tensor = P, b_tensor = b1)

        A_petsc = A.mat()
        b_petsc = b.vec()
        P_petsc = P.mat()

        u = b_petsc.duplicate()

        ksp = PETSc.KSP()
        ksp.create(comm=PETSc.COMM_WORLD)
        pc = ksp.getPC()
        ksp.setType('minres')
        ksp.setTolerances(1e-10)
        ksp.max_it = 100

        pc.setType(PETSc.PC.Type.FIELDSPLIT)

        pc.setFieldSplitIS(('u', is0), ('p', is1))
        pc.setFieldSplitType(0)
        ksp.setOperators(A_petsc, P_petsc)

        subksps = pc.getFieldSplitSubKSP()
        subksps[0].setType("preonly")
        subksps[0].getPC().setType("gamg")
        subksps[1].setType("preonly")
        subksps[1].getPC().setType("jacobi")

        reshist = {}
        def monitor(ksp, its, rnorm):
            reshist[its] = rnorm
        ksp.setMonitor(monitor)
        scale = b_petsc.norm()
        b_petsc = b_petsc/scale
        ksp.solve(b_petsc,u)
        u = scale*u
        b_petsc = b_petsc*scale
        r = b_petsc - A_petsc*u
    else:
        bc = DirichletBC(W.sub(0), Expression(("0.0", "0.0")), boundary)
        assemble_system(a, l-RHSform, bc, A_tensor = A, b_tensor = b)
        assemble_system(precNS, l, bc, A_tensor = P, b_tensor = b1)

        A_petsc = A.mat()
        b_petsc = b.vec()
        P_petsc = P.mat()

        ksp = PETSc.KSP()
        ksp.create(comm=PETSc.COMM_WORLD)
        pc = ksp.getPC()
        ksp.setType('gmres')
        ksp.setTolerances(1e-5)

        pc.setType(PETSc.PC.Type.FIELDSPLIT)
        ksp.setOperators(A_petsc,P_petsc)

        is0 = PETSc.IS().createGeneral(W.sub(0).dofmap().dofs())
        is1 = PETSc.IS().createGeneral(W.sub(1).dofmap().dofs())
        pc.setFieldSplitIS(('u', is0), ('p', is1))
        pc.setFieldSplitType(0) # 0=additive
        subksps = pc.getFieldSplitSubKSP()
        subksps[0].setType("preonly")
        subksps[0].getPC().setType("gamg")
        subksps[1].setType("preonly")
        PC = subksps[1].getPC()
        PC.setType("python")
        Fp = as_backend_type(assemble(fp)).mat()
        PC.setPythonContext(Approx(Fp, Mp, Lp))

        reshist = {}
        def monitor(ksp, its, rnorm):
            reshist[its] = rnorm
        ksp.setMonitor(monitor)
        scale = b_petsc.norm()
        b_petsc = b_petsc/scale
        ksp.solve(b_petsc,u)
        u = scale*u
        b_petsc = b_petsc*scale
        print 'its: ', ksp.its, 'norm:  ', (b_petsc-A_petsc*u).norm()
        print reshist


    u1 = u.getSubVector(is0)
    p1 = u.getSubVector(is1)
    print u1.norm()
    print p1.norm()
    eps = u.norm()
    # if iter != 0:
    print '\n\n\niter=%d: norm=%g' % (iter, eps)

    u2 = Function(V)
    u2.vector()[:] = u1.array + u_k.vector().array()
    p2 = Function(Q)
    p2.vector()[:] = p1.array + p_k.vector().array()
    # p2.vector()[:] += - assemble(p2*dx)/assemble(ones*dx)
    print np.linalg.norm(u2.vector().array())

    # sss
    u_k.assign(u2)
    p_k.assign(p2)
    iter += 1
    ksp.destroy(), A_petsc.destroy(), P_petsc.destroy()


# p = plot(p_k)




