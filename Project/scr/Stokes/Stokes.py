import petsc4py
import sys
petsc4py.init(sys.argv)
from petsc4py import PETSc
from dolfin import *

import MatrixOperations as MO
import numpy as np

parameters['linear_algebra_backend'] = 'PETSc'
parameters["form_compiler"]["optimize"]     = True
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "quadrature"
parameters['reorder_dofs_serial'] = False


# class Approx(object):

#     def __init__(self, W, A):
#         self.W = W
#         self.A = A
#         IS = MO.IndexSet(W)
#         self.u_is = IS[0]
#         self.p_is = IS[1]
#     def create(self, pc):
#         kspL = PETSc.KSP()
#         kspL.create(comm=PETSc.COMM_WORLD)
#         pcL = kspL.getPC()
#         kspL.setType('cg')
#         pcL.setType('jacobi')
#         kspL.setFromOptions()
#         self.kspL = kspL

#         kspM = PETSc.KSP()
#         kspM.create(comm=PETSc.COMM_WORLD)
#         pcM = kspM.getPC()
#         kspM.setType('cg')
#         pcM.setType('jacobi')
#         kspM.setFromOptions()
#         self.kspM = kspM


#     def setUp(self, pc):
#         A, P = pc.getOperators()
#         L = A.getSubMatrix(self.u_is,self.u_is)
#         M = P.getSubMatrix(self.p_is,self.p_is)
#         self.kspM.setOperators(M,M)
#         self.kspL.setOperators(L,L)


#     def apply(self, pc, x, y):
#         x1 = x.getSubVector(self.u_is)
#         y1 = x1.duplicate()
#         x2 = x.getSubVector(self.p_is)
#         y2 = x2.duplicate()
#         self.kspL.solve(x1, y1)
#         self.kspM.solve(x2, y2)

#         y.array = (np.concatenate([y1.array, y2.array]))


n = int(2**5)
mesh = UnitSquareMesh(n,n)
V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)
W = V*Q
F = Function(Q)
print F.vector().size()
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

nu = 1.0
a11 = nu*inner(grad(v), grad(u))*dx
a12 = -div(v)*p*dx
a21 = -div(u)*q*dx
l  = inner(v, f)*dx
a = a11+a12+a21
p = nu*inner(grad(v), grad(u))*dx + (1./nu)*inner(q,p)*dx

def boundary(x, on_boundary):
    return on_boundary

A = PETScMatrix()
P = PETScMatrix()
b = PETScVector()
b1 = PETScVector()
bc = DirichletBC(W.sub(0), f, boundary)

assemble_system(a, l, bc, A_tensor = A, b_tensor = b)
assemble_system(p, l, bc, A_tensor = P, b_tensor = b1)

A_petsc = A.mat()
b_petsc = b.vec()
P_petsc = P.mat()
u = b_petsc.duplicate()

ksp = PETSc.KSP()
ksp.create(comm=PETSc.COMM_WORLD)
pc = ksp.getPC()
ksp.setType('minres')
# pc.setType('python')
# pc.setPythonContext(Approx(W,P_petsc))
pc.setType(PETSc.PC.Type.FIELDSPLIT)
is0 = PETSc.IS().createGeneral(W.sub(0).dofmap().dofs())
is1 = PETSc.IS().createGeneral(W.sub(1).dofmap().dofs())
pc.setFieldSplitIS(('u', is0), ('p', is1))
pc.setFieldSplitType(0) # 0=additive
subksps = pc.getFieldSplitSubKSP()
subksps[0].setType("preonly")
subksps[0].getPC().setType("gamg")
subksps[1].setType("preonly")
subksps[1].getPC().setType("jacobi")
print subksps[0].view()
print subksps[1].view()
reshist = {}
def monitor(ksp, its, rnorm):
    reshist[its] = rnorm
ksp.setMonitor(monitor)
ksp.setOperators(A_petsc,P_petsc)

scale = b_petsc.norm()
b_petsc = b_petsc/scale
ksp.solve(b_petsc,u)
u = scale*u
r = b_petsc*scale-A_petsc*u
print ksp.its, r.norm()
print reshist






