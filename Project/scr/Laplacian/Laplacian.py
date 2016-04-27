from dolfin import *
import petsc4py
import sys
petsc4py.init(sys.argv)
from petsc4py import PETSc

parameters['linear_algebra_backend'] = 'PETSc'
parameters["form_compiler"]["optimize"]     = True
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "quadrature"
parameters['reorder_dofs_serial'] = False

N = int(sys.argv[1])
n = int(2**N)
mesh = UnitCubeMesh(n,n,n)


V = VectorFunctionSpace(mesh, "CG", 2)
v = TestFunction(V)
u = TrialFunction(V)
print V.dim()
class r0(Expression):
    def __init__(self):
        self.M = 1
    def eval_cell(self, values, x, ufc_cell):
        values[0] = 1.
        values[1] = 1.
        values[2] = 1.
    def value_shape(self):
        return (3,)
f = r0()
a = inner(grad(u),grad(v))*dx
l = inner(f,v)*dx
def boundary(x, on_boundary):
    return on_boundary

A = PETScMatrix()
b = PETScVector()
bc = DirichletBC(V,f, boundary)
t0 = Timer("assemble_system")
assemble_system(a, l, bc, A_tensor = A, b_tensor = b)
del(t0)

A_petsc = A.mat()
b_petsc = b.vec()
u = b_petsc.duplicate()
ksp = PETSc.KSP()
ksp.create(comm=PETSc.COMM_WORLD)
pc = ksp.getPC()
ksp.setType('cg')
pc.setType('gamg')
ksp.setFromOptions()
scale = b_petsc.norm()
b_petsc = b_petsc/scale
ksp.setOperators(A_petsc,A_petsc)
t0 = Timer("ksp.solve")
ksp.solve(b_petsc,u)
del(t0)
list_timings()
