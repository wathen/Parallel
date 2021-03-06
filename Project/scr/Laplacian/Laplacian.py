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


V = VectorFunctionSpace(mesh, "CG", 1)
comm = mpi_comm_world()
mpiRank = MPI.rank(comm)

# p = plot(mesh)
# p.write_png('mesh'+str(mpiRank))
# file = File("Mesh/mesh"+str(mpiRank)+".xdmf")
# file << mesh
meshname = "Mesh/mesh"+str(mpiRank)
# meshfunction = MeshFunction("size_t", mesh, meshname + "_meshfunction.xml")
# f = HDF5File(meshfunction.mpi_comm(), meshname+"_meshfunction.hdf5", 'w')
# f.write(meshfunction, meshname_meshfunction)
# f = HDF5File("Mesh/mesh"+str(mpiRank)+'.h5', 'w')
# f.write(mesh, 'mesh')
f = HDF5File(mesh.mpi_comm(), meshname+".hdf5", 'w')
f.write(mesh, meshname)

print V.dofmap().dofs()
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
dump_timings_to_xml('test.xml', True)
