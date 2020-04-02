from dolfin import *
import numpy

# script to plot result of mesh convergence study

# Reference solution 30
print("read 30")
mesh_file_30 = HDF5File(MPI.comm_world, "../mesh/box-10on10on10-res-156.h5", "r")
mesh30 = Mesh()
mesh_file_30.read(mesh30, "mesh", False)
V30 = FunctionSpace(mesh30, "Lagrange", 2)
u30 = Function(V)
times30 = ["0", "1.0001", "2.00205", "2.99966", "3.99353", "4.99856", "5.99878", "6.99723", "8.00391"]

print("read 10")
mesh_file_10 = HDF5File(MPI.comm_world, "../mesh/box-10on10on10-res-108.h5", "r")
mesh10 = Mesh()
mesh_file_10.read(mesh10, "mesh", False)
V10 = FunctionSpace(mesh10, "Lagrange", 2)
u10 = Function(V)
times10 = ["0", "1.00101", "2.00353", "3.00161", "3.99587", "5.00106", "6.0012", "6.99926", "8.00531"]

print("read 1")
mesh_file_1 = HDF5File(MPI.comm_world, "../mesh/box-10on10on10-res-50.h5", "r")
mesh1 = Mesh()
mesh_file_1.read(mesh1, "mesh", False)
V1 = FunctionSpace(mesh1, "Lagrange", 2)
u1 = Function(V)
times1 = ["0", "1.00055", "2.00204", "3.00102", "3.99718", "5.00022", "6.00874", "6.9997", "8.00835"]



# Other solutions
for t in range(0, 1):
    print("time ", t)
    # read in 30 solution of this timestep
    f30 = HDF5File(mesh30.mpi_comm(), "dofs-30664297-u-at-{}.h5".format(times30[t]), "r")
    f30.read(u30, "/u")
    coords30 = mesh30.coordinates()

    # read in 10 solution of this timestep
    f10 = HDF5File(mesh10.mpi_comm(), "dofs-10218313-u-at-{}.h5".format(times10[t]), "r")
    f10.read(u10, "/u")
    coords10 = mesh10.coordinates()

    # read in 10 solution of this timestep
    f1 = HDF5File(mesh1.mpi_comm(), "dofs-1030301-u-at-{}.h5".format(times1[t]), "r")
    f1.read(u1, "/u")
    coords1 = mesh1.coordinates()

    # compare 30 to 30
    value10 = numpy.zeros(len(coords1), dtype=numpy.float_)
    value11 = numpy.zeros(len(coords1), dtype=numpy.float_)
    # compare 10 to 30
    value20 = numpy.zeros(len(coords1), dtype=numpy.float_)
    value21 = numpy.zeros(len(coords1), dtype=numpy.float_)
    # compare 1 to 30
    value30 = numpy.zeros(len(coords1), dtype=numpy.float_)
    value31 = numpy.zeros(len(coords1), dtype=numpy.float_)

    for v in vertices(mesh1):
        vertex_coordinate = coords1[v.index()]

        # evaluate at the same point and store difference
        value10[v.index()] = u30(vertex_coordinate)
        value11[v.index()] = u30(vertex_coordinate)
        value20[v.index()] = u30(vertex_coordinate)
        value21[v.index()] = u10(vertex_coordinate)
        value30[v.index()] = u30(vertex_coordinate)
        value31[v.index()] = u1(vertex_coordinate)

    print("norm: ", numpy.linalg.norm(value11 - value10))
    print("norm: ", numpy.linalg.norm(value20 - value21))
    print("norm: ", numpy.linalg.norm(value30 - value31))


# output
# out = File("u{}.pvd".format(res))
# out << u
