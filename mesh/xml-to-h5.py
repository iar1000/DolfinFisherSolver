
import sys
from dolfin import *

mesh_name = sys.argv[1]

try:
    mesh = Mesh(mesh_name)
    f = HDF5File(mesh.mpi_comm(), '{}.h5'.format(mesh_name.split(".")[0]), 'w')
    f.write(mesh, "mesh/")

except Exception as e:
    print(e)

