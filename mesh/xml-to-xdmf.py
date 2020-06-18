
import sys
from dolfin import *

mesh_name = sys.argv[1]

try:
    mesh = Mesh(mesh_name)
    f = XDMFFile(mesh.mpi_comm(), '{}.xdmf'.format(mesh_name.split(".")[0]))
    f.write(mesh)
    f.close()

except Exception as e:
    print(e)