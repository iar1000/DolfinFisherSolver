

import argparse
from dolfin import *
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument("--mesh", type=str, help="name of the mesh to analyze", default="lh-plial-ball-9mio.xdmf")
args = parser.parse_args()

mesh_name = args.mesh

# load mesh
try:
    mesh = Mesh()
    f = XDMFFile(mesh.mpi_comm(), '{}.xdmf'.format(mesh_name.split(".")[0]))
    f.read(mesh)
    f.close()
except Exception as e:
    print("Error loading mesh!")
    print(e)
    quit()

print("Analyze mesh...")
V = FunctionSpace(mesh, "Lagrange", 2)

dof_coords = V.tabulate_dof_coordinates()
min_dist = 99999
max_dist = -99999
counter = 1
max_counter = len(dof_coords)
for dof in dof_coords:
    counter = counter + 1
    norms = np.linalg.norm(dof_coords - np.array(dof), axis=1)
    norms = [x for x in dof_coords if x != 0]
    min_dist = min(min_dist, np.min(norms))
    print("{} ({}%)".format(min_dist, counter/max_counter * 100))
    print("\t{} for pairs of {}".format(dof, np.min(worker_dist)))

print("Mesh Statistics of {}".format(mesh_name))
print("\tDofs: ", max_counter)
print("\tMax vertex distance: ", mesh.hmax())
print("\tMin vertex distance: ", mesh.hmin())