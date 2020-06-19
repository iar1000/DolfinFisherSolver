

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
unchecked_coords = V.tabulate_dof_coordinates()
min_dist = 99999
max_dist = -99999
counter = 1
max_counter = len(dof_coords)
print("\tTotal nodes {}".format(len(dof_coords)))
# remove all nodes which are not in the high resolution part of the mesh
unchecked_coords = [x for x in unchecked_coords if sqrt((x[0]-46)**2 + (x[1]-132)**2 + (x[2]-67)**2) < 9]
print("\tCompare {} nodes...".format(len(unchecked_coords)))
# compare each dof to all others but itself (very slow) to find closest dofs
for dof in unchecked_coords:
    unchecked_coords = [x for x in unchecked_coords if (x[0] != dof[0] or x[1] != dof[1] or x[2] != dof[2])]
    counter = counter + 1
    worker = np.min(np.linalg.norm(unchecked_coords - np.array(dof), axis=1))
    min_dist = min(min_dist, worker)
    print("{} ({}%)".format(min_dist, counter/max_counter * 100))
    print("\t{} for pairs of {}".format(dof, worker))

print("Mesh Statistics of {}".format(mesh_name))
print("\tDofs: ", max_counter)
print("\tMax vertex distance: ", mesh.hmax())
print("\tMin vertex distance: ", mesh.hmin())
print("\tMin node distance: ", min_dist)