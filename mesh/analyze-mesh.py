

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
print("\tTotal nodes {}".format(len(dof_coords)))

subs = 1
subareas = [[46, 132, 67], [40, 132, 67], [52, 132, 67], [46, 138, 67], [46, 126, 67], [46, 132, 61], [46, 132, 73]]
subarea_dofs = []
# fill in dofs
for subarea in subareas:
    subarea_dofs.append([subarea, [x for x in dof_coords if sqrt((x[0] - subarea[0]) ** 2 + (x[1] - subarea[1]) ** 2 + (x[2] - subarea[2]) ** 2) < 3]])
# check
print("Subarea sizes : ")
for subdofs in subarea_dofs:
    print("\t", subdofs[0], len(subdofs[1]))

for subarea in subarea_dofs:
    worker_dofs = subarea[1]
    counter = 1
    max_counter = len(worker_dofs)
    # compare each dof to all others but itself (very slow) to find closest dofs
    for dof in worker_dofs:
        worker_dofs = [x for x in worker_dofs if (x[0] != dof[0] or x[1] != dof[1] or x[2] != dof[2])]
        counter = counter + 1
        worker = np.min(np.linalg.norm(sub - np.array(dof), axis=1))
        min_dist = min(min_dist, worker)
        print("{} @ ({}%)".format(min_dist, counter/max_counter * 100))
        print("\tcurrently round {}/{}, point {} @ {}".format(subs, len(subareas), dof, worker))
    subs = subs + 1

print("Mesh Statistics of {}".format(mesh_name))
print("\tDofs: ", max_counter)
print("\tMax vertex distance: ", mesh.hmax())
print("\tMin vertex distance: ", mesh.hmin())
print("\tMin node distance: ", min_dist)