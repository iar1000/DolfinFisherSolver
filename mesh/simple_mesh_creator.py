# script to create simple h5 meshes

from dolfin import *

# generate w on b rectangle mesh with spatial resolutio of res
# @output: save as h5 file in ../mesh/rect-(w)on(b)-res-(res).h5
def rect_mesh(w, b, res):
    mesh = RectangleMesh(Point(0, 0), Point(w, b), res, res, diagonal="right")
    path = "../mesh/rect-" + str(w) + "on" + str(b) + "-res-" + str(res) + ".h5"
    file = HDF5File(MPI.comm_world, path, 'w')
    file.write(mesh, "/mesh")

# generate w on l on d box mesh with spatial resolution in each direction of res
# @output: save as h5 file in ../mesh/box-(w)on(l)on(d)-res-(res).h5
def box_mesh(w, l, d, res):
    mesh = BoxMesh(Point(0, 0, 0), Point(w, l, d), res, res, res)
    path = "box-" + str(w) + "on" + str(l) + "on" + str(d) + "-res-" + str(res) + ".h5"
    file = HDF5File(MPI.comm_world, path, 'w')
    file.write(mesh, "/mesh")

# generate standart test meshes
rect_mesh(10, 10, 25)     # 1.25k elements, 2.5k dof
rect_mesh(10, 10, 35)     # 2.5k elements, 5k dof
rect_mesh(10, 10, 50)     # 5k elements, 10k dof
rect_mesh(10, 10, 70)     # 10k elements, 20k dof
rect_mesh(10, 10, 85)     # 15k elements, 30k dof
rect_mesh(10, 10, 100)    # 20k elements, 40k dof
rect_mesh(10, 10, 120)    # 30k elements, 60k dof
rect_mesh(10, 10, 140)    # 40k elements, 80k dof
rect_mesh(10, 10, 170)    # 60k elements, 120k dof
rect_mesh(10, 10, 244)    # 120k elements, 240k dof
rect_mesh(10, 10, 300)    # 180k elements, 360k dof
rect_mesh(10, 10, 340)    # 240k elements, 480k dof

box_mesh(10,10,10,15)       # 20'000 elements