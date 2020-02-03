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

rect_mesh(100, 100, 200)