# script to create simple h5 meshes

from dolfin import *

# generate w on b rectangle mesh with spatial resolutio of res
# @output: save as h5/xml file in ../mesh/rect-(w)on(b)-res-(res).h5
def rect_mesh(w, b, res):
    mesh = RectangleMesh(Point(0, 0), Point(w, b), res, res, diagonal="right")
    # create h5 mesh
    path = "../mesh/rect-" + str(w) + "on" + str(b) + "-res-" + str(res) + ".h5"
    file = HDF5File(MPI.comm_world, path, 'w')
    file.write(mesh, "/mesh")

# generate w on l on d box mesh with spatial resolution in each direction of res
# @output: save as h5/xml file in ../mesh/box-(w)on(l)on(d)-res-(res).h5
def box_mesh(w, l, d, res):
    mesh = BoxMesh(Point(0, 0, 0), Point(w, l, d), res, res, res)
    path = "box-" + str(w) + "on" + str(l) + "on" + str(d) + "-res-" + str(res) + ".h5"
    file = HDF5File(MPI.comm_world, path, 'w')
    file.write(mesh, "/mesh")
  

# generate standart test meshes
rect_mesh(10, 10, 35)     # 2.5k elements, 5k dof
rect_mesh(10, 10, 100)    # 20k elements, 40k dof
rect_mesh(10, 10, 340)    # 240k elements, 480k dof

# n = pow((1.0 * totalDofs / 1.30) / 6, 1.0/3.0)
box_mesh(10, 10, 10, 23)  # 100k dof
box_mesh(10, 10, 10, 29)  # 200k dof
box_mesh(10, 10, 10, 40)  # 500k dof
box_mesh(10, 10, 10, 50)  # 1 Mio dof
box_mesh(10, 10, 10, 108) # 10 Mio dof
box_mesh(10, 10, 10, 136) # 20 Mio dof
box_mesh(10, 10, 10, 156) # 30 Mio dof
box_mesh(10, 10, 10, 172) # 40 Mio dof
box_mesh(10, 10, 10, 185) # 50 Mio dof

#box_mesh(10,10,10,67)	  # 1.8mio elements, 2.4mio dof
#box_mesh(10,10,10,115)	  # 12'000'000 DOFs


