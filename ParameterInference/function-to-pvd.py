# script to convert function files to pvd files
# the more processors are being used with this script, the more memory is being consumed, since pvd stores a file per function per processor

import argparse
from fenics import *
import os

parser = argparse.ArgumentParser()

parser.add_argument("--meshpath", type=str, help="path to mesh folder", default="../mesh/")
parser.add_argument("--mesh", type=str, help="name of the mesh", default="lh-plial-3mio.xdmf")
parser.add_argument("--functionpath", type=str, help="path to function folder", default="")
parser.add_argument("--functioncount", type=int,
                    help="number of functions from 0 to functionscount that have to be converted", default=60)
args = parser.parse_args()

meshpath = args.meshpath
meshname = args.mesh
functioncount = args.functioncount
functionpath = args.functionpath

if not functionpath:
    print("ERROR: no path to function folder provided")
    quit()

print("read in mesh...", end=" ")
mesh = Mesh(MPI.comm_world)
with XDMFFile(MPI.comm_world, meshpath + meshname) as meshfile:
    meshfile.read(mesh)
print("success!")

print("create function space...", end=" ", flush=True)
V = FunctionSpace(mesh, "Lagrange", 2)
u = Function(V)
print("success!")

print("convert {} functions to pvd...".format(functioncount), end=" ", flush=True)
outpath = functionpath + "../convertedpvd/"
try:
    os.mkdir(outpath)
except:
    pass

outfile = File(outpath + "out.pvd")
for i in range(0, functioncount + 1):
    print("\tconvert function {}...".format(i), end=" ", flush=True)
    try:
        with HDF5File(MPI.comm_world, functionpath + "function{}.h5".format(i), "r") as infile:
            infile.read(u, "u")
            outfile << (u, i)
    except:
        print("\time point ", i, " doesn't exist")
    print("success!")
