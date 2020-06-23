import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# 2D plot of the quantities of interesst extracted from simulation runs
# current quantities:
#   Volume of tumor
#   Total tumor concentration

# read in data
data = pd.DataFrame()
data = pd.read_csv("ParameterInference-results.csv", sep=",")

# format data
data.columns = ["t", "mesh", "D_white", "D_grey", "rho", "volume", "total concentration"]
data["t"] = data["t"].astype("float")
data["D_white"] = data["D_white"].astype("float")
data["D_grey"] = data["D_grey"].astype("float")
data["rho"] = data["rho"].astype("float")
data["volume"] = data["volume"].astype("float")
data["total concentration"] = data["total concentration"].astype("float")

# create dataframes for mesh cases
print("Find distinct meshes...")
mesh_names = data.mesh.unique()
mesh_frames = {}
for mesh_name in mesh_names:
    mesh_frames[mesh_name] = data[data["mesh"] == mesh_name]
    print("\t{} cases for mesh {}".format(len(mesh_frames[mesh_name]), mesh_name))
print()

# handle every mesh on it's own
for mesh_name in mesh_names:
    print("Format data for mesh {}".format(mesh_name))
    # find distinct cases
    curr_frame = mesh_frames[mesh_name]
    dist_cases = curr_frame.drop_duplicates(['D_white','rho'])[['D_white','rho']]
    dist_cases.index = list(range(0, len(dist_cases)))
    print("\tfound {} distinct D-rho combos".format(len(dist_cases)))
    # split data of distinct cases
    dist_case_frames = {}
    for i in list(range(0, len(dist_cases))):
        D = dist_cases.loc[i][0]
        rho = dist_cases.loc[i][1]
        print("\tseparate data for pair {} - {}...".format(D, rho), end=" ")
        dist_case_frames["{}-{}".format(D, rho)] = curr_frame[(curr_frame["D_white"] == D) & (curr_frame["rho"] == rho)]
        print("found {} timepoint(s)!".format(len(dist_case_frames["{}-{}".format(D, rho)])))

# plot qoi's
