import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.ticker as ticker

# 2D plot of the quantities of interesst extracted from simulation runs
# current quantities:
#   Volume of tumor
#   Total tumor concentration

# setup for thesis import
plt.rc('pgf', texsystem='pdflatex')
plt.rc('text', usetex=True)
nice_fonts = {
    # Use LaTeX to write all text
    "text.usetex": True,
    "font.family": "serif",
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 12,
    "font.size": 12,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
}
mpl.rcParams.update(nice_fonts)

# read in data
col_names = ["t", "mesh", "D_white", "rho", "volume", "total concentration"]
data = pd.DataFrame(columns = col_names)

filenames = ["paraview-volume-concentration-lh-plial-3mio-0-1-0-25-0-001",
             "paraview-volume-concentration-lh-plial-3mio-0-1-0-025-0-001"]
for file in filenames:
    temp = pd.read_csv(file + ".csv", sep=",", header=None, names=col_names)
    data = data.append(temp)


# format data
data["t"] = data["t"].astype("float")
data["D_white"] = data["D_white"].astype("float")
data["rho"] = data["rho"].astype("float")
data["volume"] = data["volume"].astype("float")
data["total concentration"] = data["total concentration"].astype("float")

# create dataframes for mesh cases
print("Find distinct meshes...")
mesh_names = data.mesh.unique()
mesh_frames = {}
for mesh_name in mesh_names:
    mesh_frames[mesh_name] = data[data["mesh"] == mesh_name]
    print("\t{} time-points for mesh {}".format(len(mesh_frames[mesh_name]), mesh_name))
print()

# handle every mesh on it's own
mesh_dist_cases = {}  # holds name of case : frame of case data
for mesh_name in mesh_names:
    print("Format data for mesh {}".format(mesh_name))
    # find distinct cases
    curr_frame = mesh_frames[mesh_name]
    dist_cases = curr_frame.drop_duplicates(['D_white', 'rho'])[['D_white', 'rho']]
    dist_cases.index = list(range(0, len(dist_cases)))
    print("\tfound {} distinct D-rho combos".format(len(dist_cases)))
    # split data of distinct cases
    dist_case_frames = {}
    for i in list(range(0, len(dist_cases))):
        D = dist_cases.loc[i][0]
        rho = dist_cases.loc[i][1]
        print("\tseparate data for pair {} - {}...".format(D, rho), end=" ")
        case_frame = curr_frame[(curr_frame["D_white"] == D) & (curr_frame["rho"] == rho)]
        case_frame["volume"] = case_frame["volume"] / case_frame["volume"].iloc[0]
        case_frame["total concentration"] = case_frame["total concentration"] / case_frame["total concentration"].iloc[0]
        dist_case_frames["{}-{}".format(D, rho)] = case_frame
        dist_case_frames["{}-{}".format(D, rho)] =case_frame
        dist_case_frames["{}-{}".format(D, rho)].index = list(
            range(1, 1 + len(dist_case_frames["{}-{}".format(D, rho)])))
        print("found {} timepoint(s)!".format(len(dist_case_frames["{}-{}".format(D, rho)])))

    mesh_dist_cases[mesh_name] = dist_case_frames
print()

# setup colors to use
thesis_case_to_color = {"0.1-0.025": "green", "0.1-0.25": "royalblue", "0.6-0.025": "orange",
                        "0.6-0.25": "crimson"}
thesis_mesh_to_alpha = {"lh-plial-3mio": 1}

# plot volume over time
fig, ax = plt.subplots()
ax.set_xlabel("t [$days$]", fontsize=12)
ax.set_ylabel("Volume growth factor", fontsize=12)
for mesh_name in mesh_dist_cases:
    print("Plot Volume data for mesh ", mesh_name)
    cases = mesh_dist_cases[mesh_name]
    for case in cases:
        print("\tPlot case", case)
        D = case.split("-")[0]
        rho = case.split("-")[1]
        cases[case]["volume"].plot(title="Tumor growth over time",
                                   label="D\_white={}, $\\rho$={}, on {}".format(D, rho, mesh_name), style="x-",
                                   color=list(
                                       colors.to_rgba(thesis_case_to_color[case], alpha=thesis_mesh_to_alpha[mesh_name])),
                                   legend=True)
plt.show()

# plot total concentration over time
fig, ax = plt.subplots()
ax.set_xlabel("t [$days$]", fontsize=12)
ax.set_ylabel("Summed up concentration growth factor", fontsize=12)
for mesh_name in mesh_dist_cases:
    print("Plot total concentration data for mesh ", mesh_name)
    cases = mesh_dist_cases[mesh_name]
    for case in cases:
        print("\tPlot case", case)
        D = case.split("-")[0]
        rho = case.split("-")[1]
        cases[case]["total concentration"].plot(title="Summed up tumor concentration over time",
                                   label="D\_white={}, $\\rho$={}, on {}".format(D, rho, mesh_name), style="x-",
                                   color=list(
                                       colors.to_rgba(thesis_case_to_color[case], alpha=thesis_mesh_to_alpha[mesh_name])),
                                   legend=True)
plt.show()
