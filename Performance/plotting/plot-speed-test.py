import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import pandas as pd
import numpy as np
from os import listdir
from os.path import isfile, join
import re

# script for plotting type 1 generated data of Performance-FisherSolver
template_filename = "speed-type-1-tol-1e-08-nprocs-1-dofpr-600905.csv"

# collect data
output_files = [f for f in listdir("../output") if isfile(join("../output", f))]
data = pd.DataFrame()
print(output_files)
for procs in [36, 48, 96, 120, 180, 240, 360, 480, 540]:
    reg = re.compile("speed-type-1-tol-.+-nprocs-{}-dofpr-.+.csv".format(procs))
    filename = [fn for fn in output_files if bool(re.match(reg, fn))]
    if filename:
        try:
            current = pd.read_csv("../output/" + filename[0], sep=",",
                                  names=(["ls", "pc", "dw", "dg", "rho", "newton iterations", "krylov iterations",
                                          "time for solve()",
                                          "newton relative residual", "newton abs residual", "residuals"] + list(
                                      range(0, 70))))
            # drop header
            current = current.drop([0])
            current.insert(0, "procs", procs)
            data = data.append(current)
        except FileNotFoundError as e:
            print(e)

# rearrange data
data = data[data.ls != "ls"]
data = data.sort_values(["ls", "pc", "procs"])
data = data.append(pd.DataFrame(np.zeros((1, data.shape[1]))))
data.index = range(0, data.shape[0])
print("shape of collected performance data: " + str(data.shape))

# plot
ls_to_sp = {"gmres": [0, 0], "cg": [1, 0], "minres": [0, 2],
            "bicgstab": [0, 1], "tfqmr": [1, 1], "richardson": [1, 2]}

# plot time to solve
fig, axes = plt.subplots(nrows=2, ncols=3)
start = 0
end = -1
pre = [data.iloc[0]["ls"], data.iloc[0]["pc"]]
for index, curr_series in data.iterrows():
    procs = curr_series["procs"]
    ls = curr_series["ls"]
    pc = curr_series["pc"]
    if not pre == [ls, pc]:
        end = index
        pre_old = pre
        pre = [ls, pc]
        frame = data.iloc[start:end, :]
        print("new subframe for " + str(pre_old[0]) + " " + str(pre_old[1]) + " of shape " + str(frame.shape))
        time = pd.Series(frame.loc[:, "time for solve()"])
        time.index = frame.loc[:, "procs"]
        time = pd.to_numeric(time)
        print(time)
        ax = time.plot(ax=axes[ls_to_sp[pre_old[0]][0], ls_to_sp[pre_old[0]][1]], title=pre_old[0], kind="line",
                       label=(pre_old[1]), legend=True)
        ax.set_xlabel("mpi ranks")
        ax.set_ylabel("time elapsed in seconds for solving F(u)=0")
        #ax.set_yscale("log")

        start = index
plt.show()
