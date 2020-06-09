import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from os import listdir
from os.path import isfile, join
import re

# script for plotting type 1 generated data of Performance-FisherSolver
template_filename = "weakscaling-80k-type-1-tol-1e-08-nprocs-2-dofpr-39943.csv"

# collect data
output_files = [f for f in listdir("../output") if isfile(join("../output", f))]
for dofs in [80]:
    data = pd.DataFrame()

    reg = re.compile("weakscaling-{}k-type-1-tol-.+-nprocs-.+-dofpr-.+.csv".format(dofs))
    files = [fn for fn in output_files if bool(re.match(reg, fn))]
    for filename in files:
        try:
            current = pd.read_csv("../output/" + filename, sep=",",
                                  names= (["ls", "pc", "dw", "dg", "rho", "newton iterations", "krylov iterations", "time for solve()",
                                         "newton relative residual", "newton abs residual", "residuals"] + list(range(0, 70))))
            # drop header
            current = current.drop([0])
            # extract important data
            current.insert(0, "procs", int(filename.split("-")[8]))
            current.insert(0, "real-dofs", filename.split("-")[10].split(".")[0])
            data = data.append(current)
        except FileNotFoundError as e:
            print(e)

    # rearrange data
    data = data[data.ls != "ls"]
    print("shape of collected performance data: " + str(data.shape))
    data = data.sort_values(["ls", "pc", "procs"])
    data = data.drop_duplicates(subset=["ls", "pc", "procs"], keep='first')
    data = data.append(pd.DataFrame(np.zeros((1, data.shape[1]))))
    data.index = range(0, data.shape[0])

    # plot
    ls_to_sp = {"gmres": [0, 0], "cg": [1, 0], "minres": [0, 2],
                "bicgstab": [0, 1], "tfqmr": [1, 1], "richardson": [1, 2]}

    # plot time to solve
    fig, axes = plt.subplots(nrows=2, ncols=3)
    start = 0
    end = -1
    pre = [data.loc[0]["ls"], data.loc[0]["pc"]]
    for index, curr_series in data.iterrows():
        procs = curr_series["procs"]
        ls = curr_series["ls"]
        pc = curr_series["pc"]
        if not pre == [ls, pc]:
            end = index
            pre_old = pre
            pre = [ls, pc]
            frame = data.iloc[start:end, :]
            print("new subframe for " + str(pre_old[0]) + " " + str(pre_old[1]) + " " + str(
                dofs) + " of shape " + str(frame.shape))
            time = pd.Series(frame.loc[:, "time for solve()"])
            time.index = frame.loc[:, "procs"]
            time = pd.to_numeric(time)
            ax = time.plot(ax=axes[ls_to_sp[pre_old[0]][0], ls_to_sp[pre_old[0]][1]],
                           title=(str(pre_old[0]) + " at " + str(dofs) + "k dofs per rank"), kind="line",
                           label=(pre_old[1]),
                           legend=True)

            ax.set_xlabel("MPI Ranks")
            ax.set_ylabel("Wall time elapsed")

            start = index
    plt.show()

    # plot efficiency
    fig, axes = plt.subplots(nrows=2, ncols=3)
    start = 0
    end = -1
    pre = [data.loc[0]["ls"], data.loc[0]["pc"]]
    for index, curr_series in data.iterrows():
        procs = curr_series["procs"]
        ls = curr_series["ls"]
        pc = curr_series["pc"]
        if not pre == [ls, pc]:
            end = index
            pre_old = pre
            pre = [ls, pc]
            frame = data.iloc[start:end, :]
            print("new subframe for " + str(pre_old[0]) + " " + str(pre_old[1]) + " " + str(
                dofs) + " of shape " + str(frame.shape))
            time = pd.Series(frame.loc[:, "time for solve()"])
            time.index = frame.loc[:, "procs"]
            time = pd.to_numeric(time)
            onetime = pd.Series(data=time.iloc[0], index=time.index)
            time = onetime / time
            ax = time.plot(ax=axes[ls_to_sp[pre_old[0]][0], ls_to_sp[pre_old[0]][1]],
                           title=(str(pre_old[0])),
                           kind="line",
                           label=(pre_old[1]),
                           legend=True, fontsize=14)

            ax.set_xlabel("MPI Ranks", fontsize=14)
            ax.set_ylabel("Efficiency", fontsize=14)
            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)
            # plt.legend(prop={"size":15})

            start = index
    plt.show()

    plt.rc('pgf', texsystem='pdflatex')
    nice_fonts = {
        # Use LaTeX to write all text
        "text.usetex": True,
        "font.family": "serif",
        # Use 10pt font in plots, to match 10pt font in document
        "axes.labelsize": 12,
        "font.size": 11,
        # Make the legend/label fonts a little smaller
        "legend.fontsize": 12,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
    }

    mpl.rcParams.update(nice_fonts)

    # plot efficiency of cg
    pretocol = {"hypre_euclid": [0.95, 0.52, 0.1, 1], "petsc_amg": [0.25, 0.75, 0.25, 1],
                "hypre_amg": [0.04, 0.45, 1, 1]}
    start = 0
    end = -1
    pre = [data.loc[0]["ls"], data.loc[0]["pc"]]
    # plot ideal scaling
    plt.hlines(1, 1, 480, label="ideal", linestyles="dashed")
    for index, curr_series in data.iterrows():
        procs = curr_series["procs"]
        ls = curr_series["ls"]
        pc = curr_series["pc"]
        if not pre == [ls, pc]:
            end = index
            pre_old = pre
            pre = [ls, pc]
            frame = data.iloc[start:end, :]
            print("new subframe for " + str(pre_old[0]) + " " + str(pre_old[1]) + " " + str(
                dofs) + " of shape " + str(frame.shape))
            time = pd.Series(frame.loc[:, "time for solve()"])
            time.index = frame.loc[:, "procs"]
            time = pd.to_numeric(time)
            onetime = pd.Series(data=time.iloc[0], index=time.index)
            time = onetime / time
            if(pre_old[0] == "cg" and pre_old[1] != "jacobi"):
                ax = time.plot(kind="line", style="x-", label=(str(pre_old[1]).replace("_", " ")), color=pretocol[pre_old[1]])

            ax.set_xlabel("MPI Ranks", fontsize=11)
            ax.set_ylabel("Efficiency", fontsize=11)
            plt.xticks(fontsize=11)
            plt.yticks(fontsize=11)
            plt.legend(prop={"size":12})

            start = index
    plt.savefig("weakscaling-cg-80k.pgf")
    plt.show()
