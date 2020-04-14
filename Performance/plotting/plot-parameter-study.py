import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import pandas as pd
import numpy as np
from os import listdir
from os.path import isfile, join
import re

# script for plotting type 1 generated data of Performance-FisherSolver
filename = "parameter-type-1-tol-1e-08-nprocs-96-dofpr-6310.csv"

# collect data
data = pd.read_csv("../output/" + filename, sep=",",
                   names=(["ls", "pc", "dw", "dg", "rho", "newton iterations", "krylov iterations", "time for solve()",
                           "newton relative residual", "newton abs residual"] + list(range(0, 70))))
# rearrange data
data = data[data.ls != "ls"]
data = data[data.ls != "richardson"]
print("shape of collected performance data: " + str(data.shape))
data = data.sort_values(["ls", "pc", "dw", "dg", "rho"])
data.index = range(0, data.shape[0])

# plot runtime
specifier = pd.DataFrame(data=(data["ls"].astype(str) + "-" + data["pc"].astype(str) + "-dw-" + data["dw"].astype(str) + "-rho-" + data[
    "rho"].astype(str)), columns=["name"])
specifier["runtime"] = data["time for solve()"]
specifier["runtime"] = pd.to_numeric(specifier["runtime"])
specifier.index = specifier["name"]
specifier.iloc[0:100].plot(kind="bar")
specifier.iloc[101:200].plot(kind="bar")
specifier.iloc[201:300].plot(kind="bar")
specifier.iloc[301:400].plot(kind="bar")
specifier.iloc[401:500].plot(kind="bar")
specifier.iloc[501:600].plot(kind="bar")
specifier.iloc[601:697].plot(kind="bar")

plt.show()

# plot
ls_to_naxis = {"gmres": 4, "cg": 4, "minres": 2,
               "bicgstab": 2, "tfqmr": 2}
pc_to_axis = {"hypre_amg": 0, "hypre_euclid": 1, "petsc_amg": 2, "jacobi": 3}

ls_old = "start"
for index, curr_series in data.iterrows():

    ls = curr_series["ls"]
    pc = curr_series["pc"]
    dw = float(curr_series["dw"])
    dg = float(curr_series["dg"])
    rho = float(curr_series["rho"])
    time = float(curr_series["time for solve()"])
    krylov_iters = int(curr_series["krylov iterations"])
    residuals = curr_series.iloc[10:krylov_iters + 10]
    residuals = pd.to_numeric(residuals)
    residuals.index = list(range(1, krylov_iters + 1))
    print("{} {} with parameters {} {} {} has {} krylov iterations: {}".format(
        ls, pc, dw, dg, rho, krylov_iters, residuals.tolist()
    ))

    if ls != ls_old:
        plt.show()
        fig, axes = plt.subplots(nrows=1, ncols=ls_to_naxis[ls])

    if (dw == 0.1 and rho == 0.025) or (dw == 0.1 and rho == 0.15) or (dw == 0.6 and rho == 0.025) or (
            dw == 0.6 and rho == 0.15) or (dw == 1e-06 and rho == 0.15) or (dw == 1e-06 and rho == 0.025) or (
            dw == 1e-06 and rho == 1e-06) or (dw == 0.1 and rho == 1e-06) or (dw == 0.6 and rho == 1e-06):
        ax = residuals.plot(ax=axes[pc_to_axis[pc]], title=("{} {} convergence".format(ls, pc)), kind="line",
                            label=("dw= {} dg= {} rho= {}".format(dw, dg, rho)), legend=True)
    ax.set_yscale("log")
    ax.set_xlabel("krylov iteration")
    ax.set_ylabel("krylov residual")

    if ls != ls_old:
        ls_old = ls

plt.show()


