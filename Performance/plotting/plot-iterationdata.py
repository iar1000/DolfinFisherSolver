import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


data = pd.DataFrame()
data = pd.read_csv("iterationdata.csv-0", sep=",")

ax00 = plt.subplot(311)
ax10 = plt.subplot(313)
ax11 = plt.subplot(312)

# dt
dt = data.iloc[:, 0]
dt.index = range(0, data.shape[0])
ax00.set_ylabel("dt")
ax00.set_xlabel("number of time steps")
ax00.plot(range(0, data.shape[0]), dt)

# residual
res = data.iloc[:,1]
res.index = range(0, data.shape[0])
ax10.set_ylabel("newton residual")
ax10.set_xlabel("number of time steps")
ax10.set_yscale("log")
ax10.plot(range(0, data.shape[0]), res)

# t
t = pd.Series(data=data.index, index=range(0, data.shape[0]))
ax11.set_ylabel("t")
ax11.set_xlabel("number of time steps")
#ax11.set_yscale("log")
ax11.plot(range(0, data.shape[0]), t)
plt.show()



