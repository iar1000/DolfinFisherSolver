# plot brain data created by DolfinFisherSolver/Brain/compareTranslated()

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import pandas as pd

names_of_slice_data_folders_to_print = ["lh-white-hull-flood-0-1-merge-5-dof-100k-20-24-34"]
fig = plt.figure()
for matter in names_of_slice_data_folders_to_print:
    ims = []
    for s in range(0, 180):
        print(s)
        X = pd.read_csv("{}/slice-{}.txt".format(matter, s), sep=",")
        im = plt.imshow(X, cmap="gray", animated=True)
        ims.append([im])
    ani = animation.ArtistAnimation(fig, ims, interval=50, repeat=False)
    ani.save("{}-scan.mp4".format(matter))
