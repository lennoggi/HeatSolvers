import numpy as np
import h5py
from matplotlib import pyplot as plt
from matplotlib.colors import SymLogNorm
import os

# ***** Parameters *****
filename = "../Data.h5"
outdir   = "Frames"

##umin = 0.
##umax = 15. 
clb_num_orders_of_mag = 7

cmap = "rainbow"

figsize = (10., 8.)
dpi     = 200

time_x =  0.65
time_y = -0.08
# **********************


##if os.path.isdir(outdir):
##    raise RuntimeError(f"Directory '{outdir}' already exists. Please remove it or change the 'outdir' variable before re-running this script.")
os.makedirs(outdir, exist_ok = False)

with h5py.File(filename, "r") as f:
    t         = f["Time"][()]
    x         = f["Coordinates/x"][()]
    y         = f["Coordinates/y"][()]
    out_every = f["Output frequency"][()]

    NT = len(t)
    NX = len(x)
    NY = len(y)

    assert NT > 0
    assert NX > 0
    assert NY > 0

    assert len(out_every) == 1
    out_every = out_every[0]

    assert x[-1] > x[0]
    assert y[-1] > y[0]
    extent = (x[0], x[-1], y[0], y[-1])


    print("Finding global min and max in the solution...")

    umin_glob =  np.inf
    umax_glob = -np.inf

    for n in range(NT):
        it = n*out_every
        u  = f[f"Solution/Iteration {it}"][()]
        assert u.shape == (NX, NY)

        umin_loc = np.min(u)
        umax_loc = np.max(u)

        ##if umin_loc < 0.0:
        ##    raise RuntimeError(f"The minimum of the solution at time {t[n]} is negative ({umin_loc}), which is not expected for the heat equation")

        if umin_loc < umin_glob: umin_glob = umin_loc
        if umax_loc > umax_glob: umax_glob = umax_loc

    print(f"Done\tumin_glob = {umin_glob}, umax_glob = {umax_glob}")

    if umax_glob < 0.0:
        raise RuntimeError("SymLogNorm does not support umax_glob = {umax_glob} < 0")

    clb_linthresh = umax_glob*10**(-clb_num_orders_of_mag)


    # Plot the solution
    for n in range(NT):
        it = n*out_every
        u  = f[f"Solution/Iteration {it}"][()]

        fig = plt.figure(figsize = figsize, dpi = dpi)
        ax  = fig.gca()

        ax.set_xlabel("$x$", fontsize = 12.)
        ax.set_ylabel("$y$", fontsize = 12.)

        im = ax.imshow(u,
            origin        = "lower",  # Transpose
            extent        = extent,
            interpolation = "none",
            aspect        = "equal",
            ##vmin          = umin_glob,
            ##vmax          = umax_glob,
            cmap          = cmap,
            norm          = SymLogNorm(linthresh = clb_linthresh, vmin = clb_linthresh, vmax = umax_glob)
        )

        clb = fig.colorbar(im, ax = ax)
        clb.set_label(r"$u\left(x, y\right)$", fontsize = 12.)

        ax.text(time_x, time_y, f"t = {t[n]:.2e}", transform = ax.transAxes, fontsize = 12., fontweight = "bold")
        ##plt.tight_layout()
        fig.savefig(f"{outdir}/it_{n:06d}.png")
        plt.close(fig)

        print(f"Iteration {it}/{NT*out_every}, time {t[n]:.2e} done")
