# Script that can be used to plot the Voronoi grid contained in vtest.txt

import numpy as np
import matplotlib
from pathlib import Path

matplotlib.use("Agg")
import matplotlib.pyplot as pl
import argparse

argparser = argparse.ArgumentParser()
argparser.add_argument("--file", "-f", action="store", required=True)
argparser.add_argument("--zoom", "-z", action="store_true")
args = argparser.parse_args()

gs = np.fromregex(
    args.file, "G\s+(\S+)\s+(\S+)", [("x", np.float64), ("y", np.float64)]
)
ms = np.fromregex(
    args.file, "M\s+(\S+)\s+(\S+)", [("x", np.float64), ("y", np.float64)]
)
fs = np.fromregex(
    args.file, "F\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", [("x0", np.float64), ("y0", np.float64), ("x1", np.float64), ("y1", np.float64)]
)
vs = np.fromregex(
    args.file, "V\s+(\S+)\s+(\S+)", [("x", np.float64), ("y", np.float64)]
)
cs = np.fromregex(
    args.file, "C\s+(\d+)\s+(\d+)", [("v0", np.int32), ("v1", np.int32)]
)

print(gs.shape, vs.shape)

for fi in range(len(fs)):
    f = fs[fi]
    x = [f["x0"], f["x1"]]
    y = [f["y0"], f["y1"]]
    pl.plot(x, y, "r-")
xlim = pl.gca().get_xlim()
ylim = pl.gca().get_ylim()

# periodic boundary conditions
for ix in range(-1, 2):
    for iy in range(-1, 2):
        if not ix == 0 or not iy == 0:
            for fi in range(len(fs)):
                f = fs[fi]
                x = np.array([f["x0"], f["x1"]]) + 2. * ix
                y = np.array([f["y0"], f["y1"]]) + 2. * iy
                pl.plot(x, y, "b-")

pl.plot(gs["x"], gs["y"], "g.")
# pl.plot(ms["x"], ms["y"], "y.")
# pl.plot(fs["x"], fs["y"], "r.")
# pl.plot(vs["x"], vs["y"], "r.")

pl.xlim(*xlim)
pl.ylim(*ylim)
if args.zoom:
    pl.xlim(0.0, 1.0e10)
    pl.ylim(0.0, 1.0e10)
pl.gca().set_aspect("equal")
# pl.show()
pl.savefig((Path(args.file)).parent / "vtest.png", dpi=300)
