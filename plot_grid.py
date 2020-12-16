# Script that can be used to plot the Voronoi grid contained in vtest.txt

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as pl
import argparse

argparser = argparse.ArgumentParser()
argparser.add_argument("--zoom", "-z", action="store_true")
args = argparser.parse_args()

gs = np.fromregex("vtest.txt", "G\s+(\S+)\s+(\S+)",
                 [("x", np.float64), ("y", np.float64)])
ms = np.fromregex("vtest.txt", "M\s+(\S+)\s+(\S+)",
                 [("x", np.float64), ("y", np.float64)])
fs = np.fromregex("vtest.txt", "F\s+(\S+)\s+(\S+)",
                 [("x", np.float64), ("y", np.float64)])
vs = np.fromregex("vtest.txt", "V\s+(\S+)\s+(\S+)",
                 [("x", np.float64), ("y", np.float64)])
cs = np.fromregex("vtest.txt", "C\s+(\d+)\s+(\d+)",
                 [("v0", np.int32), ("v1", np.int32)])

print(gs.shape, vs.shape)

for ci in range(len(cs)):
  c = cs[ci]
  it = [c["v0"], c["v1"]]
  pl.plot(vs[it]["x"], vs[it]["y"], "r-")
xlim = pl.gca().get_xlim()
ylim = pl.gca().get_ylim()

for ix in range(-1, 2):
  for iy in range(-1, 2):
    if not ix == 0 or not iy == 0:
      for ci in range(len(cs)):
        c = cs[ci]
        it = [c["v0"], c["v1"]]
        pl.plot(vs[it]["x"]+ix*1.e10, vs[it]["y"]+iy*1.e10, "b-")

pl.plot(gs["x"], gs["y"], "g.")
#pl.plot(ms["x"], ms["y"], "y.")
#pl.plot(fs["x"], fs["y"], "r.")
#pl.plot(vs["x"], vs["y"], "r.")

pl.xlim(*xlim)
pl.ylim(*ylim)
if args.zoom:
  pl.xlim(0., 1.e10)
  pl.ylim(0., 1.e10)
pl.gca().set_aspect("equal")
#pl.show()
pl.savefig("vtest.png", dpi=300)
