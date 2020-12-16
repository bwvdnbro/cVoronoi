# Script that can be used to plot the Delaunay tessellation in test.txt.

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pl
import argparse

argparser = argparse.ArgumentParser()
argparser.add_argument("--circles", "-c", action="store_true", default=False)
argparser.add_argument("--zoom", "-z", action="store_true")
args = argparser.parse_args()

plotCircles = args.circles

vs = np.fromregex(
    "test.txt",
    "V\s+(\d+)\s+(\S+)\s+(\S+)",
    [("index", np.int32), ("x", np.float64), ("y", np.float64)],
)
ts = np.fromregex(
    "test.txt",
    "T\s+(\d+)\s+(\d+)\s+(\d+)",
    [("v0", np.int32), ("v1", np.int32), ("v2", np.int32)],
)

print(vs.shape, ts.shape)

pl.plot(vs["x"], vs["y"], "k.")
for ti in range(len(ts)):
    t = ts[ti]
    it = [t["v0"], t["v1"], t["v2"], t["v0"]]
    pl.plot(vs[it]["x"], vs[it]["y"], "k-")

    if plotCircles:
        p0 = vs[it[0]]
        p1 = vs[it[1]]
        p2 = vs[it[2]]
        a = [p1["x"] - p0["x"], p1["y"] - p0["y"]]
        b = [p2["x"] - p0["x"], p2["y"] - p0["y"]]
        D = 2.0 * (a[0] * b[1] - a[1] * b[0])
        a2 = a[0] * a[0] + a[1] * a[1]
        b2 = b[0] * b[0] + b[1] * b[1]
        R = [(b[1] * a2 - a[1] * b2) / D, (a[0] * b2 - b[0] * a2) / D]
        m = [p0["x"] + R[0], p0["y"] + R[1]]
        pl.plot(m[0], m[1], "g.")
        R = np.sqrt(R[0] * R[0] + R[1] * R[1])
        th = np.linspace(0.0, 2.0 * np.pi, 1000)
        pl.plot(m[0] + R * np.cos(th), m[1] + R * np.sin(th), "-")

if False:
    p0 = (0.250128, 0.206265)
    p1 = (0.863622, 0.636559)
    p2 = (0.0249239, 0.301656)
    a = [p1[0] - p0[0], p1[1] - p0[1]]
    b = [p2[0] - p0[0], p2[1] - p0[1]]
    D = 2.0 * (a[0] * b[1] - a[1] * b[0])
    a2 = a[0] * a[0] + a[1] * a[1]
    b2 = b[0] * b[0] + b[1] * b[1]
    R = [(b[1] * a2 - a[1] * b2) / D, (a[0] * b2 - b[0] * a2) / D]
    m = [p0[0] + R[0], p0[1] + R[1]]
    pl.plot(m[0], m[1], "g.")
    R = np.sqrt(R[0] * R[0] + R[1] * R[1])
    th = np.linspace(0.0, 2.0 * np.pi, 1000)
    pl.plot(m[0] + R * np.cos(th), m[1] + R * np.sin(th), "-")

if args.zoom:
    pl.xlim(0.0, 1.0e10)
    pl.ylim(0.0, 1.0e10)
pl.gca().set_aspect("equal")
# pl.show()
pl.savefig("test.png", dpi=300)
