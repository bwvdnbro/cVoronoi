# Script that can be used to plot the Voronoi grid contained in vtest.txt

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pl
import argparse


def check_vertices(v0, v1, v2, gs):
    return (
        v0 >= 0
        and v0 < len(gs)
        and v1 >= 0
        and v1 < len(gs)
        and v2 >= 0
        and v2 < len(gs)
    )


argparser = argparse.ArgumentParser()
argparser.add_argument("--zoom", "-z", action="store_true")
args = argparser.parse_args()

do_delaunay = False

paths = [np.zeros((100, 4)), np.zeros((100, 4))]
for i in range(100):
    print(i)
    vname = "vtest{0:03d}.txt".format(i)
    gs = np.fromregex(
        vname, "G\s+(\S+)\s+(\S+)", [("x", np.float64), ("y", np.float64)]
    )
    vs = np.fromregex(
        vname, "V\s+(\S+)\s+(\S+)", [("x", np.float64), ("y", np.float64)]
    )
    cs = np.fromregex(
        vname, "C\s+(\d+)\s+(\d+)", [("v0", np.int32), ("v1", np.int32)]
    )

    pl.plot(
        [0.0, 0.0, 1.0e10, 1.0e10, 0.0], [0.0, 1.0e10, 1.0e10, 0.0, 0.0], "k--"
    )

    if do_delaunay:
        dname = "test{0:03d}.txt".format(i)
        ts = np.fromregex(
            dname,
            "T\s+(\d+)\s+(\d+)\s+(\d+)",
            [("v0", np.int32), ("v1", np.int32), ("v2", np.int32)],
        )
        ts["v0"] -= 3
        ts["v1"] -= 3
        ts["v2"] -= 3

        for ti in range(len(ts)):
            t = ts[ti]
            if check_vertices(t["v0"], t["v1"], t["v2"], gs):
                it = [t["v0"], t["v1"], t["v2"], t["v0"]]
                pl.plot(gs[it]["x"], gs[it]["y"], "k-")

    for ci in range(len(cs)):
        c = cs[ci]
        it = [c["v0"], c["v1"]]
        pl.plot(vs[it]["x"], vs[it]["y"], "r-")

    pl.plot(gs["x"], gs["y"], "g.")

    #  r = np.sqrt((gs["x"]-5.e9)**2 + (gs["y"]-5.e9)**2)
    #  isel = r < 1.e9
    #  print(isel.nonzero())
    #  pl.plot(gs[isel]["x"], gs[isel]["y"], "y.")

    isel = [12, 41, 67, 73]
    pl.plot(gs[isel]["x"], gs[isel]["y"], "b.")

    paths[0][i, :] = gs[isel]["x"]
    paths[1][i, :] = gs[isel]["y"]

    pl.plot(
        paths[0][: (i + 1), :],
        paths[1][: (i + 1), :],
        "b:",
        zorder=-99,
        alpha=0.4,
    )

    pl.xlim(-1.0e9, 1.1e10)
    pl.ylim(-1.0e9, 1.1e10)

    pl.gca().set_aspect("equal")
    pl.gca().axis("off")
    pl.tight_layout()
    pl.savefig("vtest{0:03d}.png".format(i), dpi=300)
    pl.close()
