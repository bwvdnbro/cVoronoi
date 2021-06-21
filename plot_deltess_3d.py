import mpl_toolkits.mplot3d as a3
from matplotlib import pylab as pl
import numpy as np


def tetrahedra_to_triangles(ts):
    tris = np.zeros((4*len(ts), 3), int)
    for i, tetr in enumerate(ts):
        tris[4*i] = tetr[:3]
        tris[4*i + 1] = tetr[[0, 3, 1]]
        tris[4*i + 2] = tetr[[0, 3, 2]]
        tris[4*i + 3] = tetr[1:]
    return tris


def plot(vs, ts):
    axes = a3.Axes3D(pl.figure())
    vts = vs[tetrahedra_to_triangles(ts), :]
    tri = a3.art3d.Poly3DCollection(vts)
    tri.set_alpha(0.1)
    tri.set_color('grey')
    axes.add_collection3d(tri)
    axes.plot(vs[:, 0], vs[:, 1], vs[:, 2], 'ko')
    # axes.set_axis_off()
    axes.set_xlim([-.1, 1.1])
    axes.set_ylim([-.1, 1.1])
    axes.set_zlim([-.1, 1.1])
    pl.show()


def main(fname):
    vs = np.fromregex(
        fname,
        "V\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)",
        [("index", np.int32), ("x", np.float64), ("y", np.float64), ("z", np.float64)],
    )
    vs = np.array([[i for i in arr] for arr in vs])[:, 1:]
    ts = np.fromregex(
        fname,
        "T\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)",
        [("v0", np.int32), ("v1", np.int32), ("v2", np.int32), ("v3", np.int32)],
    )
    ts = np.array([[i for i in arr] for arr in ts]).astype(int)

    print(vs.shape, ts.shape)

    # filter dummy vertices
    ts = ts[np.all((ts < 1) | (ts > 4), axis=1), :]
    plot(vs, ts)


if __name__ == "__main__":
    main("test.txt")
