import mpl_toolkits.mplot3d as a3
from matplotlib import pylab as pl
import numpy as np


def plot_voronoi(generators, vertices):
    fig = pl.figure()
    axes = fig.add_subplot(111, projection='3d')
    # axes = a3.Axes3D(pl.figure())
    poly3dcollection = a3.art3d.Poly3DCollection(vertices, facecolors="g", linewidth=1, alpha=0.3)
    poly3dcollection.set_edgecolor("k")
    # poly3dcollection.set_alpha(1)
    # poly3dcollection.set_color('grey')
    axes.add_collection3d(poly3dcollection)
    axes.plot(generators[:, 0], generators[:, 1], generators[:, 2], 'ko')
    # axes.set_axis_off()
    axes.set_xlim([-.1, 1.1])
    axes.set_ylim([-.1, 1.1])
    axes.set_zlim([-.1, 1.1])
    pl.show()


def main(fname):
    with open(fname, "r") as file:
        lines = file.readlines()
    lines = [line[:-1].split("\t") for line in lines]

    generators = np.stack([np.array(line[1:]) for line in lines if line[0] == "G"]).astype(float)

    centroids = [line for line in lines if line[0] == "C"]
    volumes = np.array([np.array(line[-2]) for line in centroids]).astype(float)
    n_neighbours = np.array([np.array(line[-1]) for line in centroids]).astype(int)
    centroids = np.stack([np.array(line[1:-2]) for line in centroids]).astype(float)

    faces = [line[1:] for line in lines if line[0] == "F"]
    sid = np.array([np.array(line[0]) for line in faces]).astype(int)
    areas = np.array([np.array(line[1]) for line in faces]).astype(float)
    midpoints = np.stack([np.array(line[2:5]) for line in faces]).astype(float)
    vertices = [np.stack([np.array(c[1:-1].split(", ")) for c in line[5:]]).astype(float) for line in faces]

    plot_voronoi(generators, vertices)


if __name__ == "__main__":
    main("vtest.txt")
