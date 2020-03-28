import scipy.stats
import numpy as np

from .tools import *

def plot_density(locs, xmax, ymax, ax, scatter=True, alpha=0.8):
    """
    Plot a 2D kernel density estimate of the population density.
    """
    tlocs = locs.T
    kde = scipy.stats.gaussian_kde(tlocs)
    X, Y = np.meshgrid(
            np.linspace(0.0, xmax, 51),
            np.linspace(0.0, ymax, 51))
    Z = kde([X.flatten(), Y.flatten()])
    Z.shape = X.shape
    if scatter:
        ax.scatter(locs[:, 0], locs[:, 1],
                   s=10,
                   alpha=0.5,
                   c='black',
                   marker="o",
                   edgecolors='none')
    ax.contour(X, Y, Z,
               colors='c',
               alpha=alpha,
               zorder=-1)


