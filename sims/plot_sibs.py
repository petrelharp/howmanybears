import pyslim, msprime
import numpy as np
import spatial_plots as sps
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.collections as cs

import scipy.sparse as sparse

np.random.seed(32)

usage = """
Usage:
    {} (treefile name)
""".format(sys.argv[0])

if len(sys.argv) != 2:
    raise ValueError(usage)

treefile = sys.argv[1]

ts = pyslim.load(treefile)

def get_children(ts, targets=None, max_parent_time=ts.slim_generation):
    '''
    Returns a list of individual IDs, giving for each individual
    the IDs of their children.
    '''
    if targets is None:
        targets = np.arange(ts.num_individuals)
    not_targets = set(range(-1, ts.num_individuals)) - set(targets)
    out = []
    nodes = ts.tables.nodes
    edges = ts.tables.edges
    for ind in ts.individuals():
        children = []
        if ind.time < max_parent_time:
            for n in ind.nodes:
                children.extend(nodes.individual[edges.child[edges.parent == n]])
            children = set(children) - not_targets
        out.append(list(children))
    return out


sampled_bears = np.random.choice(ts.individuals_alive_at(0), 1000)
children = get_children(ts, targets=sampled_bears, max_parent_time=100)
num_children = np.array([len(x) for x in children])

parent_ids = []
child_ids = []
for par, ch in enumerate(children):
    if len(ch) > 1:
        for c in ch:
            parent_ids.append(par)
            child_ids.append(c)

parent_ids = np.array(parent_ids)
child_ids = np.array(child_ids)

locs = ts.individual_locations[:,:2]
xmax = max(locs[:, 0])
ymax = max(locs[:, 1])

if False:
    fig, ax = plt.subplots(figsize=(6, 6 * ymax / xmax))
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)

    sps.plot_density(locs[ts.individuals_alive_at(0), :], xmax, ymax, ax, scatter=False)
    ax.scatter(locs[sampled_bears, 0], locs[sampled_bears, 1], s=7, edgecolors='none', facecolors='black', alpha=0.2)

    ax.quiver(locs[parent_ids, 0],  # X
              locs[parent_ids, 1],  # Y
              locs[child_ids, 0] - locs[parent_ids, 0],  # dX
              locs[child_ids, 1] - locs[parent_ids, 1],  # dY
              color = 'black',
              alpha = 0.6,
              width = 0.75,
              units='xy', scale=1)

    fig.savefig("sibs.pdf")

xy = np.column_stack([locs[parent_ids, :], locs[child_ids, :]])
np.savetxt("sib_locs.tsv", xy)
np.savetxt("sample_locs.tsv", locs[sampled_bears,:])
