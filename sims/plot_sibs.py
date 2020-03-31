#!/usr/bin/env python3
import pyslim, msprime, tskit
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
outbase = ".".join(treefile.split(".")[:-1])

ts = pyslim.load(treefile)

def get_children(ts, targets=None, max_parent_time=ts.slim_generation):
    '''
    Returns a list of individual IDs, giving for each individual
    the IDs of their children.
    '''
    if targets is None:
        targets = np.arange(ts.num_individuals)
    # edges describing relationships between individuals
    edge_parent_indiv = ts.tables.nodes.individual[ts.tables.edges.parent]
    edge_child_indiv = ts.tables.nodes.individual[ts.tables.edges.child]
    indiv_edges = np.logical_and(edge_parent_indiv != tskit.NULL,
                                 edge_child_indiv != tskit.NULL)
    # edges where the parent was actually alive at the birth time of the child
    child_births = ts.individual_times[edge_child_indiv[indiv_edges]]
    parent_births = ts.individual_times[edge_parent_indiv[indiv_edges]]
    parent_deaths = parent_births - ts.individual_ages[edge_parent_indiv[indiv_edges]]
    alive_edges = indiv_edges.copy()
    alive_edges[indiv_edges] = (child_births >= parent_deaths)
    # individuals whose only ancestors are individuals
    has_all_parents = np.repeat(True, ts.num_individuals)
    has_all_parents[
            ts.tables.nodes.individual[
                ts.tables.edges.child[np.logical_not(indiv_edges)]]] = False
    # putting it all together
    good_edges = alive_edges.copy()
    good_edges[np.logical_not(has_all_parents)[edge_child_indiv]] = False
    # and hit only the targets
    is_target = np.repeat(False, ts.num_individuals)
    is_target[targets] = True
    good_edges[np.logical_not(is_target)[edge_child_indiv]] = False
    good_parents = ts.tables.nodes.individual[ts.tables.edges.parent[good_edges]]
    good_children = ts.tables.nodes.individual[ts.tables.edges.child[good_edges]]
    out = []
    for k in range(ts.num_individuals):
        out.append(list(set(good_children[good_parents == k])))
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

    fig.savefig(f"{outbase}.sibs.pdf")

xy = np.column_stack([locs[parent_ids, :], locs[child_ids, :]])
np.savetxt(f"{outbase}.sib_locs.tsv", xy)
np.savetxt(f"{outbase}.sample_locs.tsv", locs[sampled_bears,:])
