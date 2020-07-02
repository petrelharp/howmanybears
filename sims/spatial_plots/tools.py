import scipy.sparse
import numpy as np

def relatedness_matrix(ts, left=0.0, right=None):
    """
    Constructs the sparse matrix whose [i,j]th entry gives the amount that node j
    inherited *directly* from node i, i.e., the sum of the length of all edges
    that have i as a parent and j as a child.
    """
    if right is None:
        right = ts.sequence_length
    edges = ts.tables.edges
    R = scipy.sparse.coo_matrix((np.fmin(right, edges.right) - np.fmax(left, edges.left), 
                           (edges.parent, edges.child)), 
                           shape = (ts.num_nodes, ts.num_nodes), dtype = 'float')
    return R.tocsc()


def relatedness(ts, focal_nodes, max_hops):
    """
    For each node, find the smallest number of genealogical hops to one of focal_nodes.
    """
    X = (ts.relatedness_matrix() > 0)
    Xt = X.transpose()
    out = np.repeat(np.inf, ts.num_nodes)
    out[focal_nodes] = 0
    x = np.repeat(0.0, ts.num_nodes)
    x[focal_nodes] = 1.0
    for n in range(1, max_hops + 1):
        # n is the number of up-hops
        x = X.dot(x)
        y = x.copy()
        out[y > 0] = np.fmin(out[y > 0], n)
        for k in range(1, max_hops + 1 - n):
            # k is the number of down-hops
            y = Xt.dot(y)
            # now y[j] is the number of paths of length n + k 
            #  that go from any focal node to j.
            out[y > 0] = np.fmin(out[y > 0], n + k)
    return out

