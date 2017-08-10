#!/usr/bin/env python3
# pylint: disable=missing-docstring,invalid-name,import-error
"""Hierarchical clustering computes clusters where leaves remain partially unordered.

After clustering leaf ordering can be applied within the constraints set by
hierarchical clustering. This file provides three implementations of leaf ordering:
o) knn: k-nearest neighbor with computational complexity O(n log n)
o) optimal: optimal leaf ordering with computational complexity O(n^4)
o) simulated_annealing: simulated annealing.

"""
import copy
import random
import numpy as np
import scipy

# k-nearest neighbor


def _swap_children(Z, i, j):
    """Swap left and right child nodes."""
    if int(Z[i][j]) >= len(Z) + 1:
        row = Z[int(Z[i][j]) - len(Z) - 1]
        row[0], row[1] = row[1], row[0]
        for k in range(2):
            _swap_children(Z, int(Z[i][j]) - len(Z) - 1, k)


def _get_condensed_from_square(i, j, n):
    """Get condensed (triangular) matrix index from square matrix index.

    Square matrix index is a pair (i, j). Matrix size is denoted by n.
    Condensed index is a single number.

    """
    assert i != j, "Error: no diagonal elements in condensed matrix."
    if i < j:
        i, j = j, i
    return int(n * j - j * (j + 1) / 2 + i - 1 - j)


def _get_child(Z, i, j, k):
    """Get child from node in position i in linkage.

    Parameters:
    j (0 or 1): left or right child.
    k (0 or 1): left of right child of (i, j).

    """
    if Z[i][j] < len(Z) + 1:
        return int(Z[i][j])
    return _get_child(Z, int(Z[i][j]) - len(Z) - 1, k, k)


def _get_children(Z, i):
    """Get genes (leaves!) from node in position i in linkage.

    Leaves are returned as a list in the following order:
    - leftmost leaf of the left child
    - rightmost leaf of the left child
    - leftmost leaf of the right child
    - rightmost leaf of the right child

    """
    return [[_get_child(Z, i, j, k) for k in range(2)] for j in range(2)]


def _is_ordered(Z, D, i):
    """Check if node in position i in linkage is correctly ordered."""
    children = _get_children(Z, i)
    d = [D[_get_condensed_from_square(children[0][j], children[1][k], len(Z) + 1)]
         for j in range(2) for k in range(2)]
    # ll_rl ll_rr lr_rl lr_rr
    for j in range(2):
        for k in range(2):
            if d[j + 2 * k] == min(d):
                return [k == 1, j == 0]


def knn(Z, D):
    """Order genes in linkage using k-nearest neighbor method.

    Computational complexity of this implementation is O(n log n).

    """
    for i in range(len(Z)):
        ordered = _is_ordered(Z, D, i)
        for k in range(2):
            if not ordered[k]:
                _swap_children(Z, i, k)
    return Z


# optimal leaf ordering

def _get_orderings(L, D, i, keep=None):
    """Compute optimal orderings from position i in linkage.

    Orderings are dicts with the following structure: keys of orderings are
    tuples (x, y) where x is the leftmost leaf in the subscluster and y the
    rightmost leaf. Values v of orderings are dicts where
    - v['dL'] is a list of linkage indices where nodes need to be swapped
    compared to the original linkage L to obtain the lowest cost
    permutation in which the leftmost leaf is x and the rightmost leaf is y
    - v['cost'] is the associated cost of this permutation.
    Parameter keep (int) denotes the number of the lowest cost permutations
    to keep in each round. Setting keep=2 is equivalent to knn method.

    """
    n_genes = len(L) + 1
    if L[i][0] < n_genes:
        left_orderings = {(L[i][0], L[i][0]): {'dL': [], 'cost': 0.0}}
    else:
        left_orderings = _get_orderings(L, D, int(L[i][0]) - len(L) - 1, keep)
    if L[i][1] < n_genes:
        right_orderings = {(L[i][1], L[i][1]): {'dL': [], 'cost': 0.0}}
    else:
        right_orderings = _get_orderings(L, D, int(L[i][1]) - n_genes, keep)
    orderings = {}
    for left_key, left_value in left_orderings.items():
        for right_key, right_value in right_orderings.items():
            key = (left_key[0], right_key[1])
            cost = left_value['cost'] + right_value['cost'] + \
                D[_get_condensed_from_square(left_key[1], right_key[0], n_genes)]
            if key not in orderings or cost < orderings[key]['cost']:
                orderings[key] = {
                    'dL': left_value['dL'] + right_value['dL'],
                    'cost': cost
                }
            key = (right_key[0], left_key[1])
            cost = left_value['cost'] + right_value['cost'] + \
                D[_get_condensed_from_square(right_key[1], left_key[0], n_genes)]
            if key not in orderings or cost < orderings[key]['cost']:
                orderings[key] = {
                    'dL': left_value['dL'] + right_value['dL'] + [i],
                    'cost': cost
                }
    if keep and len(orderings) > keep:
        return dict(sorted(orderings.items(), key=lambda x: x[1]['cost'])[:keep])
    return orderings


def optimal(L, D, keep=None):
    """Order genes in linkage using optimal leaf ordering.

    The method is described in Bar-Joseph et al., Fast optimal leaf ordering for
    hierarchical clustering, Bioinformatics (2001), 17: S22-S29.

    """
    orderings = _get_orderings(L, D, len(L) - 1, keep)
    best_ordering = sorted(orderings.values(), key=lambda item: item['cost'])[0]
    return np.array([L[i] if i not in best_ordering['dL'] else L[i][np.array([1, 0, 2, 3])]
                     for i in range(len(L))])


# simulated annealing

def _get_cost(Z, D):
    """Get cost of leaf ordering Z."""
    leaf_order = _get_leaf_order(Z)
    costs = [D[_get_condensed_from_square(leaf_order[i], leaf_order[i + 1], len(Z) + 1)]
             for i in range(len(Z))]
    return sum(costs)


def _get_random_permutation(Z):
    """Swap in-place a pair of children in linkage in a random node."""
    i = int(random.random() * len(Z))
    z = Z[i]
    z[0], z[1] = z[1], z[0]
    return Z


def _get_leaf_order(Z):
    """Get order of leaves in dendrogram."""
    return scipy.cluster.hierarchy.dendrogram(Z, no_plot=True)['leaves']


def simulated_annealing(Z, D, n_trials=1000):
    """Order genes in linkage using simulated annealing."""
    good_cost = _get_cost(Z, D)
    good_Z = copy.deepcopy(Z)
    for trial in range(n_trials):
        Z = _get_random_permutation(copy.deepcopy(good_Z))
        c = _get_cost(Z, D)
        if c / good_cost > 1.0 + 1.0 * (n_trials - trial) / n_trials:
            continue
        else:
            good_cost = c
            good_Z = copy.deepcopy(Z)
    return good_Z
