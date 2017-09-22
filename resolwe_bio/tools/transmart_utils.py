# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""TranSMAT utility functions."""
from __future__ import absolute_import, division, print_function

import itertools
import numpy as np  # pylint: disable=import-error
from numpy.lib.recfunctions import merge_arrays  # pylint: disable=import-error

import utils


def format_annotations(annfile, treefile, ann_ids_file):
    """Format annotations."""
    # treefile = open('UBIOPRED_tree.txt', 'rb')
    # annfile = open('UBIOPRED_annotations_short.tab', 'rb')

    # R replaces some characters to . when it creates variables
    escape_chars = {
        '(': '\\uff08',
        ')': '\\uff09',
        '[': '\\uff3b',
        ']': '\\uff3d',
        ',': '\\uff0c',
        '-': '\\uff0d',
        '<': '\\uff1c',
        '>': '\\uff1e',
        '=': '\\uff1d',
        ':': '\\uff1a',
        '/': '\\uff0f',
        '^': '\\uff3e',
        '?': '\\uff1f',
        '\'': '\\uff07'
    }

    def esc_tree(l):
        return l.strip()[1:-1].replace(' ', '.').replace('\\', '_')

    def esc_dot(s):
        for k in escape_chars:
            s = s.replace(k, '.')
        return s

    tree = [esc_tree(l) for l in treefile if l.strip() != '']
    treefile.close()
    tree_original_names = {esc_dot(l): l for l in tree}

    ann_ids = [esc_dot(esc_tree(l)) for l in ann_ids_file.next()[:-1].split(';')]
    ann_ids_file.close()

    tree = sorted(tree_original_names.iterkeys())
    tree_attributes = set(tree)

    # some attributes are binarized in tranSMART and should be groupped together to a parent attribute
    # we detect the parent attributes automatically with 3 rules:
    #   1. the parent attribute is a node in the attribute tree with at least 2 leaves as childrens
    #   2. the parent attribute node's childrens should only be just leaves (not nodes)
    #   3. all children leaves are of string type
    # group leaves of the attribute tree
    group_candidates = {k: list(g) for k, g in itertools.groupby(tree, lambda x: x[:x.rfind('_')])}
    group_keys = set(group_candidates.keys())

    def no_subpaths(path, paths):
        for _path in paths:
            if _path.startswith(path):
                return False
        return True

    # filter-out attributes: with a single leaf and with sub-attributes
    group_candidates = {k: g for k, g in group_candidates.iteritems() if
                        len(g) > 1 and no_subpaths(k, group_keys.difference([k]))}

    attrs = annfile.next()[:-1].split('\t')

    def long_name(a):
        matching_long_names = [name for name in ann_ids if name.endswith('_' + a) or name == a]  # if short or long

        if len(matching_long_names) < 1:
            print('{{"proc.warn": "Attribute {} not found in the attribute tree and will be ignored"}}'.format(a))
            return ''
        if len(matching_long_names) > 1:
            print('{{"proc.error": "Attribute {} matched with multiple long names {}"}}'.format(
                a, matching_long_names))
        return matching_long_names[0]

    attrs = [long_name(a) for a in attrs[1:]]

    annfile.next()
    annfile.next()
    annp = np.genfromtxt(annfile, dtype=None, delimiter='\t', names=['f{}'.format(i) for i in range(len(attrs) + 1)])
    annfile.close()

    annp_ncols = len(annp[0])
    if annp_ncols != len(attrs) + 1:
        print('{{"proc.error": "Severe problem with column ID matching"}}')

    # set column indexes
    attri = {a: i + 1 for i, a in enumerate(attrs)}
    attri.pop('', None)

    # use numpy for smart indexing and column merging
    attrs = set(attrs)
    if '' in attrs:
        attrs.remove('')

    if len(attrs.difference(tree_attributes)) > 0:
        print('{{"proc.error": "All attributes should be in the attribute tree at this point"}}')

    attrs_merged = set()

    # find final collumns
    leaf_node_map = {}
    for k, v in group_candidates.items():
        for vv in v:
            leaf_node_map[vv] = k

    attr_group_candidates = set(leaf_node_map[a] for a in attrs if a in leaf_node_map)

    # merge groupped columns
    seqarrays = [annp]
    for k in attr_group_candidates:
        g = group_candidates[k]

        col_ndx = 'f{}'.format(attri[g[0]])
        col_merged = annp[col_ndx]

        # check rule 3 (all should be of string type)
        if not all('S' in str(annp.dtype['f{}'.format(attri[c])]) for c in g):
            print("Attribute with mixed types: {}".format(k))
            continue

        # merge columns
        for col_name in g[1:]:
            col_to_merge = annp['f{}'.format(attri[col_name])]
            col_merged = np.core.defchararray.add(col_merged, col_to_merge)

        seqarrays.append(col_merged)

        attri[k] = annp_ncols + len(seqarrays) - 2
        attrs = attrs.difference(g)
        attrs.add(k)
        attrs_merged.add(k)
        original_name = tree_original_names[g[0]]
        original_name_ndx = original_name.rfind('_')
        if original_name_ndx >= 0:
            original_name = original_name[:original_name_ndx]
        tree_original_names[k] = original_name

    annp = merge_arrays(tuple(seqarrays), flatten=True)

    # create var template
    var_template = []
    attrs = sorted(attrs)
    attrs_final_keys = []
    dtype_final = []
    dtype = dict(annp.dtype.descr)

    # format attribute type and var_template
    for attr_name in attrs:
        # format field label
        label = tree_original_names[attr_name]
        label_unsc_ind = label.rfind('_')
        attr_key = 'f{}'.format(attri[attr_name])

        if label_unsc_ind >= 0:
            label = label[label_unsc_ind + 1:]

        label = label.replace('.', ' ')

        field, field_type = None, None
        field = {
            'name': attr_name,
            'label': label,
        }

        if attr_name in attrs_merged:
            field['type'] = 'basic:string:'
            field['choices'] = [{'value': v, 'label': v} for v in np.unique(annp['f{}'.format(attri[attr_name])])]

            if 'S' not in dtype[attr_key]:
                field_type = '|S24'

        elif 'i' in dtype[attr_key] or 'u' in dtype[attr_key]:
            field['type'] = 'basic:integer:'
        elif 'f' in dtype[attr_key]:
            field['type'] = 'basic:decimal:'
        elif 'S' in dtype[attr_key]:
            field['type'] = 'basic:string:'

        if 'type' in field:
            var_template.append(field)
            attrs_final_keys.append(field['name'])
            dtype_final.append((attr_key, field_type or dtype[attr_key]))
        else:
            print('Unknown type {} of attribute {}'.format(dtype[attr_key], attr_name))

    # transform original data to final attributes
    annp_final = annp.astype(dtype_final)

    # create final data set
    def not_nan(val):
        return (isinstance(val, str) and val != '') or (isinstance(val, str) and not np.isnan(val))

    var_samples = {key: dict([(attr, val) for attr, val in zip(attrs_final_keys, map(np.asscalar, row)) if
                              not_nan(val)]) for key, row in zip(annp['f0'], annp_final)}

    return var_samples, var_template
