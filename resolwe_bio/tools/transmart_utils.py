import itertools
import numpy as np
from numpy.lib.recfunctions import merge_arrays

import utils


def format_annotations(annfile):
    attrs = annfile.next()[:-1].split('\t')
    annfile.next()
    annfile.next()
    annp = np.genfromtxt(annfile, dtype=None, delimiter='\t', names=['f{}'.format(i) for i in range(len(attrs))])

    def no_subpaths(path, paths):
        for _path in paths:
            if _path.startswith(path):
                return False
        return True

    # set column indexes
    attri = {a: i for i, a in enumerate(attrs)}

    # group leaves of the attribute tree
    attrg = {k: list(g) for k, g in itertools.groupby(sorted(attrs), lambda x: x[:x.rfind('_')])}
    attrg_keys = set(attrg.keys())

    # filter-out attributes: with a single leaf and with sub-attributes
    attrg = {k: g for k, g in attrg.iteritems() if len(g) > 1 and no_subpaths(k, attrg_keys.difference([k]))}

    # use numpy for smart indexing and column merging
    attrs_final = set(attrs[1:])
    attrs_merged = set()

    # merge groupped columns
    for k, g in attrg.iteritems():
        col_ndx = 'f{}'.format(attri[g[0]])
        col_merged = annp[col_ndx]

        add_func = None
        dtype_symbol = None

        if 'S' in str(annp.dtype[col_ndx]):
            add_func = np.core.defchararray.add
            dtype_symbol = 'S'
        elif 'int' in str(annp.dtype[col_ndx]):
            dtype_symbol = 'int'
            # TODO: Set add_func - no such example yet
        else:
            # TODO: No such example yet
            raise

        def stringify(col):
            if 'int' in str(col.dtype):
                return [str(x) if x > -1 else '' for x in col]
            return col

        # check if mixed types
        for col_name in g[1:]:
            if dtype_symbol not in str(annp.dtype['f{}'.format(attri[col_name])]):
                # convert mixed columns to string
                add_func = lambda x, y: np.core.defchararray.add(stringify(x), stringify(y))

        # merge columns
        for col_name in g[1:]:
            col_to_merge = annp['f{}'.format(attri[col_name])]
            col_merged = add_func(col_merged, col_to_merge)

        annp = merge_arrays((annp, col_merged), flatten=True)
        attri[k] = len(attri)
        attrs_final = attrs_final.difference(g)
        attrs_final.add(k)
        attrs_merged.add(k)

    # create var template
    var_template = []
    attrs_final = sorted(attrs_final)
    attrs_final_keys = map(utils.escape_mongokey, attrs_final)
    dtype = dict(annp.dtype.descr)
    dtype_final = []

    # format attribute type and var_template
    for attr_name in attrs_final:
        # format field label
        label = attr_name
        label_unsc_ind = attr_name.rfind('_')
        attr_key = 'f{}'.format(attri[attr_name])

        if label_unsc_ind >= 0:
            label = label[label_unsc_ind + 1:]

        label = label.replace('.', ' ')

        field, field_type = None, None
        field = {
            'name': utils.escape_mongokey(attr_name),
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

        if field is not None:
            var_template.append(field)
            dtype_final.append((attr_key, field_type or dtype[attr_key]))

    # transform original data to final attributes
    annp_final = annp.astype(dtype_final)

    # create final data set
    var_samples = {key: dict(zip(attrs_final_keys, map(np.asscalar, row))) for key, row in zip(annp['f0'], annp_final)}

    return var_samples, var_template
