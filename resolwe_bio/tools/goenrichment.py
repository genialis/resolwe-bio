import argparse
import cPickle as pickle
import json
import os
import utils

from collections import defaultdict

from orangecontrib.bio.obiGO import Annotations, Ontology


parser = argparse.ArgumentParser(description="Gene Ontology Enrichment Analysis")
parser.add_argument('ontology', help="gene ontology file")
parser.add_argument('annotation', help="GO annotation file (GAF)")
parser.add_argument('orthologues', help="File defining gene orthlogues relations")
parser.add_argument('genes', nargs='+', help="list of gene ids")
parser.add_argument('--min_genes', type=int, default=2, help="minimum number of genes in GO term")
parser.add_argument('--p_threshold', type=float, default=0.2, help="maximum p-value")
args = parser.parse_args()

if not os.path.isfile(args.ontology):
    raise ValueError("File {} does not exist".format(args.ontology))

if not os.path.isfile(args.annotation):
    raise ValueError("File {} does not exist".format(args.annotation))


def get_enriched_parents(term, terms, ontology):
    """Return parents of a term that are included in provided term subset"""
    return [id for (_, id) in ontology.terms[term].related if id in terms]


def get_siblings_from_parents(parents):
    """Return a dictionary of siblings from parent dictionary."""
    siblings = defaultdict(set)
    for term, term_parents in parents.items():
        for term_parent in term_parents:
            siblings[term_parent].add(term)
    return siblings


def get_clean_parents(terms, ontology):
    """
    For a group of terms returns a dictionary of clean parents (parents
    that are not ancestors of term's parents). Also returns
    a dictionary of parents and siblings.
    """
    parents = {term: set(get_enriched_parents(term, terms, ontology))
               for term in terms}

    clean_parents = defaultdict(list)
    for term in terms:
        clean_parents[term] = parents[term] - \
            set().union(*[parents[p] for p in parents[term]])
    return clean_parents


def tree_walk(t):
    """Return a GO tree structure."""
    children = []
    if t in siblings and siblings[t]:
        children = [tree_walk(s) for s in siblings[t]]

    return {
        "pval": terms[t][1],
        "matched": len(terms[t][0]),
        "total": terms[t][2],
        "score": (len(terms[t][0]) / cluster_size) / (terms[t][2] / ref_size),
        "term_name": annotations.ontology.terms[t].name,
        "term_id": t,
        "gene_ids": [translator.get(name, name) for name in terms[t][0]],
        "children": children
    }


def abbreviate_namespace(namespace):
    return (namespace[0] + namespace[namespace.index('_') + 1]).upper()


##########################################################################
# Main
split = os.path.split
ontology_id = split(split(args.ontology)[0])[1]
annotation_id = split(split(args.annotation)[0])[1]
annotation_cache = os.path.join('/tmp', 'GO_' + ontology_id + '_Annotation_' + annotation_id)

if os.path.isfile(annotation_cache):
    try:
        with open(annotation_cache, 'rb') as fd:
            annotations = pickle.load(fd)
    except:
        os.remove(annotation_cache)
        raise

else:
    with utils.gzopen(args.ontology) as fd:
        ontology = Ontology(fd)

    annotations = Annotations(file=args.annotation, ontology=ontology)

    with open(annotation_cache, 'wb') as fd:
        pickle.dump(annotations, fd, -1)

translator = {a.DB_Object_Symbol: a.DB_Object_ID for a in annotations}

orth = {}
genes = set()

if args.orthologues:
    orth = dict(l.strip().split("\t") for l in utils.gzopen(args.orthologues))

for g in args.genes:
    if g in annotations.aliasMapper:
        genes.add(g)
    else:
        if g in orth and orth[g] in annotations.aliasMapper:
            genes.add(orth[g])

terms = annotations.GetEnrichedTerms(genes)

# p-value filter
terms = {k: term for (k, term) in terms.items()
         if term[1] <= args.p_threshold and len(term[0]) >= args.min_genes}

print "Cluster group: %d genes" % len(args.genes)
print "Terms with genes from cluster group: %d" % len(terms)
print "Filtered terms (p<=%.3f, n>=%d): %d" % (args.p_threshold, args.min_genes, len(terms))

parents = get_clean_parents(terms, annotations.ontology)
siblings = get_siblings_from_parents(parents)  # siblings from clean parents
roots = [term for term in parents if not parents[term]]

ref_size = float(len(annotations.geneNames))
cluster_size = float(len(args.genes))

tree = {'BP': [], 'CC': [], 'MF': []}
for root in roots:
    tree[abbreviate_namespace(annotations.ontology.terms[root].namespace)].append(tree_walk(root))

print json.dumps({"terms": tree}, separators=(',', ':'))
