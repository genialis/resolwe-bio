#!/usr/bin/env python2
# pylint: disable=missing-docstring,invalid-name
# XXX: Refactor to a comand line tool and remove pylint disable
"""Build gene info helper module."""
import re

gene_id_re = re.compile(r'gene_id "([\w\-\.]*)"')
gene_name_re = re.compile(r'gene_name "([\w\-\.]*)"')
transcript_id_re = re.compile(r'transcript_id "([\w\-\.]*)"')
transcript_name_re = re.compile(r'transcript_name "([\w\-\.]*)"')
ensembl_id_re = re.compile(r'Ensembl:([\w]*)')
mgi_id_re = re.compile(r'MGI:([\w\:]*)')
omim_id_re = re.compile(r'MIM:([\w]*)')


def _search(regex, string):
    match = regex.search(string)
    return match.group(1) if match else 'N/A'


def get_gene_id(ids):
    """Get gene id."""
    return _search(gene_id_re, ids)


def get_gene_name(ids):
    """Get gene name."""
    return _search(gene_name_re, ids)


def get_transcript_name(ids):
    """Get transcript name."""
    return _search(transcript_name_re, ids)


def get_transcript_id(ids):
    """Get transcript id."""
    return _search(transcript_id_re, ids)


def get_ensembl_id(ids):
    """Get ensembl id."""
    return _search(ensembl_id_re, ids)


def get_omim_id(ids):
    """Get omim id."""
    return _search(omim_id_re, ids)


def get_mgi_id(ids):
    """Get mgi id."""
    return _search(mgi_id_re, ids)
