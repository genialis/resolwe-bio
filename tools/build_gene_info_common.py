#!/usr/bin/env python2
import re

gene_id_re = re.compile('gene_id "([\w\-\.]*)"')
gene_name_re = re.compile('gene_name "([\w\-\.]*)"')
transcript_id_re = re.compile('transcript_id "([\w\-\.]*)"')
transcript_name_re = re.compile('transcript_name "([\w\-\.]*)"')
ensembl_id_re = re.compile('Ensembl:([\w]*)')
mgi_id_re = re.compile('MGI:([\w\:]*)')
omim_id_re = re.compile('MIM:([\w]*)')


def _search(regex, string):
    match = regex.search(string)
    return match.group(1) if match else 'N/A'


def get_gene_id(ids):
    return _search(gene_id_re, ids)


def get_gene_name(ids):
    return _search(gene_name_re, ids)


def get_transcript_name(ids):
    return _search(transcript_name_re, ids)


def get_transcript_id(ids):
    return _search(transcript_id_re, ids)


def get_ensembl_id(ids):
    return _search(ensembl_id_re, ids)


def get_omim_id(ids):
    return _search(omim_id_re, ids)


def get_mgi_id(ids):
    return _search(mgi_id_re, ids)
