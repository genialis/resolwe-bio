""".. Ignore pydocstyle D400.

==============================
Utility Functions for Commands
==============================

"""
from __future__ import absolute_import, division, print_function, unicode_literals

import gzip
import io
import os
import zipfile


def decompress(file_name):
    """Compression-agnostic iterator.

    Iterate over files on the archive and return a tuple of
    file name, line count and file descriptor.

    Supported file formats are .tab, .gz and .zip.

    """
    if not os.path.isfile(file_name):
        raise ValueError("Can not find file '{}'".format(file_name))

    line_count = -1
    _, ext = os.path.splitext(file_name)

    _open = None
    if ext == '.tab':
        _open = open
    elif ext == '.gz':
        _open = gzip.open

    if _open:
        with _open(file_name, 'rt') as tsv_file:
            line_count = sum(1 for row in tsv_file)

        with _open(file_name, 'rt') as tsv_file:
            yield (os.path.basename(file_name), line_count, tsv_file)
    elif ext == '.zip':
        with zipfile.ZipFile(file_name) as archive:
            for entry in archive.infolist():
                if not entry.filename.endswith('.tab'):
                    continue

                if entry.filename.startswith('__MACOSX'):
                    continue

                with archive.open(entry) as tsv_file:
                    line_count = sum(1 for row in tsv_file)

                with archive.open(entry) as tsv_file:
                    yield (entry.filename, line_count, io.TextIOWrapper(tsv_file))
    else:
        raise ValueError("Unsupported file format")
