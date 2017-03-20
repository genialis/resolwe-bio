#!/usr/bin/env python3
"""Run iCount demultiplexing and annotate samples."""
import argparse
import gzip
import resdk
import re

from time import sleep


GENOME_MAP = {
    'Hs': 'star-index-ens-hs',
    'Mm': 'star-index-ens-mm'
}

SEGMENTATION_MAP = {
    'Hs': 'icount-segmentation-hs',
    'Mm': 'icount-segmentation-mm'
}

SPECIES_MAP = {
    'Hs': 'Homo sapiens',
    'Mm': 'Mus musculus'
}


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Run iCount demultiplexing.")
    parser.add_argument('annotation_file', help="iCount sample annotation file.")
    parser.add_argument('annotation_id', type=int, help="Annotation file object ID")
    parser.add_argument('reads_id', type=int, help="Reads file object ID")
    parser.add_argument('data_id', type=int, help="Current data ID")
    return parser.parse_args()


def _tsv_to_dict(table_file):
    """Parse *.tsv file into dict."""
    data = {}
    with gzip.open(table_file, 'rt') as tfile:
        header = next(tfile).strip().split('\t')
        column_picker = {label: index for index, label in enumerate(header)}
        for line in tfile:
            line = line.split('\t')
            brc = re.sub(r"[^A-Za-z]+", '', line[column_picker["5' barcode"]])
            data[brc] = {key: line[column_picker[key]] for key in column_picker}
    return data


def _get_barcode(name):
    """Get barcode name."""
    return name.split('.')[0].split('_')[1]


def create_sample_descriptor(annotation_data, bc):
    """Create sample descriptor."""
    sample_descriptor = {
        'sample': {
            'annotator': annotation_data[bc]['Scientist'],
            'pi': annotation_data[bc]['PI'],
            'source': annotation_data[bc]['cells/tissue'],
            'organism': SPECIES_MAP[annotation_data[bc]['mapto']]
        }
    }
    return sample_descriptor


def create_reads_descriptor(annotation_data, bc):
    """Create reads descriptor."""
    exp_name = annotation_data[bc]['experiment name'] if annotation_data[bc]['experiment name'] != 'missing' else ''
    reads_descriptor = {
        'protocols': {
            'method': annotation_data[bc]['Method'],
            'protein': annotation_data[bc]['Protein'],
            'mapto': annotation_data[bc]['mapto']
        },
        'reads_info': {
            '5_barcode': bc,
            '3_adapter': ''.join([adpr.split('_')[0] for adpr in annotation_data[bc]["3' adapter"].split(',')]),
            'sequencer_type': annotation_data[bc]['sequencer type'],
            'sequencing_length': annotation_data[bc]['sequencing length (optional)'],
            'linker': annotation_data[bc]['linker'],
            'rt_primer': annotation_data[bc]['RT primer'],
            'antibody': annotation_data[bc]['antibody used']
        },
        'other': {
            'resequencing_id': annotation_data[bc]['resequencing ID (optional)'],
            'experiment_name': exp_name,
            'condition': annotation_data[bc]['condition (optional)'],
            'replicate': annotation_data[bc]['replicate (optional)'],
            'pmid': annotation_data[bc]['PMID of published data (optional)'],
            'mw_band_cut': annotation_data[bc]['MW band cut (kDa) (optional)'],
            'gel_image_date': annotation_data[bc]['date of gel images in lab notebook (optional)'],
            'notebook_page': annotation_data[bc]['page of gel images in lab notebook (optional)'],
            'comments': annotation_data[bc]['comments (optional)']
        }
    }
    return reads_descriptor


def main():
    """Invoked when run directly as a program."""
    args = parse_arguments()

    res = resdk.Resolwe()

    parent = res.data.get(args.data_id)
    parent_user = parent.contributor['username']

    permissions = {
        'users': {
            'add': {
                parent_user: ['download', 'edit', 'owner', 'share', 'view']
            }
        }
    }

    demux_data = _tsv_to_dict(args.annotation_file)
    barcodes = demux_data.keys()

    adapters_temp = list(set([demux_data[barcode]["3' adapter"] for barcode in barcodes]))
    if len(adapters_temp) != 1:
        print('{"proc.error":"Only one adapter type is allowed per demultiplexing job."}')
        exit(1)

    adapter = ''.join([adpr.split('_')[0] for adpr in adapters_temp[0].split(',')])

    demux_inputs = {
        'reads': res.data.get(args.reads_id),
        'annotation': res.data.get(args.annotation_id)
    }

    demux = res.run('icount-demultiplex-samples', input=demux_inputs, collections=parent.collections)

    # set demux permissions
    demux.api(demux.id).permissions.post(permissions)

    while demux.status not in ['ER', 'OK']:
        sleep(10)
        demux.update()

    if demux.status == 'ER':
        print('{"proc.error":"Only one adapter type is allowed per demultiplexing job."}')
        exit(1)

    demux_reads = res.data.filter(parents=demux.id)

    for read in demux_reads:
        read.api(read.id).permissions.post(permissions)
        if read.sample:
            barcode = _get_barcode(read.name)
            if barcode in demux_data.keys():
                print("ANNOTATING SAMPLE:", read.sample.name)
                read.descriptor_schema = 'reads-iclip'
                read.descriptor = create_reads_descriptor(demux_data, barcode)
                read.save()
                read.sample.descriptor_schema = 'sample-iclip'
                read.sample.descriptor = create_sample_descriptor(demux_data, barcode)
                read.sample.name = demux_data[barcode]['Renamed fastq on iCount']
                read.sample.slug = demux_data[barcode]['Renamed fastq on iCount']
                read.sample.save()
                read.sample.confirm_is_annotated()

                for pc in parent.collections:
                    pc.add_samples(read.sample)

                print("RUNNING PRIMARY ANALYSIS FOR SAMPLE:", read.sample.name)

                genome = res.data.filter(slug=GENOME_MAP[demux_data[barcode]['mapto']])
                segmentation = res.data.filter(slug=SEGMENTATION_MAP[demux_data[barcode]['mapto']])

                if len(genome) == 0:
                    print('{"proc.warning":"STAR index not found in the database."}')

                if len(segmentation) == 0:
                    print('{"proc.warning":"Segmentation map not found in the database."}')

                if len(genome) == 1 and len(segmentation) == 1:
                    icount_inputs = {
                        'reads': read.id,
                        'index': res.data.get(slug=GENOME_MAP[demux_data[barcode]['mapto']]),
                        'segmentation': res.data.get(SEGMENTATION_MAP[demux_data[barcode]['mapto']])
                    }

                    icount_analysis = res.run('workflow-icount', input=icount_inputs, collections=parent.collections)
                    icount_analysis.api(icount_analysis.id).permissions.post(permissions)

                    for child in res.data.filter(parents=icount_analysis.id):
                        child.api(child.id).permissions.post(permissions)

                else:
                    print('{"proc.warning":"Primary analysis not triggered."}')


if __name__ == "__main__":
    main()
