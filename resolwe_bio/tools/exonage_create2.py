"""
    arg 1: samtools path
    arg 2: gene name
    arg 3: GTF file of exons
    arg 4: BAM file path
    arg 5: output write file
    arg 6: warning messages
"""
import sys
from subprocess import PIPE, Popen


class gene:
    """
    For storing exon information of GTF FILE!
    """
    def __init__(self, ID):
        self.id = ID
        self.chr = None
        self.exons = {}

    def add_exon(self, chr, transcript_id, start, end):
        self.exons[(start, end)] = self.exons.get((start, end), []) + [transcript_id]
        self.chr = chr

    def __str__(self):
        return self.id


def cmdline(command):
    process = Popen(
        args=command,
        stdout=PIPE,
        shell=True
    )
    return process.communicate()[0]


def parse_bam_header(mapping):
    seqs = {}
    handle = mapping.split('\n')
    for x in handle:
        if not x:
            continue
        values = x.split()
        if values[0] == '@SQ':
            seq = values[1].split(':')[1]
            seqs[seq] = 1
    return seqs


def parse_and_build(gtf_file, gene_name, contigs):
    """
    :param gtf_file: gtf file location
    :return: list of genes (class gene) with exon locations
    """
    gene_construct = gene(gene_name)
    genes_pool = []
    handler = open(gtf_file, "r")
    for line in handler.readlines():
        line = line.split("\t")
        d = dict(map(str.strip, v.strip().split(" ")) for v in line[-1].split(";") if " " in v)
        gene_id = d["gene_name"][1:-1] if "gene_name" in d else d["gene_id"][1:-1]
        if gene_id != gene_name:
            continue
        if not line[0] in contigs:
            global error_message
            error_message = "Contig {} not found in BAM file (for gene {})".format(line[0], gene_id)
            continue
        transcript_id = d["transcript_id"][1:-1]
        gene_construct.add_exon(line[0], transcript_id, int(line[3]), int(line[4]))
    return gene_construct

error_message = ""
samtools_path = sys.argv[1]
bam_file = sys.argv[4]
bam_header = cmdline("{} view -H {}".format(samtools_path, bam_file))
contigs = parse_bam_header(bam_header)

gene_object = parse_and_build(sys.argv[3], sys.argv[2], contigs)
print 'analysis for gene {}'.format(gene_object.id)

counter = 0
# determinating start and stop of gene
if not gene_object.exons:
    exit()

start = min([(x[0]) for x in gene_object.exons])
end = max([x[1] for x in gene_object.exons])

# calling samtools for output in gene range
command = "{} depth -r {}:{}-{} {}".format(samtools_path, gene_object.chr, start, end, bam_file)

my_file = file('temporary_gene_coverage.txt', 'w')
my_file.write(cmdline(command))
my_file.close()

output_file = file(sys.argv[5], 'w')

# loop for each exon
for crange, transcripts in sorted(gene_object.exons.items()):
    # getting output for exon
    printer = cmdline("awk ' $2>={} && $2<={}' temporary_gene_coverage.txt".format(crange[0], crange[1]))
    # print for each (exon containg) transcript
    for transript in transcripts:
        counter += 1
        if printer.strip():
            output_file.write(printer)
        output_file.write("EXON\t{}\t{}\t{}\t{}-{}\t{}\n".format(
            gene_object.id, transript, counter, crange[0], crange[1], crange[1] - crange[0]))

output_file.close()
my_file = file(sys.argv[6], 'w')
my_file.write(error_message)
my_file.close()
