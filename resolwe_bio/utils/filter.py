"""Test filter functions."""


def filter_vcf_variable(line):
    """Filter variable lines from the VCF file."""
    if line.startswith(b"##samtoolsVersion"):
        return True
    elif line.startswith(b"##reference"):
        return True
    elif line.startswith(b"##fileDate"):
        return True
    elif b"/data_all/" in line:
        return True
