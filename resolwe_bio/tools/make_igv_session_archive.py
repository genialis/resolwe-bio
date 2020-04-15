#!/usr/bin/env python3
"""Create IGV session for archive."""

import argparse
import os

from lxml import etree

parser = argparse.ArgumentParser(description="Create igv session.")

parser.add_argument(
    "-f", "--input_file", required=True, help="File with paths to files for IGV."
)

args = parser.parse_args()


def get_build_info(input_file):
    """Get build information in UCSC notation, if possible."""
    build_dict = {
        "hg38": ["GRCh38"],
        "hg19": ["GRCh37", "b37"],
        "mm10": ["GRCm38"],
        "mm9": ["MGSCv37"],
        "rn6": ["Rnor_6.0"],
    }

    with open(input_file, "r") as tfile:
        build = tfile.readline().strip()

    new_build = [
        k for k, v in build_dict.items() if k == build or build.startswith(tuple(v))
    ]
    if new_build:
        return new_build[0]
    return ""


def make_xml_tree(input_file):
    """Make xml tree for IGV session."""
    global_ = etree.Element("Global", genome=get_build_info(input_file), version="3",)

    resources = etree.SubElement(global_, "Resources")

    with open(input_file, "r") as tfile:
        next(tfile)  # Skip the line with the build.
        for filename in tfile:
            filename = filename.rstrip()
            #  replace None (dir folder if species and build are not defined) with other_data
            if os.path.dirname(os.path.dirname(filename)) == "None":
                filename = filename.replace("None", "other_data")
            etree.SubElement(
                resources,
                "Resource",
                name=os.path.basename(filename),
                path=os.path.join("..", filename),
            )

    doc = etree.tostring(
        global_, pretty_print=True, xml_declaration=True, encoding="UTF-8",
    )
    return doc


def write_xml_file(input_file, doc):
    """Compose a name and write output file."""
    os.makedirs("IGV", exist_ok=True)
    output_name = os.path.join("IGV", input_file.replace("temp_igv.txt", "igv.xml"))
    with open(output_name, "wb") as f:
        f.write(doc)


def main():
    """Invoke when run directly as a program."""
    doc = make_xml_tree(args.input_file)
    write_xml_file(args.input_file, doc)


if __name__ == "__main__":
    main()
