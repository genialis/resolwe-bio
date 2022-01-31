#!/usr/bin/env python3

"""A script that checks if files in a directory exceed a size limit."""

import argparse
import os


def find_large_files(directory, limit):
    """Find files in directory that exceed the size limit.

    :param str directory: path of directory where to find large files

    :param float limit: size limit (in MBs) that determines whether a
        file is a large file

    :return: list of large files
    :rtype: list

    """
    large_files = []
    for f in sorted(os.listdir(directory)):
        file_path = os.path.join(directory, f)
        size = os.path.getsize(file_path)
        if size > limit * 1024**2:
            large_files.append(file_path)
    return large_files


def main():
    """Run script's main function."""
    parser = argparse.ArgumentParser(
        description="Checks if files in the given directory exceed the given size limit"
    )
    parser.add_argument(
        "directory", metavar="DIRECTORY", help="Directory in which to check the files"
    )
    parser.add_argument(
        "-l",
        "--limit",
        type=float,
        default="1",
        help="Size limit (in MBs) which the files should not exceed (default: %(default)s)",
    )
    args = parser.parse_args()

    large_files = find_large_files(args.directory, args.limit)
    if large_files:
        print(
            "The following files exceed the size limit of {:.2f} MB:".format(args.limit)
        )
        print("\n".join(large_files))
        exit(1)


if __name__ == "__main__":
    main()
