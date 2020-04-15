#!/usr/bin/env python3
"""Tool to download files from BaseSpace."""

import argparse
import atexit
import gzip
import os
import sys
import time
import traceback

import requests


class BaseSpaceDownloadError(Exception):
    """BaseSpace download error."""

    pass


def main():
    """Entry point."""
    session = requests.Session()
    atexit.register(on_exit, session)

    parser = argparse.ArgumentParser(
        description="Download file from Illumina BaseSpace."
    )

    def check_tries_within_range(arg):
        try:
            value = int(arg)
        except ValueError:
            raise argparse.ArgumentTypeError("Not an integer")

        if not 1 <= value <= 10:
            raise argparse.ArgumentTypeError("Not between 1 and 10")

        return value

    parser.add_argument(
        "--file-id",
        dest="file_ids",
        action="append",
        required=True,
        help="BaseSpace file ID. This argument can be repeated to specify multiple files.",
    )
    parser.add_argument(
        "--access-token-secret-path",
        dest="access_token_secret_path",
        required=True,
        help="BaseSpace access token secret path.",
    )
    parser.add_argument(
        "--output",
        dest="output",
        type=str,
        choices=["full", "filename"],
        default="full",
        help="Sets what is printed to standard output. "
        "Argument 'full' outputs everything, "
        "argument 'filename' outputs only file names of downloaded files.",
    )
    parser.add_argument(
        "--tries",
        dest="tries",
        type=check_tries_within_range,
        default=3,
        help="Number of tries to download a file before giving up.",
    )
    parser.add_argument(
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Print detailed exception information when error occurs. "
        "Output argument had no effect on this argument.",
    )
    args = parser.parse_args()

    try:
        file_ids = args.file_ids
        access_token = get_token_from_secret_file(args.access_token_secret_path)
        headers = {"x-access-token": access_token}

        for file_id in file_ids:
            file_name, file_size = get_file_properties(session, file_id, headers)
            download_file_repeatedly(
                args.tries, session, file_id, file_name, file_size, headers
            )
            output(args.output, f"filename={file_name}")
    except Exception as error:
        if args.verbose:
            traceback.print_exc()
        else:
            print(str(error))
            print("Use --verbose to see traceback.")

        sys.exit(1)


def on_exit(session):
    """Clean up function called on exit."""
    session.close()


def output(output_option, value):
    """
    Print to standard output.

    This function should be used instead of ``print`` function. Printing errors is exempted
    and can be printed without using this function.

    """
    if output_option == "full":
        print(value)
    elif output_option == "filename":
        if value.startswith("filename="):
            print(value[len("filename=") :])
    else:
        print(
            f"Internal error: output argument {output_option} handling not implemented"
        )

        sys.exit(1)


def get_token_from_secret_file(secret_file_path):
    """Read secret file to obtain access token."""
    try:
        with open(secret_file_path, "r") as f:
            return f.readline()
    except FileNotFoundError:
        raise BaseSpaceDownloadError("Secret file not found")
    except PermissionError:
        raise BaseSpaceDownloadError("No permissions to read secret file")


def make_get_request(session, url, headers, stream=False):
    """Make a get request."""
    response = session.get(url, headers=headers, stream=stream, timeout=60)

    if response.status_code == 401:
        raise BaseSpaceDownloadError(f"Authentication failed on URL {url}")
    elif response.status_code == 404:
        raise BaseSpaceDownloadError(f"BaseSpace file {url} not found")
    elif response.status_code != 200:
        raise BaseSpaceDownloadError(f"Failed to retrieve content from {url}")

    return response


def get_api_url():
    """Get base BaseSpace API URL."""
    return "https://api.basespace.illumina.com/v1pre3"


def get_api_file_url(file_id):
    """Get BaseSpace API file URL."""
    api_url = get_api_url()
    return f"{api_url}/files/{file_id}"


def get_api_file_content_url(file_id):
    """Get BaseSpace API file contents URL."""
    return f"{get_api_file_url(file_id)}/content"


def get_file_properties(session, file_id, request_headers):
    """Get file name and size (in bytes)."""
    response = make_get_request(session, get_api_file_url(file_id), request_headers)
    info = response.json()["Response"]
    return info["Name"], info["Size"]


def raise_for_file_corruption(file_name, expected_file_size):
    """Raise an error if file does not pass integrity check."""
    # Check file size.
    actual_file_size = os.path.getsize(file_name)
    if expected_file_size != actual_file_size:
        raise BaseSpaceDownloadError(
            f"File's ({file_name}) expected size ({expected_file_size}) "
            f"does not match its actual size ({actual_file_size})"
        )

    # Check gzip integrity.
    if file_name.split(".")[-1] == "gz":
        try:
            with gzip.open(file_name, "rb") as f:
                chunk_size = 1024 * 1024 * 10
                while bool(f.read(chunk_size)):
                    pass
        except OSError:
            raise BaseSpaceDownloadError(
                f"File {file_name} did not pass gzip integrity check"
            )


def download_file_repeatedly(
    tries, session, file_id, file_name, expected_file_size, request_headers
):
    """Attempt to download BaseSpace file numerous times in case of errors."""
    for index in range(tries):
        try:
            download_file(session, file_id, file_name, request_headers)
            raise_for_file_corruption(file_name, expected_file_size)
            break
        except BaseSpaceDownloadError as error:
            if index + 1 >= tries:
                raise error

        time.sleep(3)


def download_file(session, file_id, file_name, request_headers):
    """Download BaseSpace file."""
    response = make_get_request(
        session, get_api_file_content_url(file_id), request_headers, stream=True
    )

    try:
        with open(file_name, "wb") as f:
            chunk_size = 1024 * 1024 * 10
            for chunk in response.iter_content(chunk_size=chunk_size):
                f.write(chunk)
    except FileNotFoundError:
        raise BaseSpaceDownloadError(
            f"Could not save file to {file_name}, due to directory not being found"
        )
    except PermissionError:
        raise BaseSpaceDownloadError(
            f"Could not save file to {file_name}, due to insufficient permissions"
        )
    except requests.RequestException as error:
        raise BaseSpaceDownloadError(
            f"Could not save file to {file_name}, due to a network error"
        ) from error


if __name__ == "__main__":
    main()
