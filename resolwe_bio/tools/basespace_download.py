#!/usr/bin/env python3
"""Tool to download files from BaseSpace."""

import sys
import traceback
import atexit
import argparse
import requests


class BaseSpaceDownloadError(Exception):
    """BaseSpace download error."""

    pass


def main():
    """Entry point."""
    session = requests.Session()
    atexit.register(on_exit, session)

    parser = argparse.ArgumentParser(description='Download file from Illumina BaseSpace.')
    parser.add_argument('--file-id',
                        dest='file_ids',
                        action='append',
                        required=True,
                        help="BaseSpace file ID. This argument can be repeated to specify multiple files.")
    parser.add_argument('--access-token-secret-path',
                        dest='access_token_secret_path',
                        required=True,
                        help="BaseSpace access token secret path.")
    parser.add_argument('--output',
                        dest='output',
                        type=str,
                        choices=['full', 'filename'],
                        default='full',
                        help="Sets what is printed to standard output. "
                             "Argument 'full' outputs everything, "
                             "argument 'filename' outputs only file names of downloaded files.")
    parser.add_argument('--verbose',
                        dest='verbose',
                        action='store_true',
                        default=False,
                        help="Print detailed exception information when error occurs. "
                             "Output argument had no effect on this argument.")
    args = parser.parse_args()

    try:
        file_ids = args.file_ids
        access_token = get_token_from_secret_file(args.access_token_secret_path)
        headers = {'x-access-token': access_token}

        for file_id in file_ids:
            file_name = get_file_name(session, file_id, headers)
            download_file(session, file_id, file_name, headers)
            output(args.output, 'filename={}'.format(file_name))
    except:
        if args.verbose:
            traceback.print_exc()
        else:
            print("An error occurred while processing the Basespace download request. Use --verbose to see details.")

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
    if output_option == 'full':
        print(value)
    elif output_option == 'filename':
        if value.startswith('filename='):
            print(value[len('filename='):])
    else:
        print("Internal error: output argument {} handling not implemented".format(output_option))

        sys.exit(1)


def get_token_from_secret_file(secret_file_path):
    """Read secret file to obtain access token."""
    try:
        with open(secret_file_path, 'r') as f:
            return f.readline()
    except FileNotFoundError:
        raise BaseSpaceDownloadError('Secret file not found')
    except PermissionError:
        raise BaseSpaceDownloadError('No permissions to read secret file')


def make_get_request(session, url, headers, stream=False):
    """Make a get request."""
    response = session.get(url, headers=headers, stream=stream)

    if response.status_code == 401:
        raise BaseSpaceDownloadError('Authentication failed on URL {}'.format(url))
    elif response.status_code == 404:
        raise BaseSpaceDownloadError('BaseSpace file {} not found'.format(url))

    return response


def get_basespace_api_url():
    """Get base BaseSpace API URL."""
    return 'https://api.basespace.illumina.com/v1pre3'


def get_basespace_api_file_url(file_id):
    """Get BaseSpace API file URL."""
    return '{}/files/{}'.format(get_basespace_api_url(), file_id)


def get_basespace_api_file_content_url(file_id):
    """Get BaseSpace API file contents URL."""
    return '{}/content'.format(get_basespace_api_file_url(file_id))


def get_file_name(session, file_id, request_headers):
    """Get file name."""
    response = make_get_request(session, get_basespace_api_file_url(file_id), request_headers)
    return response.json()['Response']['Name']


def download_file(session, file_id, file_name, request_headers):
    """Download BaseSpace file."""
    response = make_get_request(session, get_basespace_api_file_content_url(file_id), request_headers, True)

    try:
        with open(file_name, 'wb') as f:
            for chunk in response.iter_content(chunk_size=1024):
                f.write(chunk)
    except FileNotFoundError:
        raise BaseSpaceDownloadError('Could not save file to {}, due to directory not being found'.format(file_name))
    except PermissionError:
        raise BaseSpaceDownloadError('Could not save file to {}, due to insufficient permissions'.format(file_name))


if __name__ == "__main__":
    main()
