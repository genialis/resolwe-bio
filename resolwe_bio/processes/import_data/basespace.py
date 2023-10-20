"""Import a file from Illumina BaseSpace."""

import atexit
import gzip
import os
import time
import traceback
from pathlib import Path

from requests import RequestException, Session

from resolwe.process import (
    BooleanField,
    FileField,
    GroupField,
    IntegerField,
    Persistence,
    Process,
    SecretField,
    StringField,
)


class BaseSpaceDownloadError(Exception):
    """BaseSpace download error."""

    pass


def download_file_repeatedly(
    tries, session, file_id, file_name, expected_file_size, request_headers, error
):
    """Attempt to download BaseSpace file numerous times in case of errors."""
    for i in range(tries):
        try:
            download_file(
                session=session,
                file_id=file_id,
                file_name=file_name,
                request_headers=request_headers,
                error=error,
            )
            raise_for_file_corruption(
                file_name=file_name, expected_file_size=expected_file_size, error=error
            )
            break
        except BaseSpaceDownloadError:
            if i + 1 == tries:
                error("Could not download file from BaseSpace.")
            else:
                time.sleep(3)


def download_file(session, file_id, file_name, request_headers, error):
    """Download BaseSpace file."""
    response = make_get_request(
        session=session,
        url=get_api_file_content_url(file_id=file_id),
        headers=request_headers,
        error=error,
        stream=True,
    )

    try:
        with open(file_name, "wb") as f:
            chunk_size = 1024 * 1024 * 10
            for chunk in response.iter_content(chunk_size=chunk_size):
                f.write(chunk)
    except FileNotFoundError:
        error(f"Could not save file to {file_name}, due to directory not being found")
    except PermissionError:
        error(f"Could not save file to {file_name}, due to insufficient permissions")
    except RequestException:
        error(f"Could not save file to {file_name}, due to a network error")


def get_file_properties(session, file_id, request_headers, error):
    """Get file name and size (in bytes)."""
    response = make_get_request(
        session=session,
        url=get_api_file_url(file_id=file_id),
        headers=request_headers,
        error=error,
    )
    info = response.json()["Response"]
    return info["Name"], info["Size"]


def make_get_request(session, url, headers, error, stream=False):
    """Make a get request."""
    response = session.get(url=url, headers=headers, stream=stream, timeout=60)

    if response.status_code == 401:
        error(f"Authentication failed on URL {url}")
    elif response.status_code == 404:
        error(f"BaseSpace file {url} not found")
    elif response.status_code != 200:
        error(f"Failed to retrieve content from {url}")

    return response


def get_api_file_url(file_id):
    """Get BaseSpace API file URL."""
    api_url = "https://api.basespace.illumina.com/v1pre3"
    return f"{api_url}/files/{file_id}"


def get_api_file_content_url(file_id):
    """Get BaseSpace API file contents URL."""
    return f"{get_api_file_url(file_id=file_id)}/content"


def output(output_option, value):
    """Print to standard output."""
    if output_option == "full":
        print(value)
    elif output_option == "filename":
        if value.startswith("filename="):
            print(value[len("filename=") :])


def get_token_from_secret_file(secret_file_path, error):
    """Read secret file to obtain access token."""
    try:
        with open(secret_file_path, "r") as f:
            return f.readline()
    except FileNotFoundError:
        error("Secret file not found")
    except PermissionError:
        error("No permissions to read secret file")


def on_exit(session):
    """Clean up function called on exit."""
    session.close()


def raise_for_file_corruption(file_name, expected_file_size, error):
    """Raise an error if file does not pass integrity check."""
    # Check file size.
    actual_file_size = os.path.getsize(file_name)
    if expected_file_size != actual_file_size:
        error(
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
            error(f"File {file_name} did not pass gzip integrity check")


class BaseSpaceImport(Process):
    """Import a file from Illumina BaseSpace."""

    slug = "basespace-file-import"
    name = "BaseSpace file"
    process_type = "data:file"
    version = "1.5.1"
    category = "Import"
    data_name = 'BaseSpace ({{ file_id|default("?") }})'
    persistence = Persistence.TEMP
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/common:4.1.1"}
        },
        "resources": {
            "cores": 1,
            "memory": 1024,
            "network": True,
            "secrets": True,
        },
    }

    class Input:
        """Input fields to process BaseSpaceImport."""

        file_id = StringField(label="BaseSpace file ID")
        access_token_secret = SecretField(
            label="BaseSpace access token",
            description="BaseSpace access token secret handle needed to download the file.",
        )

        class Advanced:
            """Advanced options."""

            output = StringField(
                label="Output",
                allow_custom_choice=False,
                choices=[("full", "Full"), ("filename", "Filename")],
                default="filename",
                description="Sets what is printed to standard output. "
                "Argument 'Full' outputs everything, "
                "argument 'Filename' outputs only file names of downloaded files.",
            )
            tries = IntegerField(
                label="Tries",
                description="Number of tries to download a file before giving up.",
                range=[1, 10],
                default=3,
            )
            verbose = BooleanField(
                label="Verbose",
                default=False,
                description="Print detailed exception information to standard output "
                "when error occurs. Output argument had no effect on this argument.",
            )

        advanced = GroupField(Advanced, label="Advanced options")

    class Output:
        """Output fields to process BaseSpaceImport."""

        file = FileField(label="File with reads")

    def run(self, inputs, outputs):
        """Run import."""

        secret_path = Path("/secrets") / inputs.access_token_secret["handle"]

        session = Session()
        atexit.register(on_exit, session)

        try:
            file_id = inputs.file_id
            access_token = get_token_from_secret_file(
                secret_file_path=secret_path, error=self.error
            )
            headers = {"x-access-token": access_token}

            file_name, file_size = get_file_properties(
                session=session,
                file_id=file_id,
                request_headers=headers,
                error=self.error,
            )
            download_file_repeatedly(
                tries=inputs.advanced.tries,
                session=session,
                file_id=file_id,
                file_name=file_name,
                expected_file_size=file_size,
                request_headers=headers,
                error=self.error,
            )
            output(inputs.advanced.output, f"filename={file_name}")
        except Exception as error:
            if inputs.advanced.verbose:
                traceback.print_exc()
                self.error(
                    "Unexpected error occurred while trying to download files from BaseSpace. "
                    "Check standard output for more details."
                )
            else:
                print(str(error))
                self.error(
                    "Unexpected error occurred while trying to download files from BaseSpace. "
                    "Set Verbose to True to see the traceback."
                )

        outputs.file = file_name
