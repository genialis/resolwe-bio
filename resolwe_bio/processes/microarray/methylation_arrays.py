"""Process methylation array data."""
from pathlib import Path
from shutil import copy2

import pandas as pd
from plumbum import TEE

from resolwe.process import (
    Cmd,
    DataField,
    FileField,
    Persistence,
    SchedulingClass,
    StringField,
)

from resolwe_bio.process.runtime import ProcessBio


class MethylationArraySesame(ProcessBio):
    """Illumina methylation array SeSAMe process.

    Implemented SeSAMe method for analyzing methylation array from
    Illumina. For more information on the pipeline, please see
    https://www.bioconductor.org/packages/release/bioc/html/sesame.html

    This process will input IDAT methylation array files and
    produce two results. One is the quality control file with
    some basic statistics, such as mean beta, fraction of
    (un)methylated, GCT, predicted ethnicity, gender and age.
    Methylation data file holds betas, mvals and pvals for probe ids.
    """

    slug = "methylation-array-sesame"
    name = "Methylation analysis (SeSAMe)"
    process_type = "data:methylation:sesame"
    version = "1.0.0"
    category = "Methylation arrays"
    data_name = 'SeSAMe array ({{ idat_file.red_channel.file|default("?") }})'
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.CACHED
    entity = {"type": "sample"}
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {
                "image": "public.ecr.aws/s4q6j6e8/resolwebio/methylation_arrays:1.0.0"
            }
        },
        "resources": {
            "cores": 8,
            "memory": 16384,
        },
    }

    class Input:
        """Input fields to process MethylationArraySesame."""

        idat_file = DataField(
            data_type="methylationarray:idat",
            label="Illumina methylation array IDAT file",
            description="Illumina methylation array BeadChip raw IDAT file.",
        )

    class Output:
        """Output fields to process MethylationArraySesame."""

        methylation_data = FileField(label="A gzipped tab delimited file (txt.gz)")
        qc_data = FileField(label="Quality control information from SeSAMe analysis")
        species = StringField(label="Species")
        platform = StringField(label="Platform used in the analysis")

    def run(self, inputs, outputs):
        """Run MethylationArraySesame process."""

        dirdata = Path("./data")
        if not dirdata.exists():
            dirdata.mkdir()

        red = inputs.idat_file.output.red_channel.path
        green = inputs.idat_file.output.green_channel.path
        [copy2(src=x, dst=dirdata.name) for x in [red, green]]

        platform = inputs.idat_file.output.platform
        manifest = f"{platform}.hg38.manifest"

        sesame_args = [
            f"--platform={platform}",
            f"--manifest={manifest}",
        ]
        rc, _, _ = Cmd["sesame.R"][sesame_args] & TEE(retcode=None)
        # Returns QC_data.txt and beta_values_annotated.txt.gz

        if rc:
            self.error(
                "An error was encountered during the running of SeSAMe pipeline."
            )

        meth_data = "beta_values_annotated.txt.gz"
        xy = pd.read_csv(
            filepath_or_buffer=meth_data,
            sep="\t",
        )

        print((Cmd["cat"]["-e"]["QC_data.txt"])())

        mapping_filters = {"source_db": f"ILLU_{platform}"}
        mappings = self.mapping.filter(**mapping_filters)
        ensembls = [
            {"source_id": mapping.source_id, "target_id": mapping.target_id}
            for mapping in mappings
        ]

        ensembls = pd.DataFrame(ensembls)

        out = xy.merge(ensembls, left_on="probe_ids", right_on="source_id", how="left")
        out = out[["probe_ids", "betas", "mvals", "pvals", "target_id"]]
        out.columns = ["probe_ids", "betas", "mvals", "pvals", "ensembl_id"]

        out.to_csv(
            path_or_buf=meth_data,
            sep="\t",
            index=False,
        )

        outputs.qc_data = "QC_data.txt"
        outputs.methylation_data = meth_data
        outputs.species = inputs.idat_file.output.species
        outputs.platform = platform
