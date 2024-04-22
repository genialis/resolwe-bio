"""GEO import pipeline."""

import re
from pathlib import Path

import GEOparse
import pandas as pd
import requests

from resolwe.process import (
    BooleanField,
    FileField,
    GroupField,
    IntegerField,
    Process,
    SchedulingClass,
    StringField,
)
from resolwe.process.models import Data


def parse_sample(gse, db_accession, gse_name):
    """Parse sample information from GEO."""
    sample = {"Database accession": db_accession}
    for k, v in gse.gsms[gse_name].metadata.items():
        if len(v) == 1:
            sample[k] = v[0]
        else:
            if all(": " in substring for substring in v):
                for meta in v:
                    key, value = meta.split(": ", 1)
                    sample[key] = value
            else:
                sample[k] = " ".join(v)
    return sample


def create_metadata(gse, run_info):
    """Create a tab-separated metadata file."""
    collection = [
        parse_sample(gse, row["Accession"], row["SampleName"])
        for _, row in run_info.iterrows()
    ]
    metadata = pd.json_normalize(collection)
    metadata.insert(0, "Sample name", metadata["title"])
    return metadata.set_index(["Sample name"], drop=False)


class GeoImport(Process):
    """Import all runs from a GEO Series.

    WARNING: Additional costs for storage and processing may be incurred
    if a very large data set is selected.

    RNA-seq ChIP-Seq, ATAC-Seq and expression microarray datasets can be
    uploaded.

    For RNA-Seq data sets this runs the SRA import process for each
    experiment (SRX) from the selected RNA-Seq GEO Series. The same
    procedure is followed for ChIP-Seq and ATAC-Seq data sets.

    If GSE contains microarray data, it downloads individual samples and
    uploads them as microarray expression objects. Probe IDs can be
    mapped to the Ensembl IDs if the corresponding GPL platform is
    supported, otherwise, a custom mapping file should be provided.
    Currently supported platforms are: GPL74, GPL201, GPL96, GPL571,
    GPL97, GPL570, GPL91, GPL8300, GPL92, GPL93, GPL94, GPL95, GPL17586,
    GPL5175, GPL80, GPL6244, GPL16686, GPL15207, GPL1352, GPL11068,
    GPL26966, GPL6848, GPL14550, GPL17077, GPL16981, GPL13497, GPL6947,
    GPL10558, GPL6883, GPL13376,GPL6884, GPL6254.

    In addition metadata table with sample information is created and
    uploaded to the same collection.
    """

    slug = "geo-import"
    name = "GEO import"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {
                "image": "public.ecr.aws/genialis/resolwebio/common:4.1.1",
            },
        },
        "resources": {
            "cores": 1,
            "memory": 16384,
            "network": True,
        },
    }
    data_name = "{{ gse_accession }}"
    version = "2.8.0"
    process_type = "data:geo"
    category = "Import"
    scheduling_class = SchedulingClass.BATCH

    class Input:
        """Input fields."""

        gse_accession = StringField(
            label="GEO accession", description="Enter a GEO series accession number."
        )

        class Advanced:
            """Advanced options."""

            prefetch = BooleanField(label="Prefetch SRA file", default=True)
            max_size_prefetch = StringField(
                label="Maximum file size to download in KB",
                default="20G",
                description="A unit prefix can be used instead of a value in KB (e.g. 1024M or 1G).",
            )
            min_spot_id = IntegerField(label="Minimum spot ID", required=False)
            max_spot_id = IntegerField(label="Maximum spot ID", required=False)
            min_read_len = IntegerField(label="Minimum read length", required=False)
            clip = BooleanField(label="Clip adapter sequences", default=False)
            aligned = BooleanField(label="Dump only aligned sequences", default=False)
            unaligned = BooleanField(
                label="Dump only unaligned sequences", default=False
            )
            mapping_file = FileField(
                label="File with probe ID mappings",
                description="The file should be tab-separated and contain two columns with their column names. The "
                "first column should contain Gene IDs and the second one should contain probe names. Supported file "
                "extensions are .tab.*, .tsv.*, .txt.*",
                required=False,
            )
            source = StringField(
                label="Gene ID source",
                description="Gene ID source used for probe mapping is required when using a custom file.",
                allow_custom_choice=True,
                required=False,
                choices=[
                    ("AFFY", "AFFY"),
                    ("DICTYBASE", "DICTYBASE"),
                    ("ENSEMBL", "ENSEMBL"),
                    ("NCBI", "NCBI"),
                    ("UCSC", "UCSC"),
                ],
            )
            build = StringField(
                label="Genome build",
                description="Genome build of mapping file is required when using a custom file.",
                required=False,
            )

        advanced = GroupField(Advanced, label="Advanced options")

    def upload_rna_gse(self, inputs, gse):
        """Upload RNA samples from GEO series.

        Find SRX accessions on a GEO sample (GSM) and fetch the
        corresponding Run Info from SRA. Use run info to retrieve
        individual run accessions (SRR) and library layouts needed for
        sra-import. Samples are renamed to their SRA experiment
        accessions (SRX).
        """
        process_inputs = {
            "sra_accession": [],
            "advanced": {
                "prefetch": inputs.advanced.prefetch,
                "max_size_prefetch": inputs.advanced.max_size_prefetch,
                "clip": inputs.advanced.clip,
                "aligned": inputs.advanced.aligned,
                "unaligned": inputs.advanced.unaligned,
                "min_spot_id": inputs.advanced.min_spot_id,
                "max_spot_id": inputs.advanced.max_spot_id,
                "min_read_len": inputs.advanced.min_read_len,
            },
        }

        sample_info = {}
        for name, gsm in gse.gsms.items():
            if "SRA" in gse.gsms[name].relations:
                # Match NCBI's SRX and EBI's ERX sample accessions.
                sample_found = re.findall(
                    r"([SE]RX\d{6,8})", str(gse.gsms[name].relations["SRA"])
                )
            else:
                sample_found = None

            if sample_found:
                for srx_id in sample_found:
                    sample_info[srx_id] = name
                    info_file = f"{gse.name}.csv"
                    run_info = requests.get(
                        url="https://eutils.ncbi.nlm.nih.gov/Traces/sra/sra.cgi",
                        params={
                            "save": "efetch",
                            "db": "sra",
                            "rettype": "runinfo",
                            "term": srx_id,
                        },
                    )

                    if run_info.status_code != 200:
                        self.error(
                            f"Failed to fetch SRA runs for project {srx_id} belonging to {gse.name}."
                        )

                    elif run_info.text.isspace():
                        self.error(
                            f"Got an empty response from SRA for SRX ID {srx_id} belonging to {gse.name}."
                        )

                    else:
                        with open(info_file, "wb") as handle:
                            handle.write(run_info.content)

                    run_info = pd.read_csv(
                        info_file, usecols=["Run", "SampleName", "LibraryLayout"]
                    )
                    run_info = run_info.set_index("Run", drop=False)

                    process_inputs["sra_accession"] = run_info.index.values.tolist()
                    assert run_info.nunique().loc["LibraryLayout"] == 1
                    lib_type = run_info["LibraryLayout"].iloc[0]

                    if lib_type == "PAIRED":
                        self.run_process("import-sra-paired", process_inputs)
                    elif lib_type == "SINGLE":
                        self.run_process("import-sra-single", process_inputs)
                    else:
                        self.error(
                            f"Unsupported library layout expected SINGLE or PAIRED but got {lib_type}."
                        )

                    entity_name = process_inputs["sra_accession"][0]
                    sra_data = Data.filter(entity__name=entity_name)[-1]
                    sra_data.entity.name = gsm.metadata["title"][0]
            else:
                self.warning(
                    f"Matching SRX accession number for {gsm.metadata['title'][0]}({name} "
                    "was not found in GEO metadata."
                )

        return pd.DataFrame(
            sample_info.items(), columns=["Accession", "SampleName"]
        ).set_index("Accession", drop=False)

    def upload_ma_gse(self, inputs, gse):
        """Upload microarray samples from a GEO series."""
        unmapped_data = []
        for name, gsm in gse.gsms.items():
            exp_path = f"{name}.tsv"
            gsm.table.to_csv(exp_path, sep="\t", index=False)

            platform_id = gsm.metadata["platform_id"][0]
            process_inputs = {
                "exp": exp_path,
                "exp_type": gsm.columns.loc["VALUE", "description"],
                "platform": gse.gpls[platform_id].metadata["title"][0],
                "platform_id": platform_id,
                "species": gsm.metadata["organism_ch1"][0],
            }

            self.run_process("upload-microarray-expression", process_inputs)
            ma_data = Data.filter(entity__name=exp_path)[-1]
            ma_data.entity.name = gsm.metadata["title"][0]
            unmapped_data.append(ma_data.id)

        names = [gsm for gsm in gse.gsms.keys()]
        info_df = {"Accession": names, "SampleName": names}

        mapping_inputs = {"expressions": unmapped_data}
        if inputs.advanced.mapping_file:
            mapping_file = inputs.advanced.mapping_file.import_file(
                imported_format="compressed"
            )

            stem = Path(mapping_file).stem
            supported_extensions = (".tab", ".tsv", ".txt")
            if not stem.endswith(supported_extensions):
                self.error(
                    "Mapping file has unsupported file name extension. "
                    f"The supported extensions are {supported_extensions}."
                )

            mapping_inputs["mapping_file"] = mapping_file
            if inputs.advanced.source:
                mapping_inputs["source"] = inputs.advanced.source
            else:
                self.error(
                    "Custom probe id mapping file was provided but no source was selected."
                )
            if inputs.advanced.build:
                mapping_inputs["build"] = inputs.advanced.build
            else:
                self.error(
                    "Custom probe id mapping file was provided but genome build was not defined."
                )

        self.run_process("map-microarray-probes", mapping_inputs)
        return pd.DataFrame(info_df).set_index("Accession", drop=False)

    def run(self, inputs, outputs):
        """Run the analysis."""

        if not re.match(r"(GSE\d{1,8})", inputs.gse_accession):
            self.error(
                f"GEO series accessions (GSE) are supported but {inputs.gse_accession} was provided."
            )

        try:
            gse = GEOparse.get_GEO(geo=inputs.gse_accession, destdir="./")
        except IOError:
            self.error(
                f"Download of {inputs.gse_accession} failed. ID could be incorrect or the data might not be "
                "public yet."
            )
        except Exception as err:
            self.error(
                f"Download of {inputs.gse_accession} failed. GEO parse failed with {err}"
            )

        supported = [
            "Expression profiling by high throughput sequencing",
            "Expression profiling by array",
            "Genome binding/occupancy profiling by high throughput sequencing",
        ]

        gse_type = gse.get_type() if type(gse.get_type()) is list else [gse.get_type()]
        if set(gse_type).intersection(set(supported)):
            if "SuperSeries of" in gse.relations:
                # This is a mixed GSE series which needs to be unpacked.
                super_series = [
                    GEOparse.get_GEO(geo=accession, destdir="./")
                    for accession in gse.relations["SuperSeries of"]
                ]
            else:
                super_series = [gse]
        else:
            self.error(
                f"No supported series types found. Got {', '.join(gse_type)} but only {', '.join(supported)} "
                "are supported."
            )

        metadata_tables = {}
        for series in super_series:
            series_type = series.get_type()
            if series_type == "Expression profiling by high throughput sequencing":
                run_info = self.upload_rna_gse(inputs, series)
                if run_info.empty:
                    self.warning(
                        f"No samples with SRA files were found for GEO series {series.name}."
                    )
                    continue

                metadata_tables[series.name] = create_metadata(series, run_info)
            elif series_type == "Expression profiling by array":
                run_info = self.upload_ma_gse(inputs, series)
                metadata_tables[series.name] = create_metadata(series, run_info)
            elif (
                series_type
                == "Genome binding/occupancy profiling by high throughput sequencing"
            ):
                run_info = self.upload_rna_gse(inputs, series)
                if run_info.empty:
                    self.warning(
                        f"No samples with SRA files were found for GEO series {series.name}."
                    )
                    continue

                metadata_tables[series.name] = create_metadata(series, run_info)
            else:
                self.warning(
                    f"The upload of {series_type} is currently not supported. Samples from {series.name} will be "
                    "skipped."
                )

        if len(metadata_tables) == 0:
            self.error("No supported data found.")

        meta_file = f"{inputs.gse_accession}_metadata.tsv"
        metadata = pd.concat(metadata_tables.values(), join="outer", ignore_index=False)
        metadata.to_csv(meta_file, sep="\t", index=False)
        self.run_process("upload-metadata-unique", {"src": meta_file})
