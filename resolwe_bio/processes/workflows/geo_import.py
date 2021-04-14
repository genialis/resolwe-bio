"""GEO import pipeline."""
import re

import GEOparse
import pandas as pd
import requests

from resolwe.process import (
    BooleanField,
    GroupField,
    IntegerField,
    Process,
    SchedulingClass,
    StringField,
)
from resolwe.process.models import Data


def parse_sample(gse, sample_name, gse_name):
    """Parse sample information from GEO."""
    sample = {"EntityName": sample_name}
    for k, v in gse.gsms[gse_name].metadata.items():
        if len(v) == 1:
            sample[k] = v[0]
        else:
            if all(": " in substring for substring in v):
                for meta in v:
                    key, value = meta.split(": ")
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
    metadata = pd.json_normalize(collection).set_index(["EntityName"], drop=False)
    return metadata


def construct_descriptor(metadata, sample_name):
    """Construct a descriptor from sample metadata.

    Dictionary with GEO metadata that matches the sample descriptor
    schema is created. Atributes under general that have no
    predetermined choices are matched with our naming if they
    exist in the metadata. Other fields with choices and the
    experimental section are filled separately.
    """

    sample_metadata = metadata.loc[sample_name]
    descriptor = {"general": {}, "experiment": {}}
    # TODO Replace fixed values with calls to descriptor schema once available.
    # Also organ / tissue should be added then.
    species = [
        "Caenorhabditis elegans",
        "Cricetulus griseus",
        "Dictyostelium discoideum",
        "Dictyostelium purpureum",
        "Drosophila melanogaster",
        "Homo sapiens",
        "Macaca mulatta",
        "Mus musculus",
        "Odocoileus virginianus texanus",
        "Rattus norvegicus",
        "Solanum tuberosum",
    ]
    molecule_choices = [
        "total_rna",
        "polya_rna",
        "cytoplasmic_rna",
        "nuclear_rna",
        "genomic_dna",
        "protein",
        "other",
    ]
    assay_types = [
        "rna-seq",
        "other",
    ]
    platform_types = [
        "nextseq_500",
        "hiseq_2500",
        "hiseq_2000",
        "novaseq_6000",
        "other",
    ]
    general_attributes = {
        "description": "description",
        "cell type": "cell_type",
        "source_name_ch1": "biosample_source",
        "growth_protocol_ch1": "growth_protocol",
        "treatment_protocol_ch1": "treatment_protocol",
    }
    if (
        "organism_ch1" in metadata.columns
        and sample_metadata["organism_ch1"] in species
    ):
        descriptor["general"]["species"] = sample_metadata["organism_ch1"]

    if "cell line" in metadata.columns:
        descriptor["general"]["biosample_type"] = "cell_line"
        descriptor["general"]["cell_line"] = sample_metadata["cell line"]
    elif "tissue" in metadata.columns:
        descriptor["general"]["biosample_type"] = "tissue"
    if "contact_name" in metadata.columns:
        descriptor["general"]["annotator"] = sample_metadata["contact_name"].replace(
            ",,", " "
        )

    for geo_attribute, attribute in general_attributes.items():
        if geo_attribute in metadata.columns:
            descriptor["general"][attribute] = sample_metadata[geo_attribute]

    if "library_strategy" in metadata.columns:
        formated_assay = sample_metadata["library_strategy"].lower().replace(" ", "-")
        if formated_assay in assay_types:
            descriptor["experiment"]["assay_type"] = formated_assay
    if "extract_protocol_ch1" in metadata.columns:
        descriptor["experiment"]["extract_protocol"] = sample_metadata[
            "extract_protocol_ch1"
        ]
    if "molecule_ch1" in metadata.columns:
        formated_molecule = sample_metadata["molecule_ch1"].lower().replace(" ", "_")
        if formated_molecule in molecule_choices:
            descriptor["experiment"]["molecule"] = formated_molecule
    if "instrument_model" in metadata.columns:
        formated_platform = (
            sample_metadata["instrument_model"]
            .replace("Illumina ", "")
            .lower()
            .replace(" ", "_")
        )
        if formated_platform in platform_types:
            descriptor["experiment"]["platform"] = formated_platform
    return descriptor


class GeoImport(Process):
    """Import all runs from a GEO Series.

    WARNING: Additional costs for storage and processing may be incurred
    if a very large data set is selected.

    This runs SRA import process for each individual experiment (SRX)
    from the selected GEO Series. In addition metadata table with sample
    information is created and uploadad to the same collection.
    Currently only RNA-seq datasets are supported.
    """

    slug = "geo-import"
    name = "GEO import"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {
                "image": "public.ecr.aws/s4q6j6e8/resolwebio/common:2.7.0",
            },
        },
        "resources": {
            "cores": 1,
            "memory": 8192,
            "network": True,
        },
    }
    data_name = "{{ gse_accession }}"
    version = "1.0.1"
    process_type = "data:geo"
    category = "Import"
    scheduling_class = SchedulingClass.BATCH

    class Input:
        """Input fields."""

        gse_accession = StringField(
            label="GEO accession", description="Enter a GEO series accession number."
        )
        show_advanced = BooleanField(label="Show advanced options", default=False)

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

        advanced = GroupField(
            Advanced, label="Advanced options", hidden="!show_advanced"
        )

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
            },
        }

        if inputs.advanced.min_spot_id:
            process_inputs["advanced"]["min_spot_id"] = inputs.advanced.min_spot_id
        if inputs.advanced.max_spot_id:
            process_inputs["advanced"]["max_spot_id"] = inputs.advanced.max_spot_id
        if inputs.advanced.min_read_len:
            process_inputs["advanced"]["min_read_len"] = inputs.advanced.min_read_len
        if inputs.advanced.min_read_len:
            process_inputs["advanced"]["min_read_len"] = inputs.advanced.min_read_len

        sample_info = {}
        for name, gsm in gse.gsms.items():
            sample_found = re.findall(
                r"(SRX\d{6,8})", str(gse.gsms[name].relations["SRA"])
            )
            if sample_found:
                for srx_id in sample_found:
                    sample_info[srx_id] = name
                    sra_url = f"http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term={srx_id}"
                    info_file = f"{gse.name}.csv"
                    run_info = requests.get(sra_url).content
                    if run_info == b"\n":
                        self.error(
                            f"Failed to fetch SRA runs for project {srx_id} belonging to {gse.name}."
                        )
                    else:
                        with open(info_file, "wb") as handle:
                            handle.write(requests.get(sra_url).content)
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
                    sra_data.entity.name = srx_id
            else:
                self.error(
                    f"Matching SRX accession number for {name} was not found in GEO metadata."
                )

        return pd.DataFrame(
            sample_info.items(), columns=["Accession", "SampleName"]
        ).set_index("Accession", drop=False)

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

        supported = set(
            [
                "Expression profiling by high throughput sequencing",
            ]
        )
        gse_type = gse.get_type() if type(gse.get_type()) is list else [gse.get_type()]
        if set(gse_type).intersection(supported):
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
                f"No supported series types found. Got {', '.join(gse_type)} but only {' and '.join(supported)} "
                "are supported."
            )

        metadata_tables = {}
        for series in super_series:
            series_type = series.get_type()
            if series_type == "Expression profiling by high throughput sequencing":
                run_info = self.upload_rna_gse(inputs, series)
                metadata_tables[series.name] = create_metadata(series, run_info)
            else:
                self.warning(
                    f"The upload of {series_type} is currently not supported. Samples from {series.name} will be "
                    "skipped."
                )
        meta_file = f"{inputs.gse_accession}_metadata.tsv"
        metadata = pd.concat(metadata_tables.values(), join="outer", ignore_index=False)
        metadata.to_csv(meta_file, sep="\t", index=False)
        self.run_process("upload-orange-metadata", {"src": meta_file})

        for entity_name in metadata["EntityName"].values:
            objects = Data.filter(entity__name=entity_name)
            if len(objects) > 1:
                self.warning(
                    f"Multiple samples with entity name {entity_name} are present, descriptor will be added only "
                    "to the last one"
                )
            obj = objects[-1]
            obj.entity.descriptor = construct_descriptor(metadata, obj.entity_name)
