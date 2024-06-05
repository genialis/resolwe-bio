"""Map microarray probe IDs."""

import re
from io import StringIO
from pathlib import Path

import pandas as pd
from bioservices import BioMart

from resolwe.process import (
    DataField,
    FileField,
    ListField,
    Process,
    SchedulingClass,
    StringField,
)

PLATFORM_MAP = {
    "GPL74": "affy_hc_g110",
    "GPL201": "affy_hg_focus",
    "GPL96": "affy_hg_u133a",
    "GPL571": "affy_hg_u133a_2",
    "GPL97": "affy_hg_u133b",
    "GPL570": "affy_hg_u133_plus_2",
    "GPL91": "affy_hg_u95a",
    "GPL8300": "affy_hg_u95av2",
    "GPL92": "affy_hg_u95b",
    "GPL93": "affy_hg_u95c",
    "GPL94": "affy_hg_u95d",
    "GPL95": "affy_hg_u95e",
    "GPL17586": "affy_hta_2_0",
    "GPL5175": "affy_huex_1_0_st_v2",
    "GPL80": "affy_hugenefl",
    "GPL6244": "affy_hugene_1_0_st_v1",
    "GPL16686": "affy_hugene_2_0_st_v1",
    "GPL15207": "affy_primeview",
    "GPL1352": "affy_u133_x3p",
    "GPL11068": "agilent_cgh_44b",
    "GPL26966": "agilent_gpl26966",
    "GPL6848": "agilent_gpl6848",
    "GPL14550": "agilent_sureprint_g3_ge_8x60k",
    "GPL17077": "agilent_sureprint_g3_ge_8x60k_v2",
    "GPL16981": "agilent_wholegenome",
    "GPL6480": "agilent_wholegenome_4x44k_v1",
    "GPL13497": "agilent_wholegenome_4x44k_v2",
    "GPL6947": "illumina_humanht_12_v3",
    "GPL10558": "illumina_humanht_12_v4",
    "GPL6883": "illumina_humanref_8_v3",
    "GPL13376": "illumina_humanwg_6_v2",
    "GPL6884": "illumina_humanwg_6_v3",
    "GPL6254": "phalanx_onearray",
}


def get_exp_table(fname, sample_name):
    """Get a formated expression table."""
    table = pd.read_csv(fname, sep="\t")
    if "ID_REF" in table.columns:
        table = table.set_index("ID_REF")
    else:
        table = table.set_index(table.columns[0])
    if "VALUE" in table.columns:
        table = table["VALUE"].rename(sample_name)
    else:
        table = table.iloc[:, 0].rename(sample_name)
    return table


def join_expressions(expressions):
    """Join expression tables of data objects in a list."""
    tables = [get_exp_table(e.output.exp.path, e.entity.name) for e in expressions]
    data = pd.concat(tables, axis=1, join="outer")
    data.index = data.index.astype(str)
    data.index.name = "probe"
    return data


class MapMicroarrayProbes(Process):
    """Map microarray probes to Gene IDs.

    Mapping can be done automatically or using a custom mapping file.
    For automatic probe mapping all 'Normalized expression' objects
    should have a GEO platform ID. If the platform is supported the
    provided probe IDs will be mapped to the corresponding Ensembl IDs.
    Currently supported platforms are: GPL74, GPL201, GPL96, GPL571,
    GPL97, GPL570, GPL91, GPL8300, GPL92, GPL93, GPL94, GPL95, GPL17586,
    GPL5175, GPL80, GPL6244, GPL16686, GPL15207, GPL1352, GPL11068,
    GPL26966, GPL6848, GPL14550, GPL17077, GPL16981, GPL13497, GPL6947,
    GPL10558, GPL6883, GPL13376,GPL6884, GPL6254.
    """

    slug = "map-microarray-probes"
    name = "Map microarray probes"
    process_type = "data:microarray:mapping"
    version = "1.1.1"
    scheduling_class = SchedulingClass.BATCH
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/common:4.1.1"}
        },
    }
    data_name = "Probe mapping"

    class Input:
        """Input fields to process MapMicroarrayProbes."""

        expressions = ListField(
            DataField("microarray:normalized"),
            label="Normalized expressions",
        )
        mapping_file = FileField(
            label="File with probe ID mappings",
            description="The file should be tab-separated and contain two columns with their column names. The first "
            "column should contain Gene IDs and the second one should contain probe names. Supported file extensions "
            "are .tab.*, .tsv.*, .txt.*",
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

    class Output:
        """Output fields to process MapMicroarrayProbes."""

        mapped_exp = FileField(label="Mapped expressions")
        probe_mapping = StringField(label="Probe to transcript mapping used")
        mapping = FileField(label="Mapping file")
        platform = StringField(label="Microarray platform type")
        platform_id = StringField(label="GEO platform ID", required=False)

    def run(self, inputs, outputs):
        """Run the analysis."""

        for exp in inputs.expressions:
            if exp.output.species != inputs.expressions[0].output.species:
                self.error(
                    "Input samples are of different Species: "
                    f"{exp.output.species} and {inputs.expressions[0].output.species}."
                )
            if exp.output.exp_type != inputs.expressions[0].output.exp_type:
                self.error(
                    "Input samples have different Normalization types: "
                    f"{exp.output.exp_type} and {inputs.expressions[0].output.exp_type}."
                )
            if exp.output.platform != inputs.expressions[0].output.platform:
                self.error(
                    "Input samples have different Microarray platform types: "
                    f"{exp.output.platform} and {inputs.expressions[0].output.platform}."
                )
            if exp.output.platform_id != inputs.expressions[0].output.platform_id:
                self.error(
                    "Input samples have different GEO platform IDs: "
                    f"{exp.output.platform_id} and {inputs.expressions[0].output.platform_id}."
                )

        species = inputs.expressions[0].output.species
        platform = inputs.expressions[0].output.platform
        platform_id = inputs.expressions[0].output.platform_id

        joined_expressions = join_expressions(inputs.expressions)
        probe_ids = joined_expressions.index.unique()

        if inputs.mapping_file:
            mapping_file = inputs.mapping_file.import_file(imported_format="compressed")
            stem = Path(mapping_file).stem
            supported_extensions = (".tab", ".tsv", ".txt")
            if not stem.endswith(supported_extensions):
                self.error(
                    "Mapping file has unsupported file name extension. "
                    f"The supported extensions are {supported_extensions}."
                )
            mapping = pd.read_csv(
                mapping_file,
                sep="\t",
                header=0,
                names=["ensembl_id", "probe"],
                dtype=str,
            )
            mapping = mapping.drop_duplicates()

            if inputs.source:
                source = inputs.source
            else:
                self.error(
                    "Custom probe id mapping file was provided but no source was selected."
                )
            if inputs.build:
                build = inputs.build
            else:
                self.error(
                    "Custom probe id mapping file was provided but genome build was not defined."
                )
            probe_mapping = "Custom"
        else:
            if not platform_id:
                self.error(
                    "Custom mapping file should be provided when samples do not have a GEO platform defined"
                )
            if platform_id not in PLATFORM_MAP:
                self.error(f"GEO platform {platform_id} is not supported.")

            species_low = species.lower()
            dataset = f"{species_low[0]}{species_low.split(' ')[1]}_gene_ensembl"
            probe_mapping = PLATFORM_MAP[platform_id]

            try:
                b = BioMart()
            except IOError:
                raise Exception("None of the ENSEMBL Biomart hosts is reachable.")
            except Exception as e:
                raise Exception(f"Unexpected biomart error: {e}")

            b.add_dataset_to_xml(dataset)
            b.add_attribute_to_xml("ensembl_gene_id")
            b.add_attribute_to_xml(probe_mapping)  # type of microarray
            b.add_filter_to_xml(probe_mapping, ",".join(probe_ids))
            xml_query = b.get_xml()
            res = b.query(xml_query)

            mapping = pd.read_csv(
                StringIO(res),
                sep="\t",
                header=None,
                names=["ensembl_id", "probe"],
                dtype=str,
            )
            mapping = mapping.drop_duplicates()
            mapping_file = f"{platform}_mapping.tsv"
            mapping.to_csv(mapping_file, sep="\t", index=False)

            dataset_names = b.get_datasets("ENSEMBL_MART_ENSEMBL")
            display_name = dataset_names.loc[dataset_names["name"] == dataset][
                "description"
            ].to_string()
            # Typical display name would be Human genes (GRCh38.p13)
            build = re.search(r"\((.+?)\)", display_name).group(1)
            source = "ENSEMBL"

        mapping = mapping.drop_duplicates(subset=["probe"], keep=False)

        data = joined_expressions.loc[mapping["probe"]]
        data["ensembl_id"] = mapping["ensembl_id"].tolist()
        data = data.reset_index()

        # For Ensembl IDs with multiple probe IDs retain the one with highest expression.
        data["mean"] = data.loc[
            :, data.columns.difference(["probe", "ensembl_id"])
        ].mean(axis=1)
        idx_max = data.groupby(["ensembl_id"])["mean"].idxmax()
        data = data.loc[idx_max].set_index("ensembl_id")

        data = data.drop(columns=["probe", "mean"])
        data.index.name = "Gene"

        mapped_file = "mapped_expressions.tsv.gz"
        data.to_csv(mapped_file, sep="\t", index=True, compression="gzip")
        for column, exp in zip(data.columns, inputs.expressions):
            mapped_column = f"{column}_mapped_exp.tsv.gz"
            data.to_csv(
                mapped_column,
                sep="\t",
                index=True,
                columns=[column],
                header=["Expression"],
                index_label="Gene",
                compression="gzip",
            )
            self.run_process(
                "mapped-microarray-expression",
                {
                    "exp_unmapped": exp.id,
                    "exp": mapped_column,
                    "source": source,
                    "build": build,
                    "probe_mapping": probe_mapping,
                },
            )

        outputs.mapped_exp = mapped_file
        outputs.mapping = mapping_file
        outputs.probe_mapping = probe_mapping
        outputs.platform = platform
        if platform_id:
            outputs.platform_id = platform_id
