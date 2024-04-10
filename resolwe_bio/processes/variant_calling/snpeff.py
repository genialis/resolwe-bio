"""Variant annotation with SnpEff."""

import gzip
import shutil
from pathlib import Path

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    FileField,
    FileHtmlField,
    GroupField,
    ListField,
    Persistence,
    Process,
    SchedulingClass,
    StringField,
)


def return_sample_count(vcf, error):
    """Count number of samples in the input VCF file."""
    try:
        with gzip.open(vcf, "rt") as vcf_in:
            for line in vcf_in:
                if line.startswith("#CHROM"):
                    headers = line.rstrip().split()
                    return len(headers[headers.index("FORMAT") + 1 :])

    except Exception as err:
        error(
            f"Unable to determine sample count in VCF file. Original error was: {err}"
        )


class SnpEff(Process):
    """Annotate variants with SnpEff.

    SnpEff is a variant annotation and effect prediction tool.
    It annotates and predicts the effects of genetic variants
    (such as amino acid changes).

    This process also allows filtering of variants with ``SnpSift
    filter`` command and extracting specific fields from the VCF
    file with ``SnpSift extractFields`` command.

    This tool works with multi-sample VCF file as an input.
    """

    slug = "snpeff"
    name = "snpEff (General variant annotation) (multi-sample)"
    process_type = "data:variants:vcf:snpeff"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/snpeff:3.0.0"},
        },
        "resources": {
            "cores": 2,
            "memory": 16384,
        },
    }
    category = "WGS"
    data_name = "Annotated variants (SnpEff)"
    version = "1.2.0"
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.CACHED

    class Input:
        """Input fields to SnpEff process."""

        variants = DataField(
            data_type="variants:vcf",
            label="Variants (VCF)",
        )
        database = StringField(
            label="snpEff database",
            default="GRCh38.109",
            choices=[
                ("GRCh37.75", "GRCh37.75"),
                ("GRCh38.99", "GRCh38.99"),
                ("GRCh38.109", "GRCh38.109"),
            ],
        )
        dbsnp = DataField(
            data_type="variants:vcf",
            label="Known variants",
            description="List of known variants for annotation.",
            required=False,
        )
        filtering_options = StringField(
            label="Filtering expressions",
            description="Filter VCF file using arbitraty expressions."
            "Examples of filtering expressions: '(ANN[*].GENE = 'PSD3')' "
            "or '( REF = 'A' )' or "
            "'(countHom() > 3) | (( exists INDEL ) & (QUAL >= 20)) | (QUAL >= 30 )'."
            "For more information checkout the official documentation of [SnpSift]"
            "(https://pcingola.github.io/SnpEff/ss_filter/)",
            required=False,
        )
        sets = ListField(
            DataField(data_type="geneset"),
            label="Files with list of genes",
            description="Use list of genes, if you only want variants reported for "
            "them. Each file must have one string per line.",
            hidden="!filtering_options",
            required=False,
        )
        extract_fields = ListField(
            StringField(),
            label="Fields to extract",
            description="Write fields you want to extract from annonated vcf file "
            "and press Enter after each one. Example of fields: `CHROM POS REF ALT "
            "'ANN[*].GENE'`. For more information follow this [link]"
            "(https://pcingola.github.io/SnpEff/ss_extractfields/).",
            required=False,
        )

        class Advanced:
            """Advanced options."""

            one_per_line = BooleanField(
                label="One effect per line",
                default=False,
                description="If there is more than one effect per variant, write them "
                "to seperate lines.",
            )

        advanced = GroupField(
            Advanced, label="Advanced options", hidden="!extract_fields"
        )

    class Output:
        """Output fields to process SnpEff."""

        vcf = FileField(
            label="Annotated variants (VCF)",
        )
        tbi = FileField(label="Index of annotated variants")
        vcf_extracted = FileField(
            label="Extracted annotated variants (VCF)",
            required=False,
        )
        tbi_extracted = FileField(
            label="Index of extracted variants",
            required=False,
        )
        species = StringField(label="Species")
        build = StringField(label="Build")
        genes = FileField(label="SnpEff genes")
        summary = FileHtmlField(
            label="Summary",
        )

    def run(self, inputs, outputs):
        """Run analysis."""

        output_variants = "snpeff_variants.vcf"
        annotated_variants = "annotated_variants.vcf"
        filtered_variants = "filtered_variants.vcf"
        extracted_variants = "extracted_variants.vcf"

        if not inputs.variants.output.build.startswith(inputs.database[:6]):
            self.error(
                "Genome build for the input variants file and "
                "SnpEff database should be the same. Input variants file is "
                f"based on {inputs.variants.output.build}, while SnpEff "
                f"database is based on {inputs.database[:6]}."
            )

        file_name = Path(inputs.variants.output.vcf.path).name
        shutil.copy(Path(inputs.variants.output.vcf.path), Path.cwd())

        # check the VCF file content
        sample_count = return_sample_count(vcf=Path(file_name), error=self.error)
        if not sample_count > 1:
            self.error(
                f"The input VCF file should contain data for multiple samples. "
                f"The input contains data for {sample_count} sample(s)."
            )

        args_snpeff = [
            inputs.database,
            inputs.variants.output.vcf.path,
        ]
        (Cmd["snpEff"][args_snpeff] > output_variants)()

        if inputs.dbsnp:
            if not inputs.dbsnp.output.build.startswith(inputs.database[:6]):
                self.error(
                    "Genome build for the DBSNP file and used database "
                    "should be the same. DBSNP file is based on "
                    f"{inputs.dbsnp.output.build}, while snpEff database "
                    f"is based on {inputs.database[:6]}."
                )
            args_annotation = [
                "annotate",
                inputs.dbsnp.output.vcf.path,
                output_variants,
            ]

            (Cmd["SnpSift"][args_annotation] > annotated_variants)()

        if inputs.filtering_options:
            args_filtering = [
                "filter",
                inputs.filtering_options,
            ]

            if inputs.dbsnp:
                args_filtering.append(annotated_variants)
            else:
                args_filtering.append(output_variants)

            if inputs.sets:
                for set in inputs.sets:
                    args_filtering.extend(["-s", set.output.geneset.path])

            (Cmd["SnpSift"][args_filtering] > filtered_variants)()
            (Cmd["bgzip"]["-c", filtered_variants] > filtered_variants + ".gz")()
            (Cmd["tabix"]["-p", "vcf", filtered_variants + ".gz"])()

            outputs.vcf = filtered_variants + ".gz"
            outputs.tbi = filtered_variants + ".gz.tbi"

        elif inputs.dbsnp:
            (Cmd["bgzip"]["-c", annotated_variants] > annotated_variants + ".gz")()
            (Cmd["tabix"]["-p", "vcf", annotated_variants + ".gz"])()

            outputs.vcf = annotated_variants + ".gz"
            outputs.tbi = annotated_variants + ".gz.tbi"
        else:
            (Cmd["bgzip"]["-c", output_variants] > output_variants + ".gz")()
            (Cmd["tabix"]["-p", "vcf", output_variants + ".gz"])()

            outputs.vcf = output_variants + ".gz"
            outputs.tbi = output_variants + ".gz.tbi"

        if inputs.extract_fields:
            args_extract = [
                "extractFields",
            ]

            if not inputs.dbsnp and not inputs.filtering_options:
                args_extract.append(output_variants)
            elif inputs.dbsnp and not inputs.filtering_options:
                args_extract.append(annotated_variants)
            else:
                args_extract.append(filtered_variants)

            for field in inputs.extract_fields:
                args_extract.append(field)

            if inputs.advanced.one_per_line:
                (
                    Cmd["SnpSift"][args_extract] | Cmd["vcfEffOnePerLine.pl"]
                    > extracted_variants
                )()
            else:
                (Cmd["SnpSift"][args_extract] > extracted_variants)()
            (Cmd["bgzip"]["-c", extracted_variants] > extracted_variants + ".gz")()
            (Cmd["tabix"]["-p", "vcf", extracted_variants + ".gz"])()

            outputs.vcf_extracted = extracted_variants + ".gz"
            outputs.tbi_extracted = extracted_variants + ".gz.tbi"

        outputs.species = inputs.variants.output.species
        outputs.build = inputs.variants.output.build
        outputs.genes = "snpEff_genes.txt"
        outputs.summary = "snpEff_summary.html"


class SnpEffSingleSample(Process):
    """Annotate variants with SnpEff.

    SnpEff is a variant annotation and effect prediction tool.
    It annotates and predicts the effects of genetic variants
    (such as amino acid changes).

    This process also allows filtering of variants with ``SnpSift
    filter`` command and extracting specific fields from the VCF
    file with ``SnpSift extractFields`` command.

    This tool works with single-sample VCF file as an input.
    """

    slug = "snpeff-single"
    name = "snpEff (General variant annotation) (single-sample)"
    process_type = "data:variants:vcf:snpeff:single"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/snpeff:3.0.0"},
        },
        "resources": {
            "cores": 2,
            "memory": 16384,
        },
    }
    entity = {
        "type": "sample",
    }
    category = "WGS"
    data_name = "{{ variants|name|default('?') }}"
    version = "1.1.0"
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.CACHED

    class Input:
        """Input fields to SnpEffSingleSample process."""

        variants = DataField(
            data_type="variants:vcf",
            label="Variants (VCF)",
        )
        database = StringField(
            label="snpEff database",
            default="GRCh38.109",
            choices=[
                ("GRCh37.75", "GRCh37.75"),
                ("GRCh38.99", "GRCh38.99"),
                ("GRCh38.109", "GRCh38.109"),
            ],
        )
        dbsnp = DataField(
            data_type="variants:vcf",
            label="Known variants",
            description="List of known variants for annotation.",
            required=False,
        )
        filtering_options = StringField(
            label="Filtering expressions",
            description="Filter VCF file using arbitraty expressions."
            "Examples of filtering expressions: '(ANN[*].GENE = 'PSD3')' "
            "or '( REF = 'A' )' or "
            "'(countHom() > 3) | (( exists INDEL ) & (QUAL >= 20)) | (QUAL >= 30 )'."
            "For more information checkout the official documentation of [SnpSift]"
            "(https://pcingola.github.io/SnpEff/ss_filter/)",
            required=False,
        )
        sets = ListField(
            DataField(data_type="geneset"),
            label="Files with list of genes",
            description="Use list of genes, if you only want variants reported for "
            "them. Each file must have one string per line.",
            hidden="!filtering_options",
            required=False,
        )
        extract_fields = ListField(
            StringField(),
            label="Fields to extract",
            description="Write fields you want to extract from annonated vcf file "
            "and press Enter after each one. Example of fields: `CHROM POS REF ALT "
            "'ANN[*].GENE'`. For more information follow this [link]"
            "(https://pcingola.github.io/SnpEff/ss_extractfields/).",
            required=False,
        )

        class Advanced:
            """Advanced options."""

            one_per_line = BooleanField(
                label="One effect per line",
                default=False,
                description="If there is more than one effect per variant, write them "
                "to seperate lines.",
            )

        advanced = GroupField(
            Advanced, label="Advanced options", hidden="!extract_fields"
        )

    class Output:
        """Output fields to process SnpEffSingleSample."""

        vcf = FileField(
            label="Annotated variants (VCF)",
        )
        tbi = FileField(label="Index of annotated variants")
        vcf_extracted = FileField(
            label="Extracted annotated variants (VCF)",
            required=False,
        )
        tbi_extracted = FileField(
            label="Index of extracted variants",
            required=False,
        )
        species = StringField(label="Species")
        build = StringField(label="Build")
        genes = FileField(label="SnpEff genes")
        summary = FileHtmlField(
            label="Summary",
        )

    def run(self, inputs, outputs):
        """Run analysis."""

        output_variants = "snpeff_variants.vcf"
        annotated_variants = "annotated_variants.vcf"
        filtered_variants = "filtered_variants.vcf"
        extracted_variants = "extracted_variants.vcf"

        if not inputs.variants.output.build.startswith(inputs.database[:6]):
            self.error(
                "Genome build for the input variants file and "
                "SnpEff database should be the same. Input variants file is "
                f"based on {inputs.variants.output.build}, while SnpEff "
                f"database is based on {inputs.database[:6]}."
            )

        file_name = Path(inputs.variants.output.vcf.path).name
        shutil.copy(Path(inputs.variants.output.vcf.path), Path.cwd())

        # check the VCF file content
        sample_count = return_sample_count(vcf=Path(file_name), error=self.error)
        if sample_count != 1:
            self.error(
                f"The input VCF should contain data for a single sample. "
                f"The input contains data for {sample_count} sample(s)."
            )

        args_snpeff = [
            inputs.database,
            inputs.variants.output.vcf.path,
        ]
        (Cmd["snpEff"][args_snpeff] > output_variants)()

        if inputs.dbsnp:
            if not inputs.dbsnp.output.build.startswith(inputs.database[:6]):
                self.error(
                    "Genome build for the DBSNP file and used database "
                    "should be the same. DBSNP file is based on "
                    f"{inputs.dbsnp.output.build}, while snpEff database "
                    f"is based on {inputs.database[:6]}."
                )
            args_annotation = [
                "annotate",
                inputs.dbsnp.output.vcf.path,
                output_variants,
            ]

            (Cmd["SnpSift"][args_annotation] > annotated_variants)()

        if inputs.filtering_options:
            args_filtering = [
                "filter",
                inputs.filtering_options,
            ]

            if inputs.dbsnp:
                args_filtering.append(annotated_variants)
            else:
                args_filtering.append(output_variants)

            if inputs.sets:
                for set in inputs.sets:
                    args_filtering.extend(["-s", set.output.geneset.path])

            (Cmd["SnpSift"][args_filtering] > filtered_variants)()
            (Cmd["bgzip"]["-c", filtered_variants] > filtered_variants + ".gz")()
            (Cmd["tabix"]["-p", "vcf", filtered_variants + ".gz"])()

            outputs.vcf = filtered_variants + ".gz"
            outputs.tbi = filtered_variants + ".gz.tbi"

        elif inputs.dbsnp:
            (Cmd["bgzip"]["-c", annotated_variants] > annotated_variants + ".gz")()
            (Cmd["tabix"]["-p", "vcf", annotated_variants + ".gz"])()

            outputs.vcf = annotated_variants + ".gz"
            outputs.tbi = annotated_variants + ".gz.tbi"
        else:
            (Cmd["bgzip"]["-c", output_variants] > output_variants + ".gz")()
            (Cmd["tabix"]["-p", "vcf", output_variants + ".gz"])()

            outputs.vcf = output_variants + ".gz"
            outputs.tbi = output_variants + ".gz.tbi"

        if inputs.extract_fields:
            args_extract = [
                "extractFields",
            ]

            if not inputs.dbsnp and not inputs.filtering_options:
                args_extract.append(output_variants)
            elif inputs.dbsnp and not inputs.filtering_options:
                args_extract.append(annotated_variants)
            else:
                args_extract.append(filtered_variants)

            for field in inputs.extract_fields:
                args_extract.append(field)

            if inputs.advanced.one_per_line:
                (
                    Cmd["SnpSift"][args_extract] | Cmd["vcfEffOnePerLine.pl"]
                    > extracted_variants
                )()
            else:
                (Cmd["SnpSift"][args_extract] > extracted_variants)()
            (Cmd["bgzip"]["-c", extracted_variants] > extracted_variants + ".gz")()
            (Cmd["tabix"]["-p", "vcf", extracted_variants + ".gz"])()

            outputs.vcf_extracted = extracted_variants + ".gz"
            outputs.tbi_extracted = extracted_variants + ".gz.tbi"

        outputs.species = inputs.variants.output.species
        outputs.build = inputs.variants.output.build
        outputs.genes = "snpEff_genes.txt"
        outputs.summary = "snpEff_summary.html"
