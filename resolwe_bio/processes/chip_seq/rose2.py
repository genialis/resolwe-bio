"""Run ROSE2."""

from pathlib import Path
from shutil import copy
from zipfile import ZipFile

from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    FileField,
    IntegerField,
    JsonField,
    Process,
    StringField,
)


class Rose2(Process):
    """Run ROSE2.

    Rank Ordering of Super-Enhancers algorithm (ROSE2) takes the acetylation
    peaks called by a peak caller (MACS, MACS2...) and based on the in-between
    distances and the acetylation signal at the peaks judges whether they can
    be considered super-enhancers. The ranked values are plotted and by
    locating the inflection point in the resulting graph, super-enhancers are
    assigned. See [here](http://younglab.wi.mit.edu/super_enhancer_code.html)
    for more information.
    """

    slug = "rose2"
    name = "ROSE2"
    process_type = "data:chipseq:rose2"
    version = "5.0.0"
    category = "ChIP-Seq:Post Process"
    entity = {
        "type": "sample",
        "input": "{{ input_macs if input_macs else input_upload }}",
    }
    requirements = {
        "expression-engine": "jinja",
        "executor": {"docker": {"image": "resolwebio/bamliquidator:1.2.0"}},
    }
    data_name = "{{ input_macs|sample_name|default('?') if input_macs else input_upload|sample_name|default('?') }}"

    class Input:
        """Input fields to process ROSE2."""

        input_macs = DataField(
            "chipseq:callpeak",
            label="BED/narrowPeak file (MACS results)",
            required=False,
            hidden="input_upload",
        )
        input_upload = DataField(
            "bed",
            label="BED file (Upload)",
            required=False,
            hidden="input_macs || use_filtered_bam",
        )
        use_filtered_bam = BooleanField(
            label="Use Filtered BAM File",
            default=False,
            hidden="input_upload",
            description=(
                "Use filtered BAM file from a MACS2 object to rank "
                "enhancers by. Only applicable if input is MACS2."
            ),
        )
        rankby = DataField(
            "alignment:bam",
            label="BAM file",
            required=False,
            hidden="use_filtered_bam",
            description="BAM file to rank enhancers by.",
        )
        control = DataField(
            "alignment:bam",
            label="Control BAM File",
            required=False,
            hidden="use_filtered_bam",
            description="BAM file to rank enhancers by.",
        )
        tss = IntegerField(
            label="TSS exclusion",
            default=0,
            description="Enter a distance from TSS to exclude. 0 = no TSS exclusion.",
        )
        stitch = IntegerField(
            label="Stitch",
            required=False,
            description=(
                "Enter a max linking distance for stitching. If not "
                "given, optimal stitching parameter will be determined"
                " automatically."
            ),
        )
        mask = DataField(
            "bed",
            label="Masking BED file",
            required=False,
            description=(
                "Mask a set of regions from analysis. Provide a BED of"
                " masking regions."
            ),
        )

    class Output:
        """Output field of the process ImportFastaNucleotide."""

        all_enhancers = FileField(label="All enhancers table")
        enhancers_with_super = FileField(label="Super enhancers table")
        plot_points = FileField(label="Plot points")
        plot_panel = FileField(label="Plot panel")
        enhancer_gene = FileField(label="Enhancer to gene")
        enhancer_top_gene = FileField(label="Enhancer to top gene")
        gene_enhancer = FileField(label="Gene to Enhancer")
        stitch_parameter = FileField(label="Stitch parameter", required=False)
        all_output = FileField(label="All output")
        scatter_plot = JsonField(label="Super-Enhancer plot")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run the analysis."""
        if not inputs.input_macs and not inputs.input_upload:
            self.error(
                "Peaks file missing. Please provide .bed peaks file as a file "
                "upload or a MACS output."
            )

        if inputs.input_macs and inputs.input_upload:
            self.error("Please provide only one .bed peaks file.")

        if inputs.control and not inputs.rankby:
            self.error(
                "A control BAM file cannot be provided without specifying "
                "a BAM file to rank by. If selecting Use Filtered Bam File"
                " option neither must be specified."
            )

        if inputs.input_macs:
            if (
                inputs.input_macs.type == "data:chipseq:callpeak:macs14:"
                and inputs.use_filtered_bam
            ):
                self.error(
                    "Use Filtered Bam File option can only be used with a "
                    "MACS2 input."
                )
            if (
                not (
                    inputs.input_macs.type == "data:chipseq:callpeak:macs2:"
                    and inputs.use_filtered_bam
                )
                and not inputs.rankby
            ):
                self.error(
                    "BAM file to rank by must be used unless a filtered BAM "
                    "from a MACS2 data oject is used to rank by."
                )

            species = inputs.input_macs.species
            build = inputs.input_macs.build
            if inputs.rankby:
                if species != inputs.rankby.species:
                    self.error(
                        f"Species of rankby bam file {inputs.rankby.species} "
                        f"and MACS bed file {species} do not match. Please "
                        "provide aligned reads and annotation with the same "
                        "species."
                    )
                if build != inputs.rankby.build:
                    self.error(
                        "Genome builds of rankby bam file "
                        f"{inputs.rankby.build} and MACS bed file {build} "
                        "do not match. Please provide aligned reads and "
                        "annotation with the same genome build."
                    )

        if inputs.input_upload:
            if not inputs.rankby:
                self.error(
                    "BAM file to rank by must be used unless a filtered BAM "
                    "from a MACS2 data oject is used to rank by."
                )

            species = inputs.input_upload.species
            if species != inputs.rankby.species:
                self.error(
                    f"Species of rankby bam file {inputs.rankby.species} and "
                    f"uploaded bed file {species} do not match. Please provide"
                    " aligned reads and annotation with the same species."
                )

            build = inputs.input_upload.build
            if build != inputs.rankby.build:
                self.error(
                    f"Genome builds of rankby bam file {inputs.rankby.build} "
                    f"and uploaded bed file {build} do not match. Please "
                    "provide aligned reads and annotation with the same genome"
                    " build."
                )

        if inputs.control:
            if inputs.control.species != inputs.rankby.species:
                self.error(
                    f"Species of rankby bam file {inputs.rankby.species} and "
                    f"control bam file {inputs.control.species} do not match. "
                    "Please provide aligned reads with the same species."
                )
            if inputs.control.build != inputs.rankby.build:
                self.error(
                    f"Genome builds of rankby bam file {inputs.rankby.build} "
                    f"and control bam file {inputs.control.build} do not "
                    "match. Please provide aligned reads with the same genome "
                    "build."
                )

        if inputs.mask:
            if inputs.mask.species != species:
                self.error(
                    f"Species of the masking bed file {inputs.mask.species} "
                    "does not match other inputs` species. Please provide a "
                    "masking file of the same species as other inputs."
                )
            if inputs.mask.build != build:
                self.error(
                    f"Genome build of the masking bed file {inputs.mask.build}"
                    " does not match other inputs` build. Please provide a "
                    "masking file of the same genome build as other inputs."
                )

        genome_list = [
            "hg19",
            "hg18",
            "mm10",
            "mm9",
            "mm8",
            "rn6",
            "rn4",
        ]
        if build not in genome_list:
            self.error(
                f"{build} is not a valid genome build. Accepted genome "
                f'builds are: {", ".join(genome_list)}.'
            )

        if inputs.input_upload:
            bed_path = Path(inputs.input_upload.bed.path)
            if bed_path.name[-4:] == ".bed":
                name = bed_path.stem
                copy(str(bed_path), bed_path.name)
            else:
                name = bed_path.name[:-7]
                cmd = Cmd["bgzip"]["-cd"][str(bed_path)]
                (cmd > f"{name}.bed")()

        elif inputs.input_macs.type == "data:chipseq:callpeak:macs14:":
            bed_path = Path(inputs.input_macs.peaks_bed.path)
            if bed_path.name[-4:] == ".bed":
                name = bed_path.stem
                copy(str(bed_path), bed_path.name)
            else:
                name = bed_path.name[:-7]
                cmd = Cmd["bgzip"]["-cd"][str(bed_path)]
                (cmd > f"{name}.bed")()

        elif inputs.input_macs.type == "data:chipseq:callpeak:macs2:":
            narrowpeak_path = Path(inputs.input_macs.narrow_peaks.path)
            if narrowpeak_path.name[-11:] == ".narrowPeak":
                name = narrowpeak_path.stem
                copy(str(narrowpeak_path), f"{name}.bed")
            else:
                name = narrowpeak_path.name[:-14]
                cmd = Cmd["bgzip"]["-cd"][str(narrowpeak_path)]
                (cmd > f"{name}.bed")()

        if inputs.input_macs and inputs.use_filtered_bam:
            rankby = inputs.input_macs.case_bam.path
        else:
            rankby = inputs.rankby.bam.path

        if (
            inputs.input_macs
            and inputs.input_macs.type == "data:chipseq:callpeak:macs2:"
            and inputs.input_macs.control_bam
            and inputs.use_filtered_bam
        ):
            control = inputs.input_macs.control_bam.path
        elif inputs.control:
            control = inputs.control.bam.path
        else:
            control = None

        cmd = Cmd["rose2"]
        cmd = cmd["--genome"][build.upper()]
        cmd = cmd["-i"][f"{name}.bed"]
        cmd = cmd["--rankby"][rankby]
        if control:
            cmd = cmd["--control"][control]
        cmd = cmd["--tss"][inputs.tss]
        if inputs.stitch or inputs.stitch == 0:
            cmd = cmd["--stitch"][inputs.stitch]
        if inputs.mask:
            cmd = cmd["--mask"][inputs.mask.bed.path]
        cmd = cmd["--out"]["."]
        return_code, _, _ = cmd & TEE(retcode=None)
        if return_code:
            self.error(f"ROSE2 run failed.")

        outputs.all_enhancers = f"{name}_AllEnhancers.table.txt"
        outputs.enhancers_with_super = f"{name}_Enhancers_withSuper.bed"
        outputs.plot_points = f"{name}_Plot_points.png"
        outputs.plot_panel = f"{name}_Plot_panel.png"
        outputs.enhancer_gene = f"{name}_SuperEnhancers_ENHANCER_TO_GENE.txt"
        outputs.enhancer_top_gene = f"{name}_SuperEnhancers_ENHANCER_TO_TOP_GENE.txt"
        outputs.gene_enhancer = f"{name}_SuperEnhancers_GENE_TO_ENHANCER.txt"
        if not (inputs.stitch or inputs.stitch == 0):
            outputs.stitch_parameter = f"{name}_stitch_parameter.pdf"

        cmd = Cmd["plot_enhancers.py"]
        cmd = cmd[f"{name}_AllEnhancers.table.txt"]
        cmd = cmd["scatter_plot.json"]
        if inputs.control:
            cmd["-c"]
        return_code, _, _ = cmd & TEE(retcode=None)
        if return_code:
            self.error(f"Plotting enhancers failed.")

        zipfile = f"{name}_output_all.zip"
        with ZipFile(zipfile, "w") as zip_file:
            for f in Path(".").glob(f"{name}_*"):
                if f.name == zipfile:
                    continue
                zip_file.write(f)

        outputs.all_output = zipfile
        outputs.scatter_plot = "scatter_plot.json"

        outputs.species = species
        outputs.build = build
