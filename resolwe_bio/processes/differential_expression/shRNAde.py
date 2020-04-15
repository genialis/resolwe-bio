"""Differential expression of shRNA."""
import gzip
import os
from shutil import copy

from resolwe.process import (
    Cmd,
    DataField,
    FileField,
    ListField,
    Process,
    SchedulingClass,
)


class ShortHairpinRNADifferentialExpression(Process):
    """
    Performing differential expression on a list of objects.

    Analysis starts by inputting a set of expression files (count matrices) and a parameter file. Parameter file is
    an xlsx file and consists of tabs:
        - `sample_key`: Should have column sample with exact sample name as input expression file(s),
        columns defining treatment and lastly a column which indicates replicate.
        - `contrasts`: Define groups which will be used to perform differential expression analysis. Model for DE
        uses these contrasts and replicate number. In R annotation, this would be ` ~ 1 + group + replicate`. Table
        should have two columns named `group_1` and `group_2`.
        - `overall_contrasts`: This is a layer "above" `contrasts`, where results from two contrasts are compared for
        lethal, beneficial and neutral species. Thresholds governing classification can be found in
        `classification_parameters` tab.
        - `classification_parameters`: This tab holds three columns, `threshold`, `value` and `description`. Only
        the first two are used in the workflow, description is for your benefit.

    This process outputs DESeq2 results, classified results based on provided thresholds and counts of beneficial
    and lethal species.
    """

    slug = "differentialexpression-shrna"
    name = "Differential expression of shRNA"
    process_type = "data:shrna:differentialexpression:"
    version = "1.1.1"
    category = "Differential Expression"
    scheduling_class = SchedulingClass.BATCH
    entity = {"type": "sample"}
    requirements = {
        "expression-engine": "jinja",
        "executor": {"docker": {"image": "resolwebio/rnaseq:4.9.0"}},
    }
    data_name = '{{ parameter_file.file|default("?") }}'

    class Input:
        """Input fields to process ShortHairpinRNADifferentialExpression."""

        parameter_file = DataField(
            data_type="file",
            label="Excel parameter file (.xlsx)",
            description="Select .xlsx file which holds parameters for analysis. "
            "See [here](https://github.com/genialis/shRNAde/blob/master/inst/extdata/template_doDE_inputs"
            ".xlsx) for a template.",
        )
        expression_data = ListField(
            DataField(
                data_type="expression:shrna2quant:",
                description="Data objects of expressions from process shrna-quant. These inputs should match sample "
                "names specified in parameter file.",
            ),
            label="List of expression files from shrna2quant",
        )

    class Output:
        """Output fields to process ShortHairpinRNADifferentialExpression."""

        deseq_results = FileField(label="DESeq2 results")
        class_results = FileField(
            label="Results classified based on thresholds provided by the user"
        )
        beneficial_counts = FileField(
            label="shRNAs considered as beneficial based on user input"
        )
        lethal_counts = FileField(
            label="shRNAs considered as lethal based on user input"
        )

    def run(self, inputs, outputs):
        """Run differential expression of shRNA.

        These are the steps for running the process:

          1.) Prepare data to be pulled into R for processing. Place data
              objects into expression_files/ folder.
          2.) Pass parameter file for R's shRNAde::doDE() and execute the
              function.
          3.) Prepare outputs.
        """
        dir_expressions = "./expression_files"
        os.mkdir(dir_expressions)

        # (1) Move expression files and extract files.
        sample_list = [
            copy(src=x.exp.path, dst=dir_expressions) for x in inputs.expression_data
        ]

        for fl in sample_list:
            base_filename = os.path.splitext(fl)[0]
            with gzip.open(fl) as in_file:
                with open(base_filename, "wt") as out_file:
                    for line in in_file:
                        out_file.write(line.decode())

        # (2)
        r_input = f'shRNAde::doDE(input = "{inputs.parameter_file.file.path}", sample_list = "{dir_expressions}")'

        run_cmd = Cmd["Rscript"]["-e"][r_input]
        run_cmd()

        # (3) Compress results before storing them.
        result_deseq = "deseq_results.txt"
        result_class = "class_results.txt"
        result_beneficial = "beneficial_counts.txt"
        result_lethal = "lethal_counts.txt"

        to_compress = [result_deseq, result_class, result_beneficial, result_lethal]

        for file in to_compress:
            with open(file, "rb") as txt, gzip.open(file + ".gz", "wb") as gz:
                try:
                    gz.writelines(txt)
                except FileNotFoundError:
                    return "Something went wrong during result compression."

        outputs.deseq_results = result_deseq + ".gz"
        outputs.class_results = result_class + ".gz"
        outputs.beneficial_counts = result_beneficial + ".gz"
        outputs.lethal_counts = result_lethal + ".gz"
