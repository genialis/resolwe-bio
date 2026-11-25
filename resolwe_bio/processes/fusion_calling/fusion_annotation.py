"""Clinical annotation of gene fusions."""

import csv
import gzip
import json
from collections import Counter
from pathlib import Path

from resolwe.process import (
    DataField,
    FileField,
    Persistence,
    Process,
    SchedulingClass,
    StringField,
)

ARRIBA_COLUMNS = (
    "#gene1",
    "gene2",
    "type",
    "reading_frame",
    "retained_protein_domains",
)

SNAPSHOT_SECTIONS = ("versions", "civic", "chimerkb", "normal_tissue")

ANNOTATION_COLUMNS = (
    "civic_match",
    "civic_best_levels",
    "civic_evidence",
    "civic_pubmed_ids",
    "civic_nct_ids",
    "civic_any_partner_profiles",
    "civic_any_partner_best_levels",
    "civic_any_partner_evidence",
    "civic_any_partner_pubmed_ids",
    "civic_any_partner_nct_ids",
    "chimerkb_diseases",
    "chimerkb_validation",
    "retained_kinase_domain",
    "in_frame_kinase_fusion",
    "normal_tissue_sources",
    "flags",
)

# Separator of the therapies of one CIViC evidence item by interaction type.
# Alternative therapies ("Substitutes") must not read as a combination.
THERAPY_SEPARATORS = {"Combination": "+", "Sequential": ">", "Substitutes": "/"}

# Arriba marks genes fused head-to-head or tail-to-tail with these suffixes
# of the event type. Such an event cannot yield a chimeric protein.
NO_CHIMERIC_PROTEIN_TYPES = ("5'-5'", "3'-3'")


def clean_gene(name):
    """Return the gene symbol of an Arriba gene field.

    Arriba reports an intergenic breakpoint as nearby genes with distances,
    e.g. ``MS4A2(45793),LINC02705(15427)``, and a breakpoint with no annotated
    gene as ``.``. Such a field has no single gene, so the function returns
    ``None``.
    """
    if name == "." or "(" in name or "," in name:
        return None
    return name.strip().upper()


def best_levels(evidence):
    """Return the strongest CIViC evidence level of each evidence type.

    Levels range from ``A`` (strongest) to ``E``. The levels are reported for
    each evidence type, e.g. ``Diagnostic:A,Predictive:C``, because a strong
    diagnostic claim says nothing about the response to a therapy. Evidence
    with direction ``Does Not Support`` argues against the claim, so it does
    not count.
    """
    levels = {}
    for item in evidence:
        if item["direction"] == "Supports":
            evidence_type = item["type"]
            levels[evidence_type] = min(
                item["level"], levels.get(evidence_type, item["level"])
            )
    return ",".join(
        f"{evidence_type}:{levels[evidence_type]}" for evidence_type in sorted(levels)
    )


def format_therapies(item):
    """Return the therapies of a CIViC evidence item.

    Therapies given in combination are joined with ``+``, sequential therapies
    with ``>`` and alternative therapies with ``/``.
    """
    separator = THERAPY_SEPARATORS.get(item["therapy_interaction_type"], "/")
    return separator.join(item["therapies"].split(","))


def summarize_evidence(evidence):
    """Summarize supporting CIViC evidence, strongest level first.

    Each distinct claim gives one entry of the form
    ``level: disease, type, significance, therapies``.
    """
    claims = {}
    for item in evidence:
        if item["direction"] != "Supports":
            continue
        claim = (
            item["disease"],
            item["type"],
            item["significance"],
            format_therapies(item),
        )
        claims[claim] = min(item["level"], claims.get(claim, item["level"]))
    entries = sorted((level, *claim) for claim, level in claims.items())
    return " | ".join(
        f"{level}: {', '.join(part for part in claim if part)}"
        for level, *claim in entries
    )


def pubmed_ids(evidence):
    """Return the sorted unique PubMed IDs of supporting CIViC evidence."""
    ids = {
        item["citation_id"]
        for item in evidence
        if item["direction"] == "Supports" and item["source_type"] == "PubMed"
    }
    return ",".join(sorted(ids, key=int))


def nct_ids(evidence):
    """Return the sorted unique clinical trial IDs of supporting CIViC evidence."""
    ids = {
        nct_id.strip()
        for item in evidence
        if item["direction"] == "Supports"
        for nct_id in item["nct_ids"].split(",")
        if nct_id.strip()
    }
    return ",".join(sorted(ids))


def retained_kinase_domain(retained_protein_domains):
    """Return the protein kinase domains that the fusion retains.

    Arriba lists the retained domains of the 5' and 3' genes separated by a
    pipe, each with the retained fraction, e.g.
    ``SH3_domain(100%)|Protein_kinase_domain(100%),SH2_domain(100%)``. The
    result names the partner that retains the domain, e.g.
    ``3':Protein_kinase_domain(100%)``.
    """
    five_prime, three_prime = retained_protein_domains.split("|")
    partners = (("5'", five_prime), ("3'", three_prime))
    return ",".join(
        f"{partner}:{domain}"
        for partner, domains in partners
        for domain in domains.split(",")
        if domain.startswith("Protein_kinase_domain(")
    )


def in_frame_kinase_fusion(reading_frame, retained_protein_domains):
    """Return whether the fusion is in frame and retains a protein kinase domain.

    The result is ``yes``, ``no`` or an empty string when it is unknown: Arriba
    ran without protein domains or could not predict the reading frame.
    """
    if retained_protein_domains == ".":
        return ""
    if not retained_kinase_domain(retained_protein_domains):
        return "no"
    if reading_frame == ".":
        return ""
    return "yes" if reading_frame == "in-frame" else "no"


def event_flags(gene1, gene2, fusion_type):
    """Return the flags that follow from the genes and the type of the event."""
    flags = []
    if gene1 is None or gene2 is None:
        flags.append("intergenic")
    elif "ITD" in fusion_type:
        flags.append("itd")
    elif gene1 == gene2:
        flags.append("same_gene")
    if "read-through" in fusion_type:
        flags.append("read_through")
    if fusion_type.endswith(NO_CHIMERIC_PROTEIN_TYPES):
        flags.append("no_chimeric_protein")
    return flags


def knowledge_base_annotation(gene1, gene2, snapshot):
    """Return the CIViC and ChimerKB columns of a fusion of two genes."""
    annotation = {}
    civic, chimerkb = snapshot["civic"], snapshot["chimerkb"]
    forward, reciprocal = f"{gene1}::{gene2}", f"{gene2}::{gene1}"

    for match, key in (("exact", forward), ("reciprocal", reciprocal)):
        if key in civic:
            annotation["civic_match"] = match
            annotation["civic_best_levels"] = best_levels(civic[key])
            annotation["civic_evidence"] = summarize_evidence(civic[key])
            annotation["civic_pubmed_ids"] = pubmed_ids(civic[key])
            annotation["civic_nct_ids"] = nct_ids(civic[key])
            break

    # CIViC profiles such as "v::RET Fusion" hold for any fusion partner.
    any_partner = civic.get(f"{gene1}::v", []) + civic.get(f"v::{gene2}", [])
    annotation["civic_any_partner_profiles"] = ",".join(
        sorted({item["profile"] for item in any_partner})
    )
    annotation["civic_any_partner_best_levels"] = best_levels(any_partner)
    annotation["civic_any_partner_evidence"] = summarize_evidence(any_partner)
    annotation["civic_any_partner_pubmed_ids"] = pubmed_ids(any_partner)
    annotation["civic_any_partner_nct_ids"] = nct_ids(any_partner)

    record = chimerkb.get(forward) or chimerkb.get(reciprocal)
    if record:
        annotation["chimerkb_diseases"] = "; ".join(record["Disease"])
        annotation["chimerkb_validation"] = "; ".join(record["Validation"])
    return annotation


def annotate_fusion(row, snapshot):
    """Return the annotation columns of one Arriba fusion record."""
    annotation = dict.fromkeys(ANNOTATION_COLUMNS, "")
    gene1, gene2 = clean_gene(row["#gene1"]), clean_gene(row["gene2"])

    # Arriba writes "." when it ran without protein domains or when the
    # fusion retains no domain.
    if row["retained_protein_domains"] != ".":
        annotation["retained_kinase_domain"] = retained_kinase_domain(
            row["retained_protein_domains"]
        )
    annotation["in_frame_kinase_fusion"] = in_frame_kinase_fusion(
        row["reading_frame"], row["retained_protein_domains"]
    )

    flags = event_flags(gene1, gene2, row["type"])
    two_genes = gene1 is not None and gene2 is not None and gene1 != gene2

    if two_genes:
        normal_tissue = snapshot["normal_tissue"]
        normal_sources = normal_tissue.get(f"{gene1}::{gene2}") or normal_tissue.get(
            f"{gene2}::{gene1}"
        )
        if normal_sources:
            annotation["normal_tissue_sources"] = ",".join(normal_sources)
            flags.append("seen_in_normal")
    annotation["flags"] = ",".join(flags)

    # Knowledge base records describe fusions of two genes that can yield a
    # chimeric protein, so they do not apply to other events.
    if two_genes and "no_chimeric_protein" not in flags:
        annotation.update(knowledge_base_annotation(gene1, gene2, snapshot))
    return annotation


def summarize(rows, snapshot):
    """Return the snapshot versions and the match counts of annotated fusions."""
    return {
        "versions": snapshot["versions"],
        "fusions": len(rows),
        "civic_match": dict(Counter(row["civic_match"] or "none" for row in rows)),
        "civic_best_levels": dict(
            Counter(row["civic_best_levels"] or "none" for row in rows)
        ),
        "civic_any_partner": sum(
            bool(row["civic_any_partner_profiles"]) for row in rows
        ),
        "chimerkb": sum(
            bool(row["chimerkb_diseases"] or row["chimerkb_validation"]) for row in rows
        ),
        "in_frame_kinase_fusion": sum(
            row["in_frame_kinase_fusion"] == "yes" for row in rows
        ),
        "seen_in_normal": sum(bool(row["normal_tissue_sources"]) for row in rows),
        "flagged": sum(bool(row["flags"]) for row in rows),
    }


class FusionAnnotation(Process):
    """Add clinical annotation to gene fusions detected by Arriba.

    Each fusion keeps its Arriba columns and gets the following annotation
    from a knowledge base snapshot:

    - Clinical evidence from [CIViC](https://civicdb.org) for the fusion and
      for CIViC profiles that hold for any fusion partner (e.g. ``v::RET
      Fusion``): the best evidence level of each evidence type (A is the
      strongest, E the weakest), a summary of the supporting evidence with
      diseases and therapies, and the PubMed and ClinicalTrials.gov
      identifiers.
    - Diseases and validation methods from
      [ChimerKB](https://www.kobic.re.kr/chimerdb/).
    - The retained protein kinase domain and whether the fusion is an in-frame
      kinase fusion. This requires that Arriba ran with a protein domains file.
    - Sources that report the fusion in healthy tissue.
    - Flags that affect the interpretation: ``intergenic``, ``same_gene``,
      ``itd`` (internal tandem duplication), ``read_through``,
      ``no_chimeric_protein`` (genes fused head-to-head or tail-to-tail) and
      ``seen_in_normal``.

    Knowledge base records describe fusions of two genes that can yield a
    chimeric protein. Events flagged ``intergenic``, ``same_gene``, ``itd`` or
    ``no_chimeric_protein`` therefore get no CIViC or ChimerKB annotation.

    The snapshot is a JSON file with four sections. ``versions`` describes the
    snapshot and must include the key ``snapshot``. ``civic`` maps a fusion
    (``GENE1::GENE2``, 5' gene first, ``v`` for any partner) to a list of
    evidence items with the keys ``profile``, ``disease``, ``therapies``,
    ``therapy_interaction_type``, ``type``, ``direction``, ``level``,
    ``significance``, ``source_type``, ``citation_id`` and ``nct_ids``.
    ``chimerkb`` maps a fusion to the lists ``Disease`` and ``Validation``.
    ``normal_tissue`` maps a fusion to a list of sources. The process reads
    only the snapshot, so the results do not change until the snapshot
    changes.
    """

    slug = "fusion-annotation"
    name = "Gene fusion annotation"
    process_type = "data:genefusions:annotation"
    version = "1.0.0"
    category = "Gene fusions"
    scheduling_class = SchedulingClass.BATCH
    entity = {"type": "sample"}
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/common:5.1.0"}
        },
        "resources": {
            "cores": 1,
            "memory": 4096,
        },
    }
    data_name = '{{ fusions|name|default("?") }}'
    persistence = Persistence.CACHED

    class Input:
        """Input fields for FusionAnnotation process."""

        fusions = DataField(
            data_type="genefusions:arriba",
            label="Gene fusions detected by Arriba",
        )

        snapshot = DataField(
            data_type="file",
            label="Knowledge base snapshot",
            description="Knowledge base snapshot in JSON format.",
        )

    class Output:
        """Output fields for FusionAnnotation process."""

        annotated_fusions = FileField(label="Annotated fusions")
        summary = FileField(label="Annotation summary")
        snapshot_version = StringField(label="Knowledge base snapshot version")

    def run(self, inputs, outputs):
        """Annotate fusions with the knowledge base snapshot."""

        fusions_file = Path(inputs.fusions.output.fusions.path)
        name = fusions_file.name[: -len(".tsv.gz")]
        annotated_file = f"{name}_annotated.tsv.gz"
        summary_file = f"{name}_annotation_summary.json"

        try:
            snapshot = json.loads(Path(inputs.snapshot.output.file.path).read_text())
        except ValueError:
            self.error("The knowledge base snapshot is not a valid JSON file.")

        missing_sections = [
            section for section in SNAPSHOT_SECTIONS if section not in snapshot
        ]
        if missing_sections:
            self.error(
                "The knowledge base snapshot has no section: "
                f"{', '.join(missing_sections)}."
            )
        if "snapshot" not in snapshot["versions"]:
            self.error("The knowledge base snapshot has no version.")

        with gzip.open(fusions_file, "rt", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            rows = list(reader)
            columns = list(reader.fieldnames or [])

        missing_columns = [column for column in ARRIBA_COLUMNS if column not in columns]
        if missing_columns:
            self.error(f"The fusions file has no column: {', '.join(missing_columns)}.")

        for row in rows:
            row.update(annotate_fusion(row, snapshot))

        with gzip.open(annotated_file, "wt", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[*columns, *ANNOTATION_COLUMNS],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)

        Path(summary_file).write_text(
            json.dumps(summarize(rows, snapshot), indent=1) + "\n"
        )

        outputs.annotated_fusions = annotated_file
        outputs.summary = summary_file
        outputs.snapshot_version = snapshot["versions"]["snapshot"]
