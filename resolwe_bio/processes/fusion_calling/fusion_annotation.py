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
    "strand1(gene/fusion)",
    "strand2(gene/fusion)",
    "type",
    "reading_frame",
    "retained_protein_domains",
)

SNAPSHOT_SECTIONS = ("versions", "civic", "normal_tissue")

EVIDENCE_KEYS = (
    "evidence_id",
    "profile",
    "disease",
    "therapies",
    "therapy_interaction_type",
    "type",
    "direction",
    "level",
    "significance",
    "source_type",
    "citation_id",
    "nct_ids",
)

ANNOTATION_COLUMNS = (
    "civic_match",
    "civic_best_levels",
    "civic_evidence",
    "civic_evidence_ids",
    "civic_pubmed_ids",
    "civic_nct_ids",
    "civic_any_partner_profiles",
    "civic_any_partner_best_levels",
    "civic_any_partner_evidence",
    "civic_any_partner_evidence_ids",
    "civic_any_partner_pubmed_ids",
    "civic_any_partner_nct_ids",
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

# Names of protein kinase domains in the retained_protein_domains column.
# The Arriba protein domains database names them "Protein kinase domain"
# (Pfam PF00069). Other annotations also use "Protein tyrosine kinase"
# (Pfam PF07714). Arriba writes the names with underscores.
KINASE_DOMAIN_NAMES = ("Protein_kinase_domain", "Protein_tyrosine_kinase")


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


def normalize_key(key):
    """Return a snapshot fusion key with upper-case gene symbols.

    The ``v`` placeholder for any fusion partner keeps its case.
    """
    return "::".join(
        gene if gene == "v" else gene.strip().upper() for gene in key.split("::")
    )


def evidence_item_error(item):
    """Return what is wrong with a CIViC evidence item, or ``None``."""
    if not isinstance(item, dict):
        return "is not an object"
    for key in EVIDENCE_KEYS:
        if not isinstance(item.get(key), str):
            return f"has no text value for '{key}'"
    if not item["evidence_id"].isdigit():
        return "has an evidence_id that is not a number"
    if item["source_type"] == "PubMed" and not item["citation_id"].isdigit():
        return "has a PubMed citation_id that is not a number"
    return None


def section_error(name, section):
    """Return what is wrong with a snapshot section, or ``None``."""
    if not isinstance(section, dict):
        return f"The section '{name}' is not an object."
    for key, value in section.items():
        if name == "civic":
            if not isinstance(value, list):
                return f"The CIViC evidence of {key} is not a list."
            for item in value:
                error = evidence_item_error(item)
                if error:
                    return f"A CIViC evidence item of {key} {error}."
        elif name == "normal_tissue" and not isinstance(value, list):
            return f"The healthy-tissue sources of {key} are not a list."
    return None


def snapshot_error(snapshot):
    """Return what is wrong with a knowledge base snapshot, or ``None``."""
    if not isinstance(snapshot, dict):
        return "The knowledge base snapshot is not a JSON object."
    missing_sections = [
        section for section in SNAPSHOT_SECTIONS if section not in snapshot
    ]
    if missing_sections:
        return (
            "The knowledge base snapshot has no section: "
            f"{', '.join(missing_sections)}."
        )
    for name in SNAPSHOT_SECTIONS:
        error = section_error(name, snapshot[name])
        if error:
            return f"The knowledge base snapshot is not valid. {error}"
    version = snapshot["versions"].get("snapshot")
    if not isinstance(version, str) or not version:
        return "The knowledge base snapshot has no version."
    return None


def normalize_snapshot(snapshot):
    """Return the snapshot with upper-case gene symbols in the fusion keys."""
    normalized = {"versions": snapshot["versions"]}
    for name in ("civic", "normal_tissue"):
        normalized[name] = {
            normalize_key(key): value for key, value in snapshot[name].items()
        }
    return normalized


def orientation_known(row):
    """Return whether Arriba predicted the transcribed strand of both genes.

    The strand columns hold the gene strand and the transcribed strand,
    separated by a slash. Arriba writes ``.`` for a transcribed strand that it
    cannot predict. The 5' gene is then unknown, and ``#gene1`` and ``gene2``
    are in arbitrary order.
    """
    return all(
        row[column].split("/")[-1] != "."
        for column in ("strand1(gene/fusion)", "strand2(gene/fusion)")
    )


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


def evidence_ids(evidence):
    """Return the sorted unique CIViC IDs of supporting evidence items."""
    ids = {item["evidence_id"] for item in evidence if item["direction"] == "Supports"}
    return ",".join(sorted(ids, key=int))


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
        if domain.rsplit("(", 1)[0] in KINASE_DOMAIN_NAMES
    )


def in_frame_kinase_fusion(reading_frame, retained_protein_domains, domains_reported):
    """Return whether the fusion is in frame and retains a full kinase domain.

    The result is ``yes``, ``no`` or an empty string when it is unknown.

    Only a fully retained protein kinase domain (``100%``) counts, because a
    truncated catalytic domain is probably not active. Arriba reports every
    domain with at least one retained base and rounds the fraction down.

    Arriba writes ``.`` in the retained protein domains column both when it
    ran without a protein domains file and when the fusion retains no domain.
    ``domains_reported`` tells whether any fusion of the sample has retained
    domains, which shows that Arriba ran with a protein domains file. The
    result is unknown when it is ``False`` or when Arriba could not predict
    the reading frame of a fusion with a full kinase domain.
    """
    if not domains_reported:
        return ""
    if retained_protein_domains == ".":
        return "no"
    kinase_domains = retained_kinase_domain(retained_protein_domains).split(",")
    if not any(domain.endswith("(100%)") for domain in kinase_domains):
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


def find_record(records, gene1, gene2, known_orientation):
    """Return how a fusion matches a knowledge base section, and the record.

    The match is ``exact`` for the record of ``gene1::gene2`` and
    ``reciprocal`` for the record of ``gene2::gene1``. When the 5' gene is
    unknown, a record in either order matches as ``unknown_orientation``.
    The result is ``("", None)`` without a record.
    """
    for match, key in (
        ("exact", f"{gene1}::{gene2}"),
        ("reciprocal", f"{gene2}::{gene1}"),
    ):
        if key in records:
            return (match if known_orientation else "unknown_orientation"), records[key]
    return "", None


def civic_annotation(gene1, gene2, known_orientation, snapshot):
    """Return the CIViC columns of a fusion of two genes."""
    annotation = {}
    civic = snapshot["civic"]

    match, evidence = find_record(civic, gene1, gene2, known_orientation)
    if evidence is not None:
        annotation["civic_match"] = match
        annotation["civic_best_levels"] = best_levels(evidence)
        annotation["civic_evidence"] = summarize_evidence(evidence)
        annotation["civic_evidence_ids"] = evidence_ids(evidence)
        annotation["civic_pubmed_ids"] = pubmed_ids(evidence)
        annotation["civic_nct_ids"] = nct_ids(evidence)

    # CIViC profiles such as "v::RET Fusion" hold for any fusion partner.
    # Either gene can be the 5' gene when Arriba did not predict it.
    any_partner_keys = [f"{gene1}::v", f"v::{gene2}"]
    if not known_orientation:
        any_partner_keys += [f"{gene2}::v", f"v::{gene1}"]
    any_partner = [item for key in any_partner_keys for item in civic.get(key, [])]
    annotation["civic_any_partner_profiles"] = ",".join(
        sorted({item["profile"] for item in any_partner})
    )
    annotation["civic_any_partner_best_levels"] = best_levels(any_partner)
    annotation["civic_any_partner_evidence"] = summarize_evidence(any_partner)
    annotation["civic_any_partner_evidence_ids"] = evidence_ids(any_partner)
    annotation["civic_any_partner_pubmed_ids"] = pubmed_ids(any_partner)
    annotation["civic_any_partner_nct_ids"] = nct_ids(any_partner)
    return annotation


def annotate_fusion(row, snapshot, domains_reported):
    """Return the annotation columns of one Arriba fusion record.

    ``domains_reported`` tells whether any fusion of the sample has retained
    protein domains.
    """
    annotation = dict.fromkeys(ANNOTATION_COLUMNS, "")
    gene1, gene2 = clean_gene(row["#gene1"]), clean_gene(row["gene2"])

    if row["retained_protein_domains"] != ".":
        annotation["retained_kinase_domain"] = retained_kinase_domain(
            row["retained_protein_domains"]
        )
    annotation["in_frame_kinase_fusion"] = in_frame_kinase_fusion(
        row["reading_frame"], row["retained_protein_domains"], domains_reported
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

    # CIViC fusion profiles describe fusions of two genes that can yield a
    # chimeric protein, so they do not apply to other events.
    if two_genes and "no_chimeric_protein" not in flags:
        annotation.update(
            civic_annotation(gene1, gene2, orientation_known(row), snapshot)
        )
    return annotation


def annotate_fusions(rows, snapshot):
    """Add the annotation columns to all Arriba fusion records of a sample."""
    domains_reported = any(row["retained_protein_domains"] != "." for row in rows)
    for row in rows:
        row.update(annotate_fusion(row, snapshot, domains_reported))


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
      diseases and therapies, and the CIViC evidence, PubMed and
      ClinicalTrials.gov identifiers.
    - The retained protein kinase domain and whether the fusion is an in-frame
      kinase fusion, which requires a fully retained kinase domain.
    - Sources that report the fusion in healthy tissue.
    - Flags that affect the interpretation: ``intergenic``, ``same_gene``,
      ``itd`` (internal tandem duplication), ``read_through``,
      ``no_chimeric_protein`` (genes fused head-to-head or tail-to-tail) and
      ``seen_in_normal``.

    CIViC fusion profiles describe fusions of two genes that can yield a
    chimeric protein. Events flagged ``intergenic``, ``same_gene``, ``itd`` or
    ``no_chimeric_protein`` therefore get no CIViC annotation.

    A CIViC profile matches as ``exact`` when its 5' and 3' genes are the
    genes of the fusion in the same order, and as ``reciprocal`` when they
    are in the reverse order. When Arriba cannot predict the transcribed
    strand, the 5' gene is unknown: a profile in either order then matches as
    ``unknown_orientation``, and profiles for any fusion partner are looked up
    for both genes.

    The kinase annotation requires that Arriba ran with a protein domains
    file that names kinase domains ``Protein kinase domain`` or ``Protein
    tyrosine kinase``, as the Arriba protein domains database does. It is
    empty when no fusion of the sample has retained domains, because Arriba
    then most likely ran without a protein domains file.

    The snapshot is a JSON file with three sections:

    - ``versions`` describes the snapshot. It must include the key
      ``snapshot`` with the snapshot version and should record the version
      and the license of each source.
    - ``civic`` maps a fusion to a list of evidence items with the text keys
      ``evidence_id``, ``profile``, ``disease``, ``therapies``,
      ``therapy_interaction_type``, ``type``, ``direction``, ``level``,
      ``significance``, ``source_type``, ``citation_id`` and ``nct_ids``.
    - ``normal_tissue`` maps a fusion to a list of sources.

    A fusion key is ``GENE1::GENE2`` with the 5' gene first and ``v`` for any
    fusion partner. Gene symbols match regardless of case. The process cannot
    check the following rules, so the snapshot builder has to follow them:

    - Include only accepted CIViC evidence.
    - Include only CIViC profiles that describe a fusion alone. Leave out
      profiles that combine a fusion with another variant.
    - Include only sources whose license allows the intended use.

    The process reads only the snapshot, so the results do not change until
    the snapshot changes.
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

        error = snapshot_error(snapshot)
        if error:
            self.error(error)
        snapshot = normalize_snapshot(snapshot)

        with gzip.open(fusions_file, "rt", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            rows = list(reader)
            columns = list(reader.fieldnames or [])

        missing_columns = [column for column in ARRIBA_COLUMNS if column not in columns]
        if missing_columns:
            self.error(f"The fusions file has no column: {', '.join(missing_columns)}.")

        annotate_fusions(rows, snapshot)

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
