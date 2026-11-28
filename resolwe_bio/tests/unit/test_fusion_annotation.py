import copy

from django.test import SimpleTestCase

from resolwe_bio.processes.fusion_calling.fusion_annotation import (
    annotate_fusion,
    annotate_fusions,
    best_levels,
    clean_gene,
    in_frame_kinase_fusion,
    normalize_snapshot,
    retained_kinase_domain,
    snapshot_error,
    summarize_evidence,
)

# Domain strings follow the Arriba output: 5' domains | 3' domains.
# BCR::ABL1 is from the test output of Arriba 2.4.0. The other strings use
# the domain names of the Arriba 2.4.0 protein domains database.
BCR_ABL1_DOMAINS = (
    "Bcr-Abl_oncoprotein_oligomerisation_domain(100%),C2_domain(100%)"
    "|Protein_kinase_domain(100%),SH2_domain(100%)"
)
NTRK1_DOMAINS = (
    "|Protein_kinase_domain(100%),"
    "Tyrosine_kinase_receptor_A_trans-membrane_domain(100%)"
)


def evidence(profile, level, evidence_type="Predictive", **fields):
    item = {
        "evidence_id": "1",
        "profile": profile,
        "disease": "Lung Adenocarcinoma",
        "therapies": "",
        "therapy_interaction_type": "",
        "type": evidence_type,
        "direction": "Supports",
        "level": level,
        "significance": "Sensitivity/Response",
        "source_type": "PubMed",
        "citation_id": "100",
        "nct_ids": "",
    }
    item.update(fields)
    return item


def fusion(
    gene1,
    gene2,
    fusion_type="translocation",
    frame="in-frame",
    domains=".",
    strand1="+/+",
    strand2="+/+",
):
    return {
        "#gene1": gene1,
        "gene2": gene2,
        "strand1(gene/fusion)": strand1,
        "strand2(gene/fusion)": strand2,
        "type": fusion_type,
        "reading_frame": frame,
        "retained_protein_domains": domains,
    }


def annotate(row, domains_reported=True):
    return annotate_fusion(row, SNAPSHOT, domains_reported)


SNAPSHOT = {
    "versions": {"snapshot": "test"},
    "civic": {
        "KIF5B::RET": [
            evidence(
                "KIF5B::RET Fusion",
                "C",
                evidence_id="698",
                therapies="Everolimus,Vandetanib",
                therapy_interaction_type="Combination",
                citation_id="900",
                nct_ids="NCT00000002, NCT00000001",
            ),
            evidence(
                "KIF5B::RET Fusion",
                "C",
                evidence_id="1622",
                therapies="Vandetanib",
                citation_id="88",
            ),
            evidence(
                "KIF5B::RET Fusion",
                "B",
                "Diagnostic",
                evidence_id="45",
                significance="Positive",
                source_type="ASCO",
                citation_id="7",
            ),
            evidence(
                "KIF5B::RET Fusion",
                "A",
                evidence_id="9",
                therapies="Cabozantinib",
                direction="Does Not Support",
                citation_id="5",
                nct_ids="NCT00000009",
            ),
        ],
        "v::RET": [
            evidence(
                "v::RET Fusion",
                "A",
                evidence_id="7066",
                therapies="Selpercatinib,Pralsetinib",
                therapy_interaction_type="Substitutes",
            )
        ],
        "KMT2A::v": [evidence("KMT2A::v Fusion", "A", "Diagnostic")],
    },
    "normal_tissue": {"CTBS::GNG5": ["GTEx_recurrent_StarF2019"]},
}


class FusionAnnotationTestCase(SimpleTestCase):
    def test_clean_gene(self):
        cases = {
            "kif5b": "KIF5B",
            "MS4A2(45793),LINC02705(15427)": None,
            "LINC02705(15427)": None,
            ".": None,
        }
        for field, symbol in cases.items():
            with self.subTest(field=field):
                self.assertEqual(clean_gene(field), symbol)

    def test_best_levels_are_reported_for_each_evidence_type(self):
        self.assertEqual(
            best_levels(SNAPSHOT["civic"]["KIF5B::RET"]), "Diagnostic:B,Predictive:C"
        )

    def test_best_levels_ignore_evidence_that_does_not_support(self):
        items = [evidence("A::B Fusion", "A", direction="Does Not Support")]

        self.assertEqual(best_levels(items), "")

    def test_evidence_summary_lists_claims_strongest_first(self):
        self.assertEqual(
            summarize_evidence(SNAPSHOT["civic"]["KIF5B::RET"]),
            "B: Lung Adenocarcinoma, Diagnostic, Positive"
            " | C: Lung Adenocarcinoma, Predictive, Sensitivity/Response,"
            " Everolimus+Vandetanib"
            " | C: Lung Adenocarcinoma, Predictive, Sensitivity/Response, Vandetanib",
        )

    def test_alternative_therapies_do_not_read_as_a_combination(self):
        self.assertEqual(
            summarize_evidence(SNAPSHOT["civic"]["v::RET"]),
            "A: Lung Adenocarcinoma, Predictive, Sensitivity/Response,"
            " Selpercatinib/Pralsetinib",
        )

    def test_ids_come_from_supporting_evidence_only(self):
        annotation = annotate(fusion("KIF5B", "RET"))

        self.assertEqual(annotation["civic_evidence_ids"], "45,698,1622")
        self.assertEqual(annotation["civic_pubmed_ids"], "88,900")
        self.assertEqual(annotation["civic_nct_ids"], "NCT00000001,NCT00000002")

    def test_any_partner_evidence_is_reported_apart_from_the_exact_match(self):
        annotation = annotate(fusion("KIF5B", "RET"))

        self.assertEqual(annotation["civic_match"], "exact")
        self.assertEqual(annotation["civic_best_levels"], "Diagnostic:B,Predictive:C")
        self.assertEqual(annotation["civic_any_partner_profiles"], "v::RET Fusion")
        self.assertEqual(annotation["civic_any_partner_best_levels"], "Predictive:A")
        self.assertEqual(annotation["civic_any_partner_evidence_ids"], "7066")

    def test_reciprocal_fusion_matches_and_is_labelled(self):
        annotation = annotate(fusion("RET", "KIF5B"))

        self.assertEqual(annotation["civic_match"], "reciprocal")
        self.assertEqual(annotation["civic_any_partner_profiles"], "")

    def test_unknown_transcribed_strand_matches_both_gene_orders(self):
        strands = [("+/.", "+/+"), ("+/+", "-/."), ("./.", "+/+")]
        for strand1, strand2 in strands:
            with self.subTest(strand1=strand1, strand2=strand2):
                annotation = annotate(
                    fusion("RET", "KIF5B", strand1=strand1, strand2=strand2)
                )

                self.assertEqual(annotation["civic_match"], "unknown_orientation")
                self.assertEqual(
                    annotation["civic_any_partner_profiles"], "v::RET Fusion"
                )

    def test_events_that_are_not_two_gene_fusions_get_no_civic_evidence(self):
        cases = {
            "no_chimeric_protein": fusion("KIF5B", "RET", "translocation/5'-5'"),
            "read_through,no_chimeric_protein": fusion(
                "KIF5B", "RET", "deletion/read-through/3'-3'"
            ),
            "same_gene": fusion("RET", "RET", "duplication"),
            "itd": fusion("KMT2A", "KMT2A", "duplication/ITD"),
            "intergenic": fusion("KIF5B(100),CCDC6(200)", "RET"),
        }
        for flags, row in cases.items():
            with self.subTest(flags=flags):
                annotation = annotate(row)

                self.assertEqual(annotation["flags"], flags)
                self.assertEqual(annotation["civic_match"], "")
                self.assertEqual(annotation["civic_any_partner_profiles"], "")

    def test_fusion_seen_in_healthy_tissue_is_flagged(self):
        for gene1, gene2 in (("CTBS", "GNG5"), ("GNG5", "CTBS")):
            with self.subTest(gene1=gene1):
                annotation = annotate(fusion(gene1, gene2, "deletion/read-through"))

                self.assertEqual(
                    annotation["normal_tissue_sources"], "GTEx_recurrent_StarF2019"
                )
                self.assertEqual(annotation["flags"], "read_through,seen_in_normal")

    def test_fusion_without_findings_has_no_flags(self):
        self.assertEqual(annotate(fusion("KIF5B", "RET"))["flags"], "")

    def test_retained_kinase_domain_names_the_partner(self):
        cases = {
            # ALK, ROS1 and RET fusions that retain the 3' kinase domain.
            "|Protein_kinase_domain(100%)": "3':Protein_kinase_domain(100%)",
            BCR_ABL1_DOMAINS: "3':Protein_kinase_domain(100%)",
            # The trans-membrane domain of NTRK1 is not a kinase domain.
            NTRK1_DOMAINS: "3':Protein_kinase_domain(100%)",
            "Protein_kinase_domain(40%)|SH2_domain(100%)": (
                "5':Protein_kinase_domain(40%)"
            ),
            "|Protein_tyrosine_kinase(100%)": "3':Protein_tyrosine_kinase(100%)",
            "SH2_domain(100%)|Diacylglycerol_kinase_catalytic_domain(100%)": "",
        }
        for domains, kinase_domain in cases.items():
            with self.subTest(domains=domains):
                self.assertEqual(retained_kinase_domain(domains), kinase_domain)

    def test_in_frame_kinase_fusion_requires_a_full_kinase_domain(self):
        cases = [
            ("in-frame", BCR_ABL1_DOMAINS, "yes"),
            ("in-frame", NTRK1_DOMAINS, "yes"),
            ("in-frame", "|Protein_tyrosine_kinase(100%)", "yes"),
            ("in-frame", "|Protein_kinase_domain(99%)", "no"),
            ("in-frame", "|Protein_kinase_domain(0%)", "no"),
            ("out-of-frame", BCR_ABL1_DOMAINS, "no"),
            ("stop-codon", BCR_ABL1_DOMAINS, "no"),
            ("in-frame", "SH2_domain(100%)|SH3_domain(100%)", "no"),
            # Arriba ran with protein domains and the fusion retains none.
            ("in-frame", ".", "no"),
            # Unknown: Arriba could not predict the reading frame.
            (".", BCR_ABL1_DOMAINS, ""),
        ]
        for reading_frame, domains, expected in cases:
            with self.subTest(reading_frame=reading_frame, domains=domains):
                self.assertEqual(
                    in_frame_kinase_fusion(reading_frame, domains, True), expected
                )

    def test_kinase_fusion_is_unknown_when_no_fusion_reports_domains(self):
        rows = [fusion("KIF5B", "RET"), fusion("EML4", "ALK")]

        annotate_fusions(rows, SNAPSHOT)

        self.assertEqual([row["in_frame_kinase_fusion"] for row in rows], ["", ""])

    def test_fusion_without_domains_is_no_kinase_fusion_when_others_report_them(self):
        rows = [fusion("KIF5B", "RET"), fusion("BCR", "ABL1", domains=BCR_ABL1_DOMAINS)]

        annotate_fusions(rows, SNAPSHOT)

        self.assertEqual([row["in_frame_kinase_fusion"] for row in rows], ["no", "yes"])

    def test_snapshot_keys_match_regardless_of_case(self):
        snapshot = copy.deepcopy(SNAPSHOT)
        snapshot["civic"] = {"C2orf44::Alk": [evidence("C2orf44::ALK Fusion", "B")]}

        annotation = annotate_fusion(
            fusion("C2orf44", "ALK"), normalize_snapshot(snapshot), True
        )

        self.assertEqual(annotation["civic_match"], "exact")

    def test_normalized_snapshot_keeps_the_any_partner_placeholder(self):
        self.assertEqual(
            sorted(normalize_snapshot(SNAPSHOT)["civic"]),
            ["KIF5B::RET", "KMT2A::v", "v::RET"],
        )

    def test_valid_snapshot_has_no_error(self):
        self.assertIsNone(snapshot_error(SNAPSHOT))

    def test_malformed_snapshot_gives_an_error(self):
        def changed(section, value):
            snapshot = copy.deepcopy(SNAPSHOT)
            snapshot[section] = value
            return snapshot

        no_civic = copy.deepcopy(SNAPSHOT)
        del no_civic["civic"]
        cases = {
            "is not a JSON object": [],
            "has no section: civic": no_civic,
            "The section 'versions' is not an object": changed("versions", "snapshot"),
            "has no version": changed("versions", {"snapshot": ""}),
            "The CIViC evidence of A::B is not a list": changed("civic", {"A::B": {}}),
            "has no text value for 'nct_ids'": changed(
                "civic", {"A::B": [evidence("A::B Fusion", "A", nct_ids=None)]}
            ),
            "has a PubMed citation_id that is not a number": changed(
                "civic", {"A::B": [evidence("A::B Fusion", "A", citation_id="")]}
            ),
            "has an evidence_id that is not a number": changed(
                "civic", {"A::B": [evidence("A::B Fusion", "A", evidence_id="EID1")]}
            ),
            "The healthy-tissue sources of A::B are not a list": changed(
                "normal_tissue", {"A::B": "GTEx"}
            ),
        }
        for message, snapshot in cases.items():
            with self.subTest(message=message):
                self.assertIn(message, snapshot_error(snapshot))
