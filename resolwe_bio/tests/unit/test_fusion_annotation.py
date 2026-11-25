from django.test import SimpleTestCase

from resolwe_bio.processes.fusion_calling.fusion_annotation import (
    annotate_fusion,
    best_levels,
    clean_gene,
    in_frame_kinase_fusion,
    retained_kinase_domain,
    summarize_evidence,
)

# Domain strings follow the Arriba output: 5' domains | 3' domains.
BCR_ABL1_DOMAINS = (
    "Bcr-Abl_oncoprotein_oligomerisation_domain(100%),C2_domain(100%)"
    "|Protein_kinase_domain(100%),SH2_domain(100%)"
)


def evidence(profile, level, evidence_type="Predictive", **fields):
    item = {
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


def fusion(gene1, gene2, fusion_type="translocation", frame="in-frame", domains="."):
    return {
        "#gene1": gene1,
        "gene2": gene2,
        "type": fusion_type,
        "reading_frame": frame,
        "retained_protein_domains": domains,
    }


SNAPSHOT = {
    "versions": {"snapshot": "test"},
    "civic": {
        "KIF5B::RET": [
            evidence(
                "KIF5B::RET Fusion",
                "C",
                therapies="Everolimus,Vandetanib",
                therapy_interaction_type="Combination",
                citation_id="900",
                nct_ids="NCT00000002, NCT00000001",
            ),
            evidence(
                "KIF5B::RET Fusion", "C", therapies="Vandetanib", citation_id="88"
            ),
            evidence(
                "KIF5B::RET Fusion",
                "B",
                "Diagnostic",
                significance="Positive",
                source_type="ASCO",
                citation_id="7",
            ),
            evidence(
                "KIF5B::RET Fusion",
                "A",
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
                therapies="Selpercatinib,Pralsetinib",
                therapy_interaction_type="Substitutes",
            )
        ],
        "KMT2A::v": [evidence("KMT2A::v Fusion", "A", "Diagnostic")],
    },
    "chimerkb": {
        "KIF5B::RET": {"Disease": ["NSCLC"], "Validation": ["FISH", "RT-PCR"]}
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
        annotation = annotate_fusion(fusion("KIF5B", "RET"), SNAPSHOT)

        self.assertEqual(annotation["civic_pubmed_ids"], "88,900")
        self.assertEqual(annotation["civic_nct_ids"], "NCT00000001,NCT00000002")

    def test_any_partner_evidence_is_reported_apart_from_the_exact_match(self):
        annotation = annotate_fusion(fusion("KIF5B", "RET"), SNAPSHOT)

        self.assertEqual(annotation["civic_match"], "exact")
        self.assertEqual(annotation["civic_best_levels"], "Diagnostic:B,Predictive:C")
        self.assertEqual(annotation["civic_any_partner_profiles"], "v::RET Fusion")
        self.assertEqual(annotation["civic_any_partner_best_levels"], "Predictive:A")
        self.assertEqual(annotation["chimerkb_diseases"], "NSCLC")

    def test_reciprocal_fusion_matches_and_is_labelled(self):
        annotation = annotate_fusion(fusion("RET", "KIF5B"), SNAPSHOT)

        self.assertEqual(annotation["civic_match"], "reciprocal")
        self.assertEqual(annotation["civic_any_partner_profiles"], "")
        self.assertEqual(annotation["chimerkb_validation"], "FISH; RT-PCR")

    def test_events_that_are_not_two_gene_fusions_get_no_knowledge_base_records(self):
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
                annotation = annotate_fusion(row, SNAPSHOT)

                self.assertEqual(annotation["flags"], flags)
                self.assertEqual(annotation["civic_match"], "")
                self.assertEqual(annotation["civic_any_partner_profiles"], "")
                self.assertEqual(annotation["chimerkb_diseases"], "")

    def test_fusion_seen_in_healthy_tissue_is_flagged(self):
        for gene1, gene2 in (("CTBS", "GNG5"), ("GNG5", "CTBS")):
            with self.subTest(gene1=gene1):
                annotation = annotate_fusion(
                    fusion(gene1, gene2, "deletion/read-through"), SNAPSHOT
                )

                self.assertEqual(
                    annotation["normal_tissue_sources"], "GTEx_recurrent_StarF2019"
                )
                self.assertEqual(annotation["flags"], "read_through,seen_in_normal")

    def test_fusion_without_findings_has_no_flags(self):
        self.assertEqual(annotate_fusion(fusion("KIF5B", "RET"), SNAPSHOT)["flags"], "")

    def test_retained_kinase_domain_names_the_partner(self):
        cases = {
            BCR_ABL1_DOMAINS: "3':Protein_kinase_domain(100%)",
            "Protein_kinase_domain(40%)|SH2_domain(100%)": (
                "5':Protein_kinase_domain(40%)"
            ),
            "SH2_domain(100%)|Diacylglycerol_kinase_catalytic_domain(100%)": "",
            "|": "",
        }
        for domains, kinase_domain in cases.items():
            with self.subTest(domains=domains):
                self.assertEqual(retained_kinase_domain(domains), kinase_domain)

    def test_in_frame_kinase_fusion(self):
        cases = [
            ("in-frame", BCR_ABL1_DOMAINS, "yes"),
            ("out-of-frame", BCR_ABL1_DOMAINS, "no"),
            ("stop-codon", BCR_ABL1_DOMAINS, "no"),
            ("in-frame", "SH2_domain(100%)|SH3_domain(100%)", "no"),
            # Unknown: no predicted reading frame, or Arriba ran without domains.
            (".", BCR_ABL1_DOMAINS, ""),
            ("in-frame", ".", ""),
        ]
        for reading_frame, domains, expected in cases:
            with self.subTest(reading_frame=reading_frame, domains=domains):
                self.assertEqual(
                    in_frame_kinase_fusion(reading_frame, domains), expected
                )
