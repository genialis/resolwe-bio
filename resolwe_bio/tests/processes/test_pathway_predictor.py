# pylint: disable=missing-docstring
import unittest


from resolwe.test import ProcessTestCase, tag_process


class PathwayPredictorTestCase(ProcessTestCase):
    @tag_process('pathways-predictor')
    def test_pathway_predictor(self):
        n_pathways = 3
        inputs = {'model': 'iMM904', 'product': 'methanol', 'n_pathways': n_pathways}
        annotation = self.run_process('pathways-predictor', inputs)
        self.assertEqual(annotation.output['pathways'], 3)
