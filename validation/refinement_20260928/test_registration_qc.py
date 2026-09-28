"""Behavior checks for a symmetric specimen with an anatomically wrong pose."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'skills/spateo-4d-pipeline/scripts'))
import numpy as np
from registration_qc import pair_qc

class AnatomicalQC(unittest.TestCase):
    def setUp(self):
        self.xyz=np.array([[x,y,z] for z in (-1,1) for x in (-1,1) for y in (-2,2)],float)
        self.labels=np.array(['ventral epidermal']*4+['dorsal epidermal']*4)
    def test_same_geometry_wrong_orientation_is_not_accepted_by_nn(self):
        rotated=self.xyz@np.diag([-1,1,-1])
        qc=pair_qc(self.xyz,rotated,self.labels,self.labels)
        self.assertEqual(qc['symmetric_nn_mean'],0)
        self.assertTrue(qc['dorsoventral_z']['opposite_signs'])
        self.assertEqual(qc['nearest_label_agreement_source'],0)
    def test_identical_specimens(self):
        qc=pair_qc(self.xyz,self.xyz,self.labels,self.labels)
        self.assertEqual(qc['symmetric_nn_mean'],0)
        self.assertFalse(qc['dorsoventral_z']['opposite_signs'])
        self.assertEqual(qc['nearest_label_agreement_source'],1)
    def test_species_without_planarian_labels(self):
        qc=pair_qc(self.xyz,self.xyz,['A']*8,['A']*8)
        self.assertNotIn('dorsoventral_z',qc)
        self.assertEqual(qc['annotations']['A']['source_n'],8)

if __name__=='__main__': unittest.main()
