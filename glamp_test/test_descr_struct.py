'''
Unittests for glamp/descr_struct.py

Run as:
export PYTHONPATH=/home/mdw/ingenza/machine_learning/amps
python test_descr_struct.py

TODO check code coverage
'''
import unittest
from glamp.descr_struct import StructuralDescriptor

class StructuralDescriptorTest(unittest.TestCase):

    def setUp(self):
        '''
        Run before tests
        '''
        print("Running StructuralDescriptorTest ...")
        self.epi_sequence='MAAFMKLIQFLATKGQKYVSLAWKHKGTILKWINAGQSFEWIYKQIKKLWA'
        self.aucn_sequence='MSWLNFLKYIAKYGKKAVSAAWKYKGKVLEWLNVGPTLEWVWQKLKKIAGL'
        self.epi_offset=0
        self.aucn_offset=1
    
    def tearDown(self):
        '''
        Run after tests
        '''
        pass

    def test_get_residue_location(self):
        d = StructuralDescriptor()
        # Ile42 of NI01 should be well buried
        self.assertEqual(d.get_residue_location(42,0),'buried')
        # Asn34 of NI01 should be on surface
        self.assertEqual(d.get_residue_location(34,0),'surface')

        # Leu45 of A53 should be well buried
        self.assertEqual(d.get_residue_location(45,1),'buried')
        # Glu39 of A53 should be on surface
        self.assertEqual(d.get_residue_location(39,1),'surface')
    
    def test_testseq_score_structure(self):
        d = StructuralDescriptor([self.epi_sequence,self.aucn_sequence],
                                 [self.epi_offset,self.aucn_offset])
        scores = d.testseq_score_structure()
        self.assertEqual(len(scores[0]),12)
        # check the NI01 scores
        self.assertEqual(scores[0][6],0.5)
        # check the Aucn scores
        self.assertEqual(scores[1][6],0.3)
    
    def test_testseq_score_core_res(self):
        d = StructuralDescriptor([self.epi_sequence,self.aucn_sequence],
                                 [self.epi_offset,self.aucn_offset])
        scores = d.testseq_score_core_res()
        self.assertEqual(len(scores[0]),7)
        # check the NI01 scores
        self.assertTrue(scores[0][1] > 0.230 and scores[0][1] < 0.231,
                        msg='core res score for NI01 wrong')
        # check the Aucn scores
        self.assertTrue(scores[1][1] > 0.307 and scores[1][1] < 0.308,
                        msg='core res score for Aucn wrong')

    def test_testseq_score_iface_res(self):
        d = StructuralDescriptor([self.epi_sequence,self.aucn_sequence],
                                 [self.epi_offset,self.aucn_offset])
        scores = d.testseq_score_iface_res()
        self.assertEqual(len(scores[0]),4)
        # check the NI01 scores
        self.assertTrue(scores[0][3] > 0.29 and scores[0][3] < 0.31,
                        msg='iface res score for NI01 wrong')
        # check the Aucn scores
        self.assertTrue(scores[1][3] > 0.59 and scores[1][3] < 0.61,
                        msg='iface res score for Aucn wrong')

    def test_testseq_score_energetics(self):
        d = StructuralDescriptor([self.epi_sequence,self.aucn_sequence],
                                 [self.epi_offset,self.aucn_offset])
        scores = d.testseq_score_energetics()
        self.assertEqual(len(scores[0]),8)
        # check the NI01 scores
        self.assertTrue(scores[0][3] > 0.200 and scores[0][3] < 0.201,
                        msg='energetics score for NI01 wrong')
        # check the Aucn scores
        self.assertTrue(scores[1][3] > 0.190 and scores[1][3] < 0.191,
                        msg='energetics score for Aucn wrong')

if __name__ == '__main__':
    unittest.main()
