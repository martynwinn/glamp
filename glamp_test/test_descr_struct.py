'''
Unittests for glamp/descr_struct.py

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
        self.epi_sequence=['MAAFMKLIQFLATKGQKYVSLAWKHKGTILKWINAGQSFEWIYKQIKKLWA']
        self.epi_offset=[0]
    
    def tearDown(self):
        '''
        Run after tests
        '''
        pass

    def test_testseq_score_structure(self):
        d = StructuralDescriptor(self.epi_sequence,self.epi_offset)
        scores = d.testseq_score_structure()
        self.assertEqual(len(scores[0]),12)
        self.assertEqual(scores[0][4],0.3)

if __name__ == '__main__':
    unittest.main()
