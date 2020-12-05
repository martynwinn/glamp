'''
Unittests for glamp/sequence_sets.py

TODO check code coverage
'''
import unittest
from glamp.sequence_sets import SequenceSets

class SequenceSetsTest(unittest.TestCase):

    def setUp(self):
        '''
        Run before tests
        '''
        print("Running SequenceSetsTest ...")
        self.sequences1 = ['ACDEF','PQRST']
        self.sequences2 = ['AGDEF','PQRST']

    def tearDown(self):
        '''
        Run after tests
        '''
        pass

    def test_add_sequences(self):
        s = SequenceSets(self.sequences1)
        self.assertEqual(s.add_sequences(self.sequences2),4)

    def test_find_duplicates(self):
        s = SequenceSets(self.sequences1)
        s.add_sequences(self.sequences2)
        self.assertEqual(s.find_duplicates(),1)

if __name__ == '__main__':
    unittest.main()
