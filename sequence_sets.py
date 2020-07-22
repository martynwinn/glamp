#
"""
.. codeauthor:: Martyn Winn
"""

import numpy as np
import difflib

class SequenceSets():
    """
    Base class for sets of peptide sequences
    """
    
    def __init__(self,seqs,offsets=None):
        """
        Defaults ...
        """

        print("oogle")
        self.seqs = seqs
        if offsets is not None:
            self.offsets = np.array(offsets)
        else:
            self.offsets = np.zeros(len(seqs))
        print("Created sequence set with %i sequences" % (len(self.seqs)))

        self.reference_seq = 0

    def add_sequences(self,seqs,offsets=None):
        """
        Defaults ...
        """

        self.seqs.extend(seqs)
        if offsets is not None:
            self.offsets = np.append(self.offsets,offsets)
        else:
            self.offsets = np.append(self.offsets,np.zeros(len(seqs)))
        print("Appended %i sequences to give %i in the set" % (len(seqs),len(self.seqs)))

    def find_duplicates(self):
        '''Check to see if there are any duplicated sequences.'''

        duplicates = False
        # this is not the most efficient O(n^2) but don't think that
        # matters in this context
        for s in self.seqs:
            if self.seqs.count(s) > 1:
                print("Sequence %s occurs %i times" % (s,self.seqs.count(s)))
                duplicates = True
                print(' '.join(list(difflib.context_diff(self.seqs[self.reference_seq],s))))

        if duplicates:
            print("Duplicates found!")
        else:
            print("No duplicates.")

    def generate_sequence_path(self,start_seq,cumulative=True,nseq=10):
        '''Generate sequences using single mutations from a starting
           sequence.'''

        new_seqs = []
        for iseq in range(nseq):
            i = randrange(51)
