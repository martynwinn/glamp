#
"""
.. codeauthor:: Martyn Winn
"""

import numpy as np
import difflib
import random

class SequenceSets():
    """
    Base class for sets of peptide sequences
    """
    
    def __init__(self,seqs,offsets=None):
        """
        Class is a set of related sequences, grouped together so we
        can do some simple comparisons.
        """

        self.seqs = seqs
        if offsets is not None:
            self.offsets = np.array(offsets)
        else:
            self.offsets = np.zeros(len(seqs))
        print("Created sequence set with %i sequences" % (len(self.seqs)))

        # one member of the set can be designated the reference
        self.reference_seq = 0

    def add_sequences(self,seqs,offsets=None):
        """
        Add a set of sequences to the existing set.
        """

        self.seqs.extend(seqs)
        if offsets is not None:
            self.offsets = np.append(self.offsets,offsets)
        else:
            self.offsets = np.append(self.offsets,np.zeros(len(seqs)))
        print("Appended %i sequences to give %i in the set" % (len(seqs),len(self.seqs)))

    def find_duplicates(self):
        '''Check to see if there are any duplicated sequences.'''

        duplicates = []
        # this is not the most efficient O(n^2) but don't think that
        # matters in this context
        for s in self.seqs:
            if self.seqs.count(s) > 1:
                if s not in duplicates:
                    print("Sequence %s occurs %i times" % (s,self.seqs.count(s)))
                    duplicates.append(s)
                    print(' '.join(list(difflib.context_diff(self.seqs[self.reference_seq],s))))

        if len(duplicates) > 0:
            print("Duplicates found!")
        else:
            print("No duplicates.")

    def generate_sequence_path(self,start_seq,nser=10,lser=10,nmut=1,scope='cons'):
        '''Generate sequences using mutations from a starting sequence.
           nser: number of series
           lser: length of a series
           nmut: number of residues mutated each time
           scope: conservative or random mutations
           TODO check sequences in set are unique
           TODO set probability threshold based on scope parameter
        '''

        aa = {'L':'h','I':'h','V':'h','F':'h','A':'h','W':'h',
              'Y':'p','T':'p','S':'p','Q':'p','N':'p',
              'E':'a','D':'a',
              'K':'c','R':'c','H':'c',
              'M':'x','C':'x','G':'x','P':'x'}
        
        new_seqs = []
        for iser in range(nser):
            # for each series, reset the current sequence
            curr_seq = start_seq
            for iseq in range(lser):
                for m in range(nmut):
                    # select place to mutate
                    i = random.randint(0,50)
                    curr_type = aa[curr_seq[i]]
                    # selecting what to mutate into
                    not_mutated = True
                    while not_mutated:
                        target = random.choice(list(aa))
                        target_type = aa[target]
                        # conservative mutation, accept with probability 0.8
                        if target_type == curr_type:
                            if random.randint(0,10) < 8:
                                print('Conservative mutation of residue %i from %s to %s' % (i+1,curr_seq[i],target))
                                curr_seq = curr_seq[:i] + target + curr_seq[i+1:]
                                not_mutated = False
                        # non-conservative mutation, accept with probability 0.2
                        else:
                             if random.randint(0,10) < 2:
                                 print('Non-conservative mutation of residue %i from %s to %s' % (i+1,curr_seq[i],target))
                                 curr_seq = curr_seq[:i] + target + curr_seq[i+1:]
                                 not_mutated = False
                            
                # after nmut mutations, append the new sequence    
                new_seqs.append(curr_seq)

        return new_seqs
                
