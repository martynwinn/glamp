#
"""
.. codeauthor:: Martyn Winn
"""

import numpy as np

class StructuralDescriptor():
    """
    Base class for structure-based peptide descriptors
    """
    
    def __init__(self,seqs,offsets):
        """
        Defaults ...
        """

        print("iggy")
        self.seqs = seqs
        self.offsets = offsets
        self.desc_struct = np.zeros((len(self.seqs),12))
        self.desc_core = np.zeros((len(self.seqs),7))
        self.desc_iface = np.zeros((len(self.seqs),4))
        self.desc_ener = np.zeros((len(self.seqs),8))
        self.num_descriptors = 0

    
    def testseq_score_structure(self):
        '''Score using known NI01 xtal structure.
        Structure is used in terms of which residues are surface or buried.
        Thus, can have not just number of Lys, but number of Lys on surface.'''
    
        # reference structure (move this out)
        # 0 = buried, 1 = intermediate, 2 = surface
        ref_structure = [2,2,1,0,2,2,0,0,2,1,
                    0,1,2,1,1,2,2,0,0,1,
                    0,0,2,2,1,1,1,1,0,1,
                    2,0,1,2,2,2,1,1,1,1,
                    1,0,0,2,1,0,2,2,2,1,
                    2]
    
        for iseq,in_seq in enumerate(self.seqs):
            # can't go beyond ref_structure description
            # TODO fix when have offset
            seq_length = min(51,len(in_seq))
            n_cationic_buried = 0
            n_anionic_buried = 0
            n_polar_buried = 0
            n_hphobic_buried = 0
            n_cationic_inter = 0
            n_anionic_inter = 0
            n_polar_inter = 0
            n_hphobic_inter = 0
            n_cationic_surface = 0
            n_anionic_surface = 0
            n_polar_surface = 0
            n_hphobic_surface = 0

            for i in range(seq_length):
                if ref_structure[i+self.offsets[iseq]] == 0:
                    if in_seq[i] in ['K','R','H']:
                        n_cationic_buried += 1
                    if in_seq[i] in ['D','E']:
                        n_anionic_buried += 1
                    if in_seq[i] in ['T','S','N','Q']:
                        n_polar_buried += 1
                    if in_seq[i] in ['A','V','L','I','F']:
                        n_hphobic_buried += 1
                if ref_structure[i+self.offsets[iseq]] == 1:
                    if in_seq[i] in ['K','R','H']:
                        n_cationic_inter += 1
                    if in_seq[i] in ['D','E']:
                        n_anionic_inter += 1
                    if in_seq[i] in ['T','S','N','Q']:
                        n_polar_inter += 1
                    if in_seq[i] in ['A','V','L','I','F']:
                        n_hphobic_inter += 1
                if ref_structure[i+self.offsets[iseq]] == 2:
                    if in_seq[i] in ['K','R','H']:
                        n_cationic_surface += 1
                    if in_seq[i] in ['D','E']:
                        n_anionic_surface += 1
                    if in_seq[i] in ['T','S','N','Q']:
                        n_polar_surface += 1
                    if in_seq[i] in ['A','V','L','I','F']:
                        n_hphobic_surface += 1

            self.desc_struct[iseq][0] = n_cationic_buried / 10.0
            self.desc_struct[iseq][1] = n_anionic_buried / 10.0
            self.desc_struct[iseq][2] = n_polar_buried / 10.0
            self.desc_struct[iseq][3] = n_hphobic_buried / 10.0
            self.desc_struct[iseq][4] = n_cationic_inter / 10.0
            self.desc_struct[iseq][5] = n_anionic_inter / 10.0
            self.desc_struct[iseq][6] = n_polar_inter / 10.0
            self.desc_struct[iseq][7] = n_hphobic_inter / 10.0
            self.desc_struct[iseq][8] = n_cationic_surface / 10.0
            self.desc_struct[iseq][9] = n_anionic_surface / 10.0
            self.desc_struct[iseq][10] = n_polar_surface / 10.0
            self.desc_struct[iseq][11] = n_hphobic_surface / 10.0
                      
        return self.desc_struct
   
    def testseq_score_core_res(self):
        '''Score using known NI01 xtal structure.
        Look at which residue types form the hydrophobic core.'''
    
        # reference structure (move this out)
        # residue numbers for core residues (starting at zero)
        # list agrees with 0=buried above, with addition of Ile33
        core_res = [3,6,7,10,17,18,20,21,28,31,32,41,42,45]
    
        for iseq,in_seq in enumerate(self.seqs):
            # TODO fix when have offset
            n_val = 0
            n_leu = 0
            n_ile = 0
            n_ala = 0
            n_phe = 0
            n_tyr = 0
            n_trp = 0

            for i in core_res:
                if in_seq[i] in ['V']:
                    n_val += 1
                if in_seq[i] in ['L']:
                    n_leu += 1
                if in_seq[i] in ['I']:
                    n_ile += 1
                if in_seq[i] in ['A']:
                    n_ala += 1
                if in_seq[i] in ['F']:
                    n_phe += 1
                if in_seq[i] in ['Y']:
                    n_tyr += 1
                if in_seq[i] in ['W']:
                    n_trp += 1

            # features for this sequence, normalise to 14 core residues
            self.desc_core[iseq][0] = n_val / 14.0
            self.desc_core[iseq][1] = n_leu / 14.0
            self.desc_core[iseq][2] = n_ile / 14.0
            self.desc_core[iseq][3] = n_ala / 14.0
            self.desc_core[iseq][4] = n_phe / 14.0
            self.desc_core[iseq][5] = n_tyr / 14.0
            self.desc_core[iseq][6] = n_trp / 14.0

        return self.desc_core
   
    def testseq_score_iface_res(self):
        '''Score using known NI01 xtal structure.
        Look at which residue types form the putative interface.'''
    
        # reference structure (move this out)
        # residue numbers for interface residues (starting at zero)
        iface_res = [23,24,27,30,31,34,36,40,44,47]
    
        for iseq,in_seq in enumerate(self.seqs):
            # TODO fix when have offset
            n_pos = 0
            n_neg = 0
            n_pol = 0
            n_hph = 0

            for i in iface_res:
                if in_seq[i] in ['K','R','H']:
                    n_pos += 1
                if in_seq[i] in ['D','E']:
                    n_neg += 1
                if in_seq[i] in ['T','S','N','Q']:
                    n_pol += 1
                if in_seq[i] in ['A','V','L','I','F','W']:
                    n_hph += 1

            # features for this sequence, normalise to 10 interface residues
            self.desc_iface[iseq][0] = n_pos / 10.0
            self.desc_iface[iseq][1] = n_neg / 10.0
            self.desc_iface[iseq][2] = n_pol / 10.0
            self.desc_iface[iseq][3] = n_hph / 10.0

        return self.desc_iface
            
    def testseq_score_energetics(self):
        '''Score input sequence according to deltaG from lookup table of
        http://jgp.rupress.org/content/129/5/371 (Tieleman 2007). The latter
        compares these values to the experimentally determined Wimley-White (1996)
        scales for deltaG to POPC interface and octanol. These can be got from 
        https://blanco.biomol.uci.edu/hydrophobicity_scales.html
        Modlamp has several hydrophobicity scales, but I assume these don't 
        consider lipid interfaces.

        Returns del_fe_interface, del_fe_center for input sequence in_seq
        Note 14/6/19: changed from deltadeltaG w.r.t. template to deltaG
        Note 13/7/20: updated to separate values for 4 helices'''
    
        water_to_interface = {'L': -14.1,
                         'I': -20.6,
                         'V': -12.2,
                         'F': -14.9,
                         'A': -6.8,
                         'W': -21.6,
                         'M': -10.5,
                         'C': -6.6,
                         'Y': -14.0,
                         'T': -4.2,
                         'S': -0.7,
                         'Q': -8.9,
                         'K': -18.6,
                         'N': -6.5,
                         'E': -1.68,
                         'D': 1.6,
                         'R': -21.2,
                         'H': 0.0,
                         'G': 0.0,
                         'P': 0.0}
        water_to_center = {'L': -15.2,
                         'I': -22.1,
                         'V': -13.8,
                         'F': -12.8,
                         'A': -8.4,
                         'W': -4.9,
                         'M': -4.4,
                         'C': -3.4,
                         'Y': 6.6,
                         'T': 13.9,
                         'S': 15.8,
                         'Q': 20.2,
                         'K': 19.9,
                         'N': 23.9,
                         'E': 21.1,
                         'D': 31.0,
                         'R': 58.1,
                         'H': 0.0,
                         'G': 0.0,
                         'P': 0.0}

        helix_boundaries = [[0,13],[15,25],[27,34],[36,50]]
        
        for iseq,in_seq in enumerate(self.seqs):

            seq_length = len(in_seq)
            helix_boundaries[3][1] = seq_length
            # accumulate separately for each helix in 4-helical bundle
            for h in range(4):
                del_fe_interface = 0.0
                del_fe_center = 0.0
                for i in range(helix_boundaries[h][0],helix_boundaries[h][1]):
                    del_fe_interface += water_to_interface[in_seq[i]]
                    del_fe_center += water_to_center[in_seq[i]]
        
                self.desc_ener[iseq][h*2+0] = del_fe_interface / 200.0
                self.desc_ener[iseq][h*2+1] = del_fe_center / 200.0

        return self.desc_ener
