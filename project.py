import Bio.PDB as pdb
from Bio import pairwise2, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB import Selection 
from Bio.PDB.Polypeptide import is_aa 
#from Bio.SubsMat.MatrixInfo import blosum62
from Bio.pairwise2 import format_alignment


import os
import copy as cp
import string




pdb_parser = pdb.PDBParser(PERMISSIVE=True, QUIET=True)
m1 = pdb_parser.get_structure("m1", "./example1/pair_his3_sc_XA.pdb.chainA.pdb")
m2 = pdb_parser.get_structure("m2", "./example1/pair_his3_sc_XA.pdb.chainB.pdb")

def get_fasta_alignment2(m1, m2):
    ppb = pdb.CaPPBuilder()
    for polypeptide in ppb.build_peptides(m1):
        sequence_ref = polypeptide.get_sequence()
    for polypeptide in ppb.build_peptides(m2):
        sequence_sample = polypeptide.get_sequence()
    align = pairwise2.align.globalxx(sequence_ref, sequence_sample)
    return list(align[0])
    for ref, sample, score, begin, end in align:
        reference_string=str(ref)
        sample_string=str(sample)
        filename = "alignment_ref_sample.fasta"
        with open(filename, "w") as handle:
            fasta_file = handle.write(">A\n%s\n>B\n%s\n" % (ref, sample)) 
            #print (fasta_file)

class StructureAlignment (object):
    from Bio.Data import SCOPData
    def __init__(self, align, m1, m2):
        """
         Attributes: 
           - fasta_align - Alignment object 
           - m1, m2 - two models 
           - si, sj - the sequences in the Alignment object that 
             correspond to the structures 
        """
        length = align[4]-align[3]
        # Get the residues in the models 
        rl1 = Selection.unfold_entities(m1, 'R') 
        rl2 = Selection.unfold_entities(m2, 'R') 
          # Residue positions 
        p1 = 0 
        p2 = 0 
          # Map equivalent residues to each other 
        map12 = {} 
        map21 = {} 
        residue_pairs = []
        for i in range(length):
            aa1 = align[0][i]
            aa2 = align[1][i]
            if aa1 != "-":
                while True:
                    r1 = rl1[p1]
                    p1 = p1 + 1
                    if is_aa(r1):
                        break
                self._test_equivalence(r1, aa1)
            else:
                r1 = None
            if aa2 != "-":
                while True:
                    r2 = rl2[p2]
                    p2 = p2 +1 
                    if is_aa(r2):
                        break
                self._test_equivalence(r2, aa2)
            else:
                r2 = None
            if r1: 
                map12[r1] = r2
            if r2:
                map21[r2] = r1

            residue_pairs.append((r1,r2))
        self.map12 = map12
        self.map21 = map21
        self.residue_pairs= residue_pairs
    def _test_equivalence(self, r1, aa1): 
        """Test if aa in sequence fits aa in structure (PRIVATE).""" 
        resname = r1.get_resname() 
        resname = SCOPData.protein_letters_3to1[resname] 
        assert(aa1 == resname)
    def get_maps(self): 
          """Map residues between the structures. 
   
          Return two dictionaries that map a residue in one structure to 
          the equivealent residue in the other structure. 
          """ 
        return self.map12, self.map21
    def get_iterator(self): 
         """Create an iterator over all residue pairs.""" 
        for i in range(0, len(self.residue_pairs)): 
            yield self.residue_pairs[i]



A=get_fasta_alignment2(m1,m2)
print(A)
print(A[0][1])
print(StructureAlignment(A,m1,m2).residue_pairs)





