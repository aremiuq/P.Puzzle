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
structure1 = pdb_parser.get_structure("structure1", "../example1/test.pdb")
structure2 = pdb_parser.get_structure("structure2", "../example1/pair_his3_sc_XA.pdb.chainB.pdb")

chain1 = structure1.get_chains()
chain2 = structure2.get_chains()



def check_residues_resnames(chain1, chain2):
    """
    input = chain1, chain2
    output = dictionary with key of (residue_chain_1, residue_chain_2)
    from these chains, we get the residues and residue names
    first check and see if residue names between chains are equal
        if equal then return output
    second, if residue names between chains are not equal and if the tuple is not in dictionary
        do pairwise align, structural align, and add to dictionary 
        return dictionary 
    
    """
    residue_chain_1 = chain1.get_residues()
    residue_chain_2 = chain2.get_residues()
    resnames_chain_1 = residue_chain_1.get_resname()
    resnames_chain_2 = residue_chain_2.get_resname()

    if resnames_chain_1 == resnames_chain_2:
        return dict[(residues_chain_1, residues_chain_2)]  #superimpose

    if resnames_chain_1 != resnames_chain_2:
        if (residues_chain_1,residues_chain_2) not in dict:
            align = get_fasta_alignment2(residue_chain_1, residue_chain_2)
            stuct_align = align_structure(align, residue_chain_1, residue_chain_2)
            dict[(residue_chain_1,residue_chain_2)] = relations_1
            dict[(residue_chain_2,residue_chain_1)] = relations_2
        
        return dict[(residues_chain_1, residues_chain_2)]


def get_pairwise_alignment(m1, m2):
    """
    input = model1 and model2
    output = pairwise alignment
    first, build peptides and get sequence for model1 and model2 
    then, do a pairwise global alignment of both sequences 
    """
    ppb = pdb.CaPPBuilder()
    for polypeptide in ppb.build_peptides(m1):
        sequence_ref = polypeptide.get_sequence()
    for polypeptide in ppb.build_peptides(m2):
        sequence_sample = polypeptide.get_sequence()
    align = pairwise2.align.globalxx(sequence_ref, sequence_sample)
    return align
    #for ref, sample, score, begin, end in align:
     #   reference_string=str(ref)
      #  sample_string=str(sample)
       # filename = "alignment_ref_sample.fasta"
       # with open(filename, "w") as handle:
        #    fasta_file = handle.write(">A\n%s\n>B\n%s\n" % (ref, sample))
         #   print (fasta_file)


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

    def remove_nones (map12, map21):
        """
        input = map12, map21 
        output = new_map12, new_map21
        here I want to remove all keys from both maps where the value = None
        a filter is created to remove nones, the original maps are cleared, updated, and renamed 
        """
        #return len(second_map) 
        #filter map 1
        filter_m1_none = {k: v for k,v in map12.items() if v is not None}
        map12.clear()
        new_map12 = map12.update(filter_m1_none)
        #return len(first_map)
        #filter map 2
        filter_m2_none = {k: v for k,v in map21.items() if v is not None}
        map21.clear()
        new_map21 = map21.update(filter_m2_none)
        #return len(second_map)
        return new_map12, new_map21


#pairwise = pairwise("P", "AX", "AB", {'AX': {'A': 'object', 'X': 'object', ('P', 'AB'): "relationship"}, 'AB': {'A': 'object', 'B': 'object'}})
#print(pairwise)
#print (check_residues_resnames)

#A=get_fasta_alignment2(m1,m2)
#print(A)
#print(A[0][1])
#aligment_map = StructureAlignment(A,m1,m2).get_maps()
#print((aligment_map[0]))
#print((aligment_map[1]))
#print(len(aligment_2))
#print(list(aligment_m1_map.keys())[0].get_id())









