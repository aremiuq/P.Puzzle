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
import itertools
import settings


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

    def without_nones (self):
        """Return the maps without the none values.

        Input = map12, map21
        output = new_map12, new_map21
        here I want to remove all keys from both maps where the value = None
        a filter is created to remove nones, the original maps are cleared, updated, and renamed
        """

        #filter map 1
        new_map12 = {k: v for k,v in self.map12.items() if v is not None}
        #filter map 2
        new_map21 = {k: v for k,v in self.map21.items() if v is not None}
        return new_map12, new_map21

def Get_Pairwise(m1, m2):#error handle different posible dictionaris
    """
    input = model1 and model2
    output = pairwise alignment
    first, build peptides and get sequence for model1 and model2
    then, do a pairwise global alignment of both sequences
    """

    try:
        ppb = pdb.CaPPBuilder()
        for polypeptide in ppb.build_peptides(m1):
            sequence_ref = polypeptide.get_sequence()
        for polypeptide in ppb.build_peptides(m2):
            sequence_sample = polypeptide.get_sequence()
        align = pairwise2.align.globalxx(sequence_ref, sequence_sample)
        max_pair = max(align,key = lambda x:x[2])
        return max_pair
    except UnboundLocalError: #this error raises when the two compared chains are of diferent kind
        return None

def Check_Similarity(chain1, chain2, percent = 95):
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
    residue_chain_1, residue_chain_1_copy = itertools.tee(residue_chain_1)
    residue_chain_2 = chain2.get_residues()
    residue_chain_2, residue_chain_2_copy = itertools.tee(residue_chain_2)

    resnames_chain_1 = []
    for residue in residue_chain_1:
        resnames_chain_1.append(residue.get_resname())
    resnames_chain_1 = tuple(resnames_chain_1)
    resnames_chain_2 = []
    for residue in residue_chain_2:
        resnames_chain_2.append(residue.get_resname())
    resnames_chain_2 = tuple(resnames_chain_2)

    if resnames_chain_1 != resnames_chain_2:
        if (resnames_chain_1,resnames_chain_2) not in settings.similarity:
            align = Get_Pairwise(chain1, chain2)
            try:
                sim_percent = (align[2]/align[4]) * 100
            except TypeError as e:
                print("Diferent kinds of chains compared")
                return None
            if sim_percent <= percent:
                print("%s% of similarity between %s and %s"%(sim_percent, chain1.get_id(), chain1.get_id()))
                return None
            relations_1, relations_2 = StructureAlignment(align, chain1, chain2).without_nones()
            settings.similarity[(resnames_chain_1,resnames_chain_2)] = relations_1
            settings.similarity[(resnames_chain_2,resnames_chain_1)] = relations_2

        id_list_1 = [residue.get_id() for residue in settings.similarity[(resnames_chain_1, resnames_chain_2)].keys()]
        atom_list_1 = []
        for id in id_list_1:
            if "CA" in chain1[id]:
                print("Residue atoms added")
                atom_list_1.append(chain1[id]["CA"])
            else:
                print("DNA/RNA atoms added")
                atom_list_1.append(chain1[id]["C4'"])

        id_list_2 = [residue.get_id() for residue in settings.similarity[(resnames_chain_2, resnames_chain_1)].values()]
        atom_list_2 = []
        for id in id_list_2:
            if "CA" in chain2[id]:
                print("Residue atoms added")
                atom_list_2.append(chain2[id]["CA"])
            else:
                print("DNA/RNA atoms added")
                atom_list_2.append(chain2[id]["C4'"])
    else:
        atom_list_1 = []
        for residue in residue_chain_1_copy:
            atom_list_1.extend(residue.get_atoms())
        atom_list_2 = []
        for residue in residue_chain_2_copy:
            atom_list_2.extend(residue.get_atoms())

    return atom_list_1, atom_list_2

if __name__ == "__main__":

    pdb_parser = pdb.PDBParser(PERMISSIVE=True, QUIET=True)
    structure1 = pdb_parser.get_structure("structure1", "../example3/chain_X.pdb")
    structure2 = pdb_parser.get_structure("structure2", "../example3/chain_X_similar.pdb")

    for chain in structure1.get_chains():
        chain1 = chain
        for chain in structure2.get_chains():
            chain2 = chain

    print(check_residues_resnames(chain1, chain2, settings.similarity))

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
