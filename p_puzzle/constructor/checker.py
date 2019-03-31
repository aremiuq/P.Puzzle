#checker.py imports

from Bio.pairwise2 import format_alignment
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2, SeqIO
from Bio.PDB import Selection
from Bio.Seq import Seq
from ..support import settings as s
import Bio.PDB as pdb
import itertools


#checker.py classes

class StructureAlignment (object):
    from Bio.Data import SCOPData
    def __init__(self, align, m1, m2):
        """Produces a structural alignment of two models
        Input:
        - fasta_align - Alignment object
        - m1, m2 - two models
        - si, sj - the sequences in the Alignment object that correspond to the structures
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
          Output:
          map12, map21 = residue maps of structures
          """
        return self.map12, self.map21
    def get_iterator(self):
        """Create an iterator over all residue pairs."""
        for i in range(0, len(self.residue_pairs)):
            yield self.residue_pairs[i]

    def without_nones (self):
        """Return the maps without the none values.
        Input:
        -map12, map21 = residue maps of structures
        output:
        -new_map12, new_map21 = maps with "Nones" removed
        """

        #filter map 1
        new_map12 = {k: v for k,v in self.map12.items() if v is not None}
        #filter map 2
        new_map21 = {k: v for k,v in self.map21.items() if v is not None}
        return new_map12, new_map21

#checker.py functions

def Get_Pairwise(m1, m2):
    """Complete a pairwise global alignment of both model's sequence.
    Input
    -model1, model2
    Output:
    max_pair = return pair of sequences with highest alignment score
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
    """Checks similarity between two chains
    Input:
    -chain1, chain2 = chains to be checked
    Output:
    -atom_list_1, atom_list_2 = atoms of the chains being checked
    """

    s.argprint("Checking similarity...", s.options.verbose, s.options.quiet, 1)

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
        s.argprint("Different chains...", s.options.verbose, s.options.quiet, 1)
        if (resnames_chain_1,resnames_chain_2) not in s.similarity:
            s.argprint("First time of this difference", s.options.verbose, s.options.quiet, 2)
            align = Get_Pairwise(chain1, chain2)
            try:
                sim_percent = (align[2]/align[4]) * 100
            except TypeError as e:
                s.argprint("Different kinds of chains compared", s.options.verbose, s.options.quiet, 2)
                return None
            if sim_percent <= percent:
                s.argprint("%s% of similarity between %s and %s"%(sim_percent, chain1.get_id(), chain1.get_id()), s.options.verbose, s.options.quiet, 2)
                return None
            relations_1, relations_2 = StructureAlignment(align, chain1, chain2).without_nones()
            s.similarity[(resnames_chain_1,resnames_chain_2)] = relations_1
            s.similarity[(resnames_chain_2,resnames_chain_1)] = relations_2

        id_list_1 = [residue.get_id() for residue in s.similarity[(resnames_chain_1, resnames_chain_2)].keys()]
        atom_list_1 = []
        for id in id_list_1:
            if "CA" in chain1[id]:
                s.argprint("Residue atoms added", s.options.verbose, s.options.quiet, 2)
                atom_list_1.append(chain1[id]["CA"])
            else:
                s.argprint("DNA/RNA atoms added", s.options.verbose, s.options.quiet, 2)
                atom_list_1.append(chain1[id]["C4'"])

        id_list_2 = [residue.get_id() for residue in s.similarity[(resnames_chain_2, resnames_chain_1)].values()]
        atom_list_2 = []
        for id in id_list_2:
            if "CA" in chain2[id]:
                s.argprint("Residue atoms added", s.options.verbose, s.options.quiet, 2)
                atom_list_2.append(chain2[id]["CA"])
            else:
                s.argprint("DNA/RNA atoms added", s.options.verbose, s.options.quiet, 2)
                atom_list_2.append(chain2[id]["C4'"])
    else:
        atom_list_1 = []
        for residue in residue_chain_1_copy:
            atom_list_1.extend(residue.get_atoms())
        atom_list_2 = []
        for residue in residue_chain_2_copy:
            atom_list_2.extend(residue.get_atoms())

    return atom_list_1, atom_list_2
