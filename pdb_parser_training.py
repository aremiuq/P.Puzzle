import Bio.PDB as pdb
from Bio import pairwise2
#from Bio.SubsMat.MatrixInfo import blosum62
from Bio.pairwise2 import format_alignment

import os

def Check_folder(folder):
    """If folder don't exists creates it in the actual path"""
    path = os.getcwd()

    new_path = os.path.join(path,folder)

    if not os.path.isdir(new_path):

        try:
            os.mkdir(new_path)
        except OSError:
            print ("Creation of the directory %s failed" % path)
            return False
        else:
            print ("Successfully created the directory %s " % path)

    return True

def Get_Chains(pdb_file,pairs):
    """Separate the two chains and return their name in a list

    pdb_file = target file to process
    pairs = pairs dictionary for see if the interaction pair is already processed
    """

    pdb_parser = pdb.PDBParser(PERMISSIVE=True, QUIET=True)

    pdb_structure = pdb_parser.get_structure("pdb_file", pdb_file)

    interaction = pdb_file[:-4].split("_")[-1] # Obtain str with the interaction from the file name, the order of the letters need match with the order in the pdb file(format= something_chains.pdb)

    if len(interaction) != 2 or interaction in pairs or interaction[::-1] in pairs: #if the length is not true or the pair already exists, something goes wrong
        return error

    pair = {}

    for model in pdb_structure:
        for chain in model:

            pair[interaction[0]] = chain

            interaction = interaction[::-1]

    return interaction, pair


def Parse_List(list):
    """Apply the Get_Chain to all the list and returns a relationship and pairs dictionary

    Relationship dictionary:
    -Keys = str with chain_name
    -Values = set of chains with relathionship with the key one

    Pairs dictionary:
    -Keys = interaction pair
    -Values = dictionary:
      -Keys = Chain name
      -Values = Chain object
    """

    if not "pairs" in locals():
        pairs = {}

    RelationDict = {}

    for pdb_file in list:

        chains,pair = Get_Chains(pdb_file, pairs)

        pairs[chains] = pair

        for letter in chains:

            if letter in RelationDict:
                RelationDict[letter].add(chains[1])
            else:
                RelationDict[letter] = set([chains[1]])

            chains = chains[::-1] #reverse the order of the chains

    return RelationDict, pairs

def Delete_Folder(folder):
    """Delete all the files in a folder and the folder"""

    for generated_file in os.listdir(folder):
        file_path = os.path.join(folder, generated_file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
        except Exception as e:
            print(e) #Error
            return False

    try:
        os.rmdir(folder)
    except Exception as e:
        print(e) #Error
        return False

    return True



def get_alignment(file_1, file_2):
    pdb_parser = pdb.PDBParser(PERMISSIVE=True, QUIET=True)
    structure_1 = pdb_parser.get_structure("file_1", file_1)
    structure_2 = pdb_parser.get_structure("file_2", file_2)


#################################
##Extract Sequence from structure
#################################
##first, you must extract the polypeptides -- using alpha carbons here
## ppb = PPBuilder() for c-n distances
##second, we will extract the sequence
    seq = []
    ppb = pdb.CaPPBuilder()
    for polypeptide in ppb.build_peptides(structure_1):
        sequence_ref = polypeptide.get_sequence()
        seq.append(sequence_ref)
    for polypeptide in ppb.build_peptides(structure_2):
        sequence_sample = polypeptide.get_sequence()
        seq.append(sequence_sample)
    return seq
    
#################################
##Align sequences
#################################
##now we have sequences from two structures
##next, let's globally align them 
    align = pairwise2.align.globalxx(sequence_ref, sequence_sample)
##format_alignment to get pretty print 
    #print(format_alignment(*align[0]))

 
def superimpose(file_1, file_2):
#################################
##Superimpose objects
#################################
    pdb_parser = pdb.PDBParser(PERMISSIVE=True, QUIET=True)
    structure_1 = pdb_parser.get_structure("file_1", file_1)
    structure_2 = pdb_parser.get_structure("file_2", file_2)
##retrieve list of atoms here for each structure 
    structure1_atoms = list(structure_1.get_atoms())
    structure2_atoms = list(structure_2.get_atoms())
##apply superimposer tool 
    sup = pdb.Superimposer()


##set both lists of atoms here to get rmsd  
    sup.set_atoms(structure1_atoms, structure2_atoms)
    #return "RMSD: " + str(sup.rms)
    return str(sup.rotran)





A,B=Parse_List(["../example1/pair_his3_sc_XA.pdb","../example1/pair_his3_sc_XB.pdb"])

print(A)
print(B)

# def Separate_Chains(pdb_file):
#     """Separate the two chains and return their name in a list"""

#     folder = "pdb_chains"

#     if not Check_folder(folder):

#         return False

#     pdb_parser = pdb.PDBParser(PERMISSIVE=True, QUIET=True)

#     pdb_structure = pdb_parser.get_structure("pdb_file", pdb_file)

#     interaction = list(pdb_file[:-4].split("_")[-1]) # Obtain 2 len list with the chain names from file name, the order of the letters need match with the order in the pdb file(format= something_chains.pdb)

#     if len(interaction) != 2: #if the length is not true, something goes wrong

#         return error

#     i=0

#     for model in pdb_structure:
#         for chain in model:

#             id=chain.get_id()

#             class chain(pdb.Select):
#                 def accept_chain(self, chain):
#                     if chain.get_id()== id:
#                         return True
#                     else:
#                         return False

#             io = pdb.PDBIO()

#             io.set_structure(pdb_structure)

#             name = "%s_chain_%s.pdb" %(interaction[0] + interaction[1], interaction[i])

#             file_name = os.path.join(folder,name)

#             io.save(file_name, chain())

#             i += 1

#     return interaction
