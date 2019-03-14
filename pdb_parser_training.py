import Bio.PDB as pdb

import os

def Check_folder(folder):
    """If folder don't exists creates it in the actual path"""
    path = os.getcwd()

    new_path = path + "/" + folder

    if not os.path.isdir(new_path):

        try:
            os.mkdir(new_path)
        except OSError:
            print ("Creation of the directory %s failed" % path)
            return False
        else:
            print ("Successfully created the directory %s " % path)

    return True

def Get_Chains(pdb_file):
    """Separate the two chains and return their name in a list"""

    folder = "pdb_chains"

    if not Check_folder(folder):

        return False

    pdb_parser = pdb.PDBParser(PERMISSIVE=True, QUIET=True)

    pdb_structure = pdb_parser.get_structure("pdb_file", pdb_file)

    interaction = list(pdb_file[:-4].split("_")[-1]) # Obtain 2 len list with the chain names from file name, the order of the letters need match with the order in the pdb file(format= something_chains.pdb)

    if len(interaction) != 2: #if the length is not true, something goes wrong

        return error

    i=0

    for model in pdb_structure:
        for chain in model:

            id=chain.get_id()

            class chain(pdb.Select):
                def accept_chain(self, chain):
                    if chain.get_id()== id:
                        return True
                    else:
                        return False

            io = pdb.PDBIO()

            io.set_structure(pdb_structure)

            io.save("./%s/chain_%s.pdb" %(folder,interaction[i]), chain())

            i += 1

    return interaction

def Parse_List(list):
    """Apply the Get_Chain to all the list and returns a relationship dictionary

    Relationship dictionary:
    -Keys = str with chain_name
    -Values = set of chains with relathionship with the key one
    """

    RelationDict = {}

    for pdb_file in list:

        chains = Get_Chains(pdb_file)

        iterative_chains = chains.copy()

        for letter in iterative_chains:

            if letter in RelationDict:
                RelationDict[letter].add(chains[1])
            else:
                RelationDict[letter] = set([chains[1]])

            chains = chains[::-1] #reverse the order of the chains

    return RelationDict





A=Parse_List(["./example1/pair_his3_sc_XA.pdb","./example1/pair_his3_sc_XB.pdb"])

print(A)
