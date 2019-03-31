
#pdb_files.py imports

from ..support import settings as s
import Bio.PDB as pdb
import copy as cp

#pdb_files.py functions

def Get_Chains(pdb_file,pairs):
    """Divide the interaction pair and return the interaction and the chains object
    Input:
    -pdb_file = target file to process
    -pairs = pairs dictionary to check if the interaction pair is already processed
    Output:
    -interaction = string containing the two chains interacting
    -pair = dictionary:
     -Keys = Chain name (extracted from the file name)
     -Values = Chain object (corresponding to the chain name by order criteria in the file)
    """
    pdb_parser = pdb.PDBParser(PERMISSIVE=s.options.strict, QUIET=True)

    interaction = pdb_file[:-4].split("_")[-1] # Obtain str with the interaction from the file name, the order of the letters need match with the order in the pdb file(format= something_chains.pdb)
    # interaction = pdb_file[:-4].split("/")[-1]
    pdb_structure = pdb_parser.get_structure(interaction, pdb_file)

    s.argprint(interaction, s.options.verbose, s.options.quiet, 1)
    if len(interaction) != 2:
        s.argprint(s.IncorrectName(interaction), s.options.verbose, s.options.quiet, 0)
        exit()
    if interaction in pairs:
        s.argprint(s.RepeatedChain(interaction), s.options.verbose, s.options.quiet, 1)

    pair = {}

    for model in pdb_structure: #In this loop a pdb object is created with only the information of 1 chain
        for chain in model:
            model_copy = cp.deepcopy(model)
            model_copy.detach_child(chain.get_id())
            for chain_copy in model_copy:
                pair[interaction[1]] = chain_copy
            interaction = interaction[::-1]

    del pdb_structure

    return interaction, pair

def Parse_List(list):
    """Apply the Get_Chain to all the list and returns a relationship and pairs dictionary
    Input:
    -list = pdb information
    Output:
    -Relationship dictionary:
        -Keys = str with chain_name
        -Values = set of chains with relathionship with the key one
    -Pairs dictionary:
       -Keys = interaction pair
        -Values = dictionary:
              -Keys = Chain name
              -Values = Chain object
    """
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
