import Bio.PDB as pdb
from Bio import pairwise2
#from Bio.SubsMat.MatrixInfo import blosum62
from Bio.pairwise2 import format_alignment

import os
import copy as cp
import string
import numpy as np

available_chain_names = set(string.ascii_letters + string.digits)

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

def Separate_Chains(pdb_file):
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

            name = "%s_chain_%s.pdb" %(interaction[0] + interaction[1], interaction[i])

            file_name = os.path.join(folder,name)

            io.save(file_name, chain())

            i += 1

    return interaction

def Get_Chains(pdb_file,pairs):
    """Divide the interaction pair and return the interaction and the chains object

    Input:
    -pdb_file = target file to process
    -pairs = pairs dictionary for see if the interaction pair is already processed
    Output:
    -interaction = string containg the two chains interacting
    -pair = dictionary:
     -Keys = Chain name (extracted from the file name)
     -Values = Chain object (corresponding to the chain name by order criteria in the file)
    """

    pdb_parser = pdb.PDBParser(PERMISSIVE=True, QUIET=True)

    pdb_structure = pdb_parser.get_structure("pdb_file", pdb_file)

    interaction = pdb_file[:-4].split("_")[-1] # Obtain str with the interaction from the file name, the order of the letters need match with the order in the pdb file(format= something_chains.pdb)

    if len(interaction) != 2 or interaction in pairs or interaction[::-1] in pairs: #if the length is not true or the pair already exists, something goes wrong
        return error

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

    Relationship dictionary:
    -Keys = str with chain_name
    -Values = set of chains with relathionship with the key one

    Pairs dictionary:
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

def Collision_Check(model_atom_list,addition_atom_list,radius):
    """Return a list of tuple containing the residue id interacting and the coordinates in conflict

    input:
    -model = list of atoms take as reference for the collitions
    -addition = list of atoms being us to check against the model
    -radius = integrer required for become the empty radious (Amstrongs) around each atom
    """

    model_ns = pdb.NeighborSearch(model_atom_list)

    collinsions_list=[]
    for atom in addition_atom_list:
        collision_list = model_ns.search(atom.get_coord(),radius,"R")
        collinsions_list.extend([tuple([str(atom.get_parent().get_id()[1]),str(x.get_id()[1]),str(atom.get_coord())]) for x in collision_list])

    return collinsions_list

def Superimpose_Chain(reference, interaction, target, pairs, radius = 7, collisions_accepted = 0):
    """Superimpose two chains with the same length and returns the related chain object moved accordingly

    Input:
    reference = chain object used as reference
    interaction = String with two chain names related between them
    target = string chain name of the common chain
    pairs = Pairs dictionary:
    -Keys = interaction pair
    -Values = dictionary:
      -Keys = Chain name
      -Values = Chain object
    """

    fixed_atoms = list(reference.get_atoms())
    mobile_atoms = list(pairs[interaction][target].get_atoms())

    sup = pdb.Superimposer()#apply superimposer tool
    sup.set_atoms(fixed_atoms, mobile_atoms)#set both lists of atoms here to get rotation
    rotran = sup.rotran

    mobile_chain_name = interaction.replace(target,"",1)

    mobile_chain = cp.deepcopy(pairs[interaction][mobile_chain_name])

    mobile_chain.transform(rotran[0],rotran[1])

    model_atom_list = []
    for chain in reference.get_parent():
        model_atom_list.extend(list(residue["CA"] for residue in chain))

    addition_atom_list = list(residue["CA"] for residue in mobile_chain)

    collisions = Collision_Check(model_atom_list, addition_atom_list, radius)

    if len(collisions) > collisions_accepted:
        return error
    else:
        return mobile_chain

def Join_Piece(model,addition):
    """Return the chain "addition" object added to the model

    input:
    -model = chain object containg the current model
    -addition = chain object moved acording a refernce chain superimposed
    -radius = integrer required for become the empty radious (Amstrongs) around each atom
    -collisions_accepted = number of atom collisions allowed in each join
    """
    character = available_chain_names.pop()#Check before run this function if the set is empty
    addition.id = character

    model.get_parent().add(addition)

    return model.get_parent()[character]

def Simulate_Model(simulation_list,relations):
    """Returns the chain with less required steps expected for build a macrocomplex with all the chains without clashes

    The chain with obtain earlier all the posible chains, whas chosen as initial chian.
    if more than one chain fullfill this condition or more than one chain are the last to fill the maximum number of chains, one of them are piked at random.
    input:
    -simulation_list = list of tuples with the following structure (initial_chain,last_added_chains,total_added_chains)
    -relations = relathionship dictionary:
     -Keys = str with chain_name
     -Values = set of chains with relathionship with the key one
    """

    possible_simulations = []
    completed_simulations = set()

    while len(completed_simulations) == 0:
        for simulation in simulation_list:
            last_added_chains = "".join(list("".join("".join(relations[chains]) for chains in simulation[1])))
            new_simulation = (simulation[0],last_added_chains,simulation[2] + last_added_chains)
            if all(chain in relations.keys() for chain in new_simulation[2]):
                completed_simulations.add(new_simulation[0])
            else:
                if new_simulation <= len(string.ascii_letters + string.digits):
                    possible_simulations.append(new_simulation)

        if len(completed_simulations) == 0:

            if len(possible_simulations) == 0:
                possible_starts = set(chain[0] for chain in simulation_list)
                return possible_starts.pop()
            else:
                return Simulate_Model(possible_simulations,relations)
        else:
            possible_starts = set(chain[0] for chain in completed_simulations)
            return possible_starts.pop()



def Chose_Start(relationships,pairs):
    """Returns a starting chain name and their respective chain object with their center of mass in the (0,0,0) coordinates

    For chose a starting point, we pick the chain/s with less relathionships making the hipotesis than will be able of make the correct macrocomplex with less steps.
    if more than one chain share the condition of have less relatihionships Simulate_Model function is applyed:

    Input:
    -relationships = Relationship dictionary:
     -Keys = str with chain_name
     -Values = set of chains with relathionship with the key one
    -pairs = Chain objects dictionary:
     -Keys = interaction pair
     -Values = dictionary:
      -Keys = Chain name
      -Values = Chain object
    """

    all_chains = relationships.keys()#Not required but put here as a note

    less_related = list(key for key,val in relationships.items() if len(val) == min((len(val) for val in relationships.values())))

    if len(less_related) == 1:

        less_related_simulation = []

        for chain in less_related: #This loop format the less related chains for the simulation
            chain_simulation = (chain,chain,chain)
            less_related_simulation.append(chain_simulation)

        starting_chain = Simulate_Model(less_related_simulation,relationships)

    else:
        starting_chain = less_related[0]

    chain_object = set()

    for chains in pairs.values():
        for object in chains.values():
            if object.get_id() == starting_chain:
                chain_object.add(object)

    starting_chain_object = chain_object.pop()

    coordinates_list = []

    for residue in starting_chain_object:
        for atom in residue:
            coordinates_list.append(atom.get_coord())

    coordinates_array = np.array(coordinates_list)
    center_of_masses = np.mean(coordinates_array,axis=0,dtype=np.float64)

    initial_model = cp.deepcopy(starting_chain_object)
    for residue in initial_model:
        for atom in residue:
            atom.coord = atom.get_vector()-center_of_masses

    character = available_chain_names.pop()#Check before run this function if the set is empty
    initial_model.id = character

    return starting_chain, initial_model

def Merge_chains(candidates, model, collisions_accepted = 0):
    """

    Input:
    -candidates = 

    """

    posible_results = {}
    for number in range(len(candidates)):
        posible_results[number] = tuple(range(len(candidates)))

    i = 0
    while i < len(candidates)-1:
        for key in range(i+1,len(candidates)):
            reference_atom_list = list(candidates[i].get_atoms())
            mobile_atom_list = list(candidates[key].get_atoms())
            collisions = Collision_Check(reference_atom_list, mobile_atom_list,radius = 7)
            if len(collisions) > collisions_accepted:
                posible_results[i] = tuple(chain for chain in posible_results[i] if chain != key)
                posible_results[key] = tuple(chain for chain in posible_results[i] if chain != i)

    results = set(posible_results.values()) #This set merge the diferent possibilities for each chain, equal models results in the same entry.

    models = []

    for combination in results:
        for chain in combination:
            models.append(Join_Piece(model,candidates[chain]))



def Build_Model(starting_chain, initial_model, relathionships, pairs):
    """

    Input:
    -starting_chain = str with the name of the model chain
    -initial_model = chain object with their center of mass at origin
    """

    relations = tuple(chain for chain in relathionships[starting_chain])
    possible_additions = []
    for chain in relations:#obtain the matching pieces

        interaction = list(pair for pair in pairs if chain in pair and starting_chain in pair)[0]


        possible_additions.append((Superimpose_Chain(initial_model, interaction, starting_chain, pairs),chain))



    if len(initial_model.get_parent().get_list()) == len(string.ascii_letters + string.digits):
        return model






if __name__ == "__main__":

    A,B=Parse_List(["../example1/pair_his3_sc_XA.pdb","../example1/pair_his3_sc_XB.pdb"])

    print(B["XA"]["A"][1]["CA"].get_coord())

    print(Superimpose_Chain(B,"XA","X",B["XB"]["X"])[1]["CA"].get_coord())

    # print(Collision_Check(B["XA"]["X"],B["XB"]["X"],8))

    C = Join_Piece(B["XB"]["X"],B["XB"]["B"])

    for chain in C.get_parent():
        print (chain)

    print(Chose_Start(A,B))


# def get_alignment(file_1, file_2):
#     pdb_parser = pdb.PDBParser(PERMISSIVE=True, QUIET=True)
#     structure_1 = pdb_parser.get_structure("file_1", file_1)
#     structure_2 = pdb_parser.get_structure("file_2", file_2)


# #################################
# ##Extract Sequence from structure
# #################################
# ##first, you must extract the polypeptides -- using alpha carbons here
# ## ppb = PPBuilder() for c-n distances
# ##second, we will extract the sequence
#     seq = []
#     ppb = pdb.CaPPBuilder()
#     for polypeptide in ppb.build_peptides(structure_1):
#         sequence_ref = polypeptide.get_sequence()
#         seq.append(sequence_ref)
#     for polypeptide in ppb.build_peptides(structure_2):
#         sequence_sample = polypeptide.get_sequence()
#         seq.append(sequence_sample)
#     return seq

# #################################
# ##Align sequences
# #################################
# ##now we have sequences from two structures
# ##next, let's globally align them
#     align = pairwise2.align.globalxx(sequence_ref, sequence_sample)
# ##format_alignment to get pretty print
#     #print(format_alignment(*align[0]))

