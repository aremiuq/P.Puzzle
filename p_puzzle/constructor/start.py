
#start.py imports

import Bio.PDB as pdb
import numpy as np
import copy as cp
import string

#start.py functions

def Simulate_Model(simulation_list,relations):
    """Returns a string with the chain that has less required steps expected to build a macrocomplex with all the chains without clashes.
    The chain obtain earlier all the posible chains, whas chosen as initial chian.
    If more than one chain fulfills this condition or more than one chain are the last to fill the maximum number of chains, one of them are piked at random.
    Input:
    -simulation_list = list of tuples with the following structure (initial_chain,last_added_chains,total_added_chains)
    -relations = relathionship dictionary:
        -Keys = str with chain_name
        -Values = set of chains with relathionship with the key one
    Output:
    -chain_chosed = return the name of one chain at random from the set
    """

    possible_simulations = []
    completed_simulations = set()

    for simulation in simulation_list:
        last_added_chains = "".join(list("".join("".join(relations[chains]) for chains in simulation[1])))
        new_simulation = (simulation[0],last_added_chains,simulation[2] + last_added_chains)
        if all(chain in new_simulation[2] for chain in relations.keys()):
            completed_simulations.add(new_simulation[0])
        else:
            if len(new_simulation[2]) <= len(string.ascii_letters + string.digits):
                possible_simulations.append(new_simulation)

    if len(completed_simulations) == 0:

        if len(possible_simulations) == 0:
            possible_starts = set(chain[0] for chain in simulation_list)
            return possible_starts.pop()
        else:
            return Simulate_Model(possible_simulations,relations)
    else:
        possible_starts = set(chain[0] for chain in completed_simulations)
        chain_chosed = possible_starts.pop()
        return chain_chosed

def Chose_Start(relationships, pairs, available_chain_names, selected_chain = None):
    """Returns a starting chain name and their respective chain object with their center of mass in the (0,0,0) coordinates
    To choose starting point, the program picks the chain/s with less relathionships making the hipotesis than will be able of make the correct macrocomplex with less steps.
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
    -available_chain_names = set with the avaiable names for the chains in the pdb format
    -selected_chain = optional, instead of use the simulation for guess the shortest past, introduce the name of the starting chain.
    Output:
    - tuple containing starting chain name, starting chain object, and the list with the remaining aviable names
    """

    if selected_chain is None:
        all_chains = relationships.keys()#Not required but put here as a note

        less_related = list(key for key,val in relationships.items() if len(val) == min((len(val) for val in relationships.values())))

        if len(less_related) != 1:

            less_related_simulation = []

            for chain in less_related: #This loop format the less related chains for the simulation
                chain_simulation = (chain,chain,chain)
                less_related_simulation.append(chain_simulation)

            starting_chain = Simulate_Model(less_related_simulation,relationships)

        else:
            starting_chain = less_related[0]
    else:
        starting_chain = selected_chain

    chain_object = set()

    for chains in pairs.values():
        for key,value in chains.items():
            if key == starting_chain:
                chain_object.add(value)

    starting_chain_object = chain_object.pop()
    initial_chian = cp.deepcopy(starting_chain_object)
    character = available_chain_names.pop()
    initial_chian.id = character
    initial_model = initial_chian.get_parent()
    initial_model.get_parent().id = "my_model"

    coordinates_list = []

    for residue in starting_chain_object:
        for atom in residue:
            coordinates_list.append(atom.get_coord())

    coordinates_array = np.array(coordinates_list)
    center_of_masses = np.mean(coordinates_array,axis=0,dtype=np.float64)

    for atom in initial_model.get_atoms():#center the chain in the point (0,0,0)
        new_coords = atom.get_vector()-center_of_masses
        atom.coord = new_coords.get_array()

    return tuple((initial_model,[(starting_chain, character)], available_chain_names))
