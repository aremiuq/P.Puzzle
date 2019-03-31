
#core.py imports

from ..support import settings as s
from .checker import *
import Bio.PDB as pdb
import copy as cp
import itertools

#core.py functions

def Collision_Check(model_atom_list,addition_atom_list, radius):
    """Return a list of tuple containing the residue id interacting and the coordinates in conflict
    Input:
    -model = list of atoms take as reference for the collitions
    -addition = list of atoms being us to check against the model
    -radius = integrer required for become the empty radious (Amstrongs) around each atom
    Output:
    -collisions_list = list of tuples of coordinates where atoms collide
    """

    model_ns = pdb.NeighborSearch(model_atom_list)

    collisions_list=[]
    for atom in addition_atom_list:
        collision_list = model_ns.search(atom.get_coord(),radius,"R")
        collisions_list.extend([tuple([str(atom.get_parent().get_id()[1]),str(x.get_id()[1]),str(atom.get_coord())]) for x in collision_list])#The list is complemented with information about the clash

    return collisions_list

def Superimpose_Chain(reference, interaction, target, pairs, collisions_accepted, radius, percent):
    """Superimpose two chains with the same length and returns the related chain object moved accordingly
    Input:
    -reference = chain object used as reference
    -interaction = String with two chain names related between them
    -target = string chain name of the common chain
    -pairs = Pairs dictionary:
        -Keys = interaction pair
        -Values = dictionary:
            -Keys = Chain name
            -Values = Chain object
    -radius = integrer required for become the empty radius (measure in Amstrongs) around each atom, user input
    -collisions_accepted = number of atom collisions allowed in each join, user input
    -percent = similarity percentage accepted, user input
    Output:
    -mobile_chain = list of atoms from the mobile chain that have had their coordinates moved for superimposing return None if the chain can't be added to the actual model
    """

    fixed_atoms, mobile_atoms = Check_Similarity(reference, pairs[interaction][target], percent)

    sup = pdb.Superimposer()#apply superimposer tool
    sup.set_atoms(fixed_atoms, mobile_atoms)#set both lists of atoms here to get rotation
    rotran = sup.rotran

    mobile_chain_name = interaction.replace(target,"",1)# Obtain the chain name of the mobile one

    mobile_chain = cp.deepcopy(pairs[interaction][mobile_chain_name])

    mobile_chain.transform(rotran[0],rotran[1])

    model_atom_list = []
    for chain in reference.get_parent():
        model_atom_list.extend(chain.get_atoms())

    addition_atom_list = list(mobile_chain.get_atoms())

    collisions = Collision_Check(model_atom_list, addition_atom_list, radius)

    if len(collisions) > 0 :
        s.argprint("Number of collisions:", s.options.verbose, s.options.quiet, 2)
        s.argprint(len(collisions), s.options.verbose, s.options.quiet, 2)

    if len(collisions) > collisions_accepted:
        e = s.CollisionAppears(target, mobile_chain_name, collisions)
        s.argprint(e, s.options.verbose, s.options.quiet, 2)
        s.argprint(e.Get_Collisions(), s.options.verbose, s.options.quiet, 3)
        return None
    else:
        return mobile_chain

def Join_Piece(model,addition,available_chain_names):
    """Return the chain "addition" object added to the model
    Input:
    -model = model object containing the current model
    -addition = chain object moved acording a refernce chain superimposed
    Output:
    -model = model object containg the original model and addition
    -character = id of model being added
    -available_chain_names = list of chain ids that are still available
    """
    character = available_chain_names.pop()
    s.argprint("Model remaining available names:", s.options.verbose, s.options.quiet, 3)
    s.argprint(len(available_chain_names), s.options.verbose, s.options.quiet, 3)
    addition.id = character

    model.add(addition)

    return model, character, available_chain_names

def Merge_chains(candidates, model_object, available_chain_names, collisions_accepted, radius):
    """Returns all the possible combinations as a list of models, of the model and the candidates
    Input:
    -candidates = list of tuples with: (chain_object,chain_name) structure
    -model_object = model object, for join the candidates to them
    -available_chain_names = set with the avaiable names for the chains in the pdb format
    -collisions_accepted = number of atom collisions allowed in each join, user input
    -radius = integrer required for become the empty radius (measure in Amstrongs) around each atom, user input
    Output:
    -resulting_models = list of tuples with (model_object, [(chain_name,chain_id), ...], set(avaiable_names))
    """

    s.argprint("Merging Chains...", s.options.verbose, s.options.quiet, 1)
    s.argprint("Candidates for merging:", s.options.verbose, s.options.quiet, 2)
    s.argprint(len(candidates), s.options.verbose, s.options.quiet, 2)
    s.argprint(candidates, s.options.verbose, s.options.quiet, 3)

    possible_results = set() #This become a dictionary containg the templates of the resulting models
    for number in range(len(candidates)):
        possible_results.add((number,tuple(range(len(candidates)))))#format example in the case of 4 candidates {(0,(0123))}

    for first,second in itertools.combinations(tuple(range(len(candidates))),2):

        reference_atom_list = list(candidates[first][0].get_atoms())
        mobile_atom_list = list(candidates[second][0].get_atoms())
        collisions = Collision_Check(reference_atom_list, mobile_atom_list, radius)
        if len(collisions) > collisions_accepted:#extract the collition pair from the posible results
            s.argprint("Number of collisions:", s.options.verbose, s.options.quiet, 2)
            s.argprint(len(collisions), s.options.verbose, s.options.quiet, 2)

            if len(collisions) >= len(reference_atom_list):#Remove the candidate if is equal
                repeated_removed = set()
                for model in possible_results:
                    if model[0] != second:
                        repeated_removed.add((model[0],tuple(chain for chain in model[1] if chain != second)))
                        possible_results = repeated_removed

            new_possible_results = set()
            for model in possible_results:
                if model[0] == first:
                    new_possible_results.add((first,tuple(chain for chain in model[1] if chain != second)))
                elif model[1] == second:
                    new_possible_results.add((second,tuple(chain for chain in model[1] if chain != first)))
                else:
                    new_possible_results.add((model[0],tuple(chain for chain in model[1] if chain != first)))
                    new_possible_results.add((model[0],tuple(chain for chain in model[1] if chain != second)))
            possible_results = new_possible_results

    possible_results_unrepeated = set([tuple(element[1]) for element in possible_results]) #remove repeated combinations
    results_unprocesed = [set(element) for element in possible_results_unrepeated] #format for remove combinations with more possibilities
    results = cp.deepcopy(results_unprocesed)
    for result_1, result_2 in itertools.permutations(results_unprocesed,2):
        if result_1.issubset(result_2):
            if result_1 in results:
                results.remove(result_1)

    s.argprint("Possible models were found! Processing now...", s.options.verbose, s.options.quiet, 2)
    s.argprint(len(results), s.options.verbose, s.options.quiet, 2)
    s.argprint(results, s.options.verbose, s.options.quiet, 3)

    resulting_models = []
    for combination in results:
        s.argprint("Combination %s"%(combination), s.options.verbose, s.options.quiet, 3)
        if len(available_chain_names) < len(combination):
            s.argprint("Available chain names have been exhausted. The following chains can't be added:", s.options.verbose, s.options.quiet, 2)
            for number in combination:
                s.argprint(candidates[number], s.options.verbose, s.options.quiet, 2)
            resulting_models.append((model_object,None,set()))
        else:
            added_chains = []
            merged_model = cp.deepcopy(model_object)
            merged_avaiable_chain_names = cp.deepcopy(available_chain_names)
            for number in combination:
                s.argprint("%s started"%(number), s.options.verbose, s.options.quiet, 3)
                candidate_copied = cp.deepcopy(candidates[number][0])
                merged_model, chain_name, merged_available_chain_names = Join_Piece(merged_model, candidate_copied, merged_avaiable_chain_names)
                added_chains.append(((candidates[number][1]), chain_name))#tupple containing (chain_name,random_name_gived)

            resulting_models.append((merged_model, added_chains, merged_available_chain_names))

    s.argprint("Resulting models:", s.options.verbose, s.options.quiet, 2)
    s.argprint(len(resulting_models), s.options.verbose, s.options.quiet, 2)
    s.argprint(resulting_models, s.options.verbose, s.options.quiet, 3)

    return resulting_models #This output contains a list of tuples with (model_object,list of tuples with : (chain_name,chain_id) of the last added chains, set(avaiable_names))

def Check_Chains(model_tupled, relationships, pairs, collisions_accepted, radius, percent):
    """Return the matching chains of one model
    Input:
    -model_tupled = tuple with (model_object, list of tuples with : (chain_name,chain_id), set(avaiable_names))
    -relationships = Relationship dictionary:
        -Keys = str with chain_name
        -Values = set of chains with relathionship with the key one
    -pairs = Chain objects dictionary:
        -Keys = interaction pair
        -Values = dictionary:
            -Keys = Chain name
            -Values = Chain object
    -collisions_accepted = number of atom collisions allowed in each join, user input
    -radius = integrer required for become the empty radius (measure in Amstrongs) around each atom, user input
    -percent = minimum similarity percentage accepted, user input

    Output:
    -possible_additions: list of tuples with the following format: (Chain_object,chain_name)
    """

    s.argprint("Checking chains...", s.options.verbose, s.options.quiet, 1)
    possible_additions = []
    for element in model_tupled[1]:
        relations = tuple(chain for chain in relationships[element[0]])
        s.argprint("Relations obtained:", s.options.verbose, s.options.quiet, 2)
        s.argprint(len(relations), s.options.verbose, s.options.quiet, 2)
        s.argprint(relations, s.options.verbose, s.options.quiet, 3)
        for chain in relations:#obtain the matching pieces
            interaction = element[0] + chain
            if not interaction in pairs:#if both orders dosn't exist
                interaction = list(pair for pair in pairs if chain in pair and element[0] in pair)[0]
            addition_tested = (Superimpose_Chain(model_tupled[0][element[1]], interaction, element[0], pairs, collisions_accepted, radius, percent),chain)
            s.argprint("Superimpose results:", s.options.verbose, s.options.quiet, 3)
            s.argprint(addition_tested, s.options.verbose, s.options.quiet, 3)
            if addition_tested[0] != None:
                possible_additions.append(addition_tested)

    if len(possible_additions) == 0:
        s.argprint("Any chain can be added", s.options.verbose, s.options.quiet, 1)
        return None
    else:
        s.argprint("Resulting chains:", s.options.verbose, s.options.quiet, 2)
        s.argprint(len(possible_additions), s.options.verbose, s.options.quiet, 2)
        s.argprint(possible_additions, s.options.verbose, s.options.quiet, 3)
        return possible_additions

def Build_model(model_tupled, relationships, pairs, collisions_accepted = 30, radius = 2, percent = 95):
    """Returns a list of tuples of all the possible models taking in account the given restrictions.
    Input:
    -model_tupled = Tuple with (model_object, (chain_name,chain_id), set(avaiable_names))
    -relationships = Relationship dictionary:
        -Keys = str with chain_name
        -Values = set of chains with relathionship with the key one
    -pairs = Chain objects dictionary:
        -Keys = interaction pair
        -Values = dictionary:
            -Keys = Chain name
            -Values = Chain object
    -collisions_accepted = number of atom collisions allowed in each join, user input
    -radius = integrer required for become the empty radius (measure in Amstrongs) around each atom, user input
    -percent = minimum similarity percentatge accepted between chains with same name

    Output:
    -finished_models = list of completed models
    """

    s.argprint("Building model...", s.options.verbose, s.options.quiet, 1)

    finished_models = []
    pieces_for_merge = Check_Chains(model_tupled, relationships, pairs, collisions_accepted, radius, percent)
    if pieces_for_merge == None:
        finished_models.append(model_tupled)
    else:
        resulting_models = Merge_chains(pieces_for_merge, model_tupled[0], model_tupled[2], collisions_accepted, radius)
        for resulting_model in resulting_models:
            if len(resulting_model[2]) == 0:
                finished_models.extend(resulting_model)
            else:
                finished_models.extend(Build_model(resulting_model, relationships, pairs, collisions_accepted, radius))

    s.argprint("Complete models:", s.options.verbose, s.options.quiet, 2)
    s.argprint(len(finished_models), s.options.verbose, s.options.quiet, 2)
    s.argprint(finished_models, s.options.verbose, s.options.quiet, 3)
    return finished_models
