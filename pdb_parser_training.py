import Bio.PDB as pdb
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

import os
import copy as cp
import string
import numpy as np
import itertools
import pickle
from project import *
import settings


def Separate_Chains(pdb_file): #destroy once you don't need nothing from it
    """Separate the two chains and return their name in a list

    Input:
    -pdb file = target file
    Output:
    -interaction = list with chain information
    """
    folder = "pdb_chains"

    if not Check_folder(folder):

        return False

    pdb_parser = pdb.PDBParser(PERMISSIVE=True, QUIET=True)

    pdb_structure = pdb_parser.get_structure("pdb_file", pdb_file)

    interaction = list(pdb_file[:-4].split("_")[-1]) # Obtain 2 length lists with the chain names from file name, the order of the letters need match with the order in the pdb file(format= something_chains.pdb)

    if len(interaction) != 2: #if the length is not true, something goes wrong

        print(settings.IncorrectName(interaction))

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

    return interaction



def Delete_Folder(folder): #support if is required finally
    """Delete all the files in a folder and the folder
    Input:
    folder = contains pdb files
    """

    for generated_file in os.listdir(folder):
        file_path = os.path.join(folder, generated_file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
        except Exception as e:
            print(e) #The error is printed for be able of know the reason of the problem
            return False

    try:
        os.rmdir(folder)
    except Exception as e:
        print(e)
        return False

    return True

if __name__ == "__main__":

    print("Program start")
    settings.init()

    unable_pickle = settings.options.pickle
    folder_arg = settings.options.folder
    verbose = settings.options.verbose
    quiet = settings.options.quiet
    initial_chain_arg = settings.options.initial
    collisions_accepted_arg = settings.options.collisions
    radius_arg = settings.options.radius
    percent_arg = settings.options.similarity
    pdb_list_arg = settings.options.pdb_list
    output_arg = settings.options.output

    if unable_pickle:
        relation_pickle = os.path.join(folder_arg, "relationships.p")
        pairs_pickle = os.path.join(folder_arg, "pairs.p")
        try:
            relationships = pickle.load(open(relation_pickle,"rb"))
            pairs = pickle.load(open(pairs_pickle,"rb"))
        except:
            settings.argprint("No pickle available! Generating a new one...", verbose, quiet, 2)
            onlyfiles = settings.Get_PdbList(folder_arg, pdb_list_arg)
            relationships, pairs = Parse_List(onlyfiles)
            out_fd = open(relation_pickle,"wb")
            pickle.dump(relationships ,out_fd)
            out_fd.close()
            out_fd = open(pairs_pickle,"wb")
            pickle.dump(pairs,out_fd)
            out_fd.close()
    else:
        onlyfiles = settings.Get_PdbList(folder_arg, pdb_list_arg)
        relationships, pairs = Parse_List(onlyfiles)


    settings.argprint("files_parsed", verbose, quiet, 0)
    print("Number of chains detected: %s" %(len(relationships.keys())))

    initial_model = Chose_Start(relationships, pairs, available_chain_names, initial_chain_arg)
    models = Build_model(initial_model, relationships, pairs, collisions_accepted_arg, radius_arg, percent_arg)

    print("Writing pdb file...")

    print(models)


    io = pdb.PDBIO()
    i = 1
    for model in models:
        filename = "%s_model_%s.pdb" %(os.path.basename(folder_arg),i)
        io.set_structure(model[0])
        io.save(os.path.join(output_arg, filename))
        i += 1
