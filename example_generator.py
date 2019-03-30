import Bio.PDB as pdb
from Bio import pairwise2
#from Bio.SubsMat.MatrixInfo import blosum62
from Bio.pairwise2 import format_alignment

import os
import copy as cp
import string
import numpy as np
import itertools
import pickle
from project import Get_Pairwise
from pdb_parser_training import *

def Create_Example(pdb_file,folder):
    """Parse pdb files into files containing individual pairwise interactions
    Input:
    -pdb_file = target pdb file 
    -folder = target folder housing pdb files for parsing 
    Output:
    -parsed files containing individual pairwise interactions 
    """

    pdb_parser = pdb.PDBParser(PERMISSIVE=True, QUIET=True)
    file_name = os.path.basename(pdb_file)
    print(file_name)
    pdb_structure = pdb_parser.get_structure(file_name, pdb_file)
    print("File(s) have been parsed.")
    chain_names = []
    Check_folder(folder)
    for model in pdb_structure:
        empty_model = cp.deepcopy(pdb_structure[model.get_id()])
        for chain in model:
            id_name = chain.get_id()
            chain_names.append(id_name)
            empty_model.detach_child(id_name)

        for number in range(len(chain_names)):
            if number % 2 == 1 and number != 0:
                interaction = chain_names[number - 1] + chain_names[number]
                new_model = cp.deepcopy(empty_model)
                new_model.add(model[chain_names[number - 1]])
                new_model.add(model[chain_names[number]])
                io = pdb.PDBIO()
                io.set_structure(new_model)
                io.save(os.path.join(folder, "%s_%s.pdb" %(file_name, interaction)))
                try:
                    interaction = chain_names[number] + chain_names[number + 1]
                    new_model = cp.deepcopy(empty_model)
                    new_model.add(model[chain_names[number]])
                    new_model.add(model[chain_names[number + 1]])
                    io = pdb.PDBIO()
                    io.set_structure(new_model)
                    io.save(os.path.join(folder, "%s_%s.pdb" %(file_name, interaction)))
                except:
                    print("Parsing complete!")
                    pass

def Create_Example2(pdb_file, folder, percent = 95):

    equivalence = {}

    pdb_parser = pdb.PDBParser(PERMISSIVE=True, QUIET=True)
    file_name = os.path.basename(pdb_file)
    print(file_name)
    pdb_structure = pdb_parser.get_structure(file_name, pdb_file)
    print("File(s) have been parsed")
    chain_names = []
    Check_folder(folder)
    for model in pdb_structure:
        empty_model = cp.deepcopy(pdb_structure[model.get_id()])
        for chain in model:
            id_name = chain.get_id()
            equivalence[id_name] = id_name
            chain_names.append(id_name)
            empty_model.detach_child(id_name)

        for chain_1, chain_2 in itertools.combinations(chain_names,2):
            residue_chain_1 = model[chain_1].get_residues()
            residue_chain_2 = model[chain_2].get_residues()
            resnames_chain_1 = []
            for residue in residue_chain_1:
                resnames_chain_1.append(residue.get_resname())
            resnames_chain_1 = tuple(resnames_chain_1)
            resnames_chain_2 = []
            for residue in residue_chain_2:
                resnames_chain_2.append(residue.get_resname())
            resnames_chain_2 = tuple(resnames_chain_2)
            if resnames_chain_1 == resnames_chain_2:
                equivalence[chain_2] = equivalence[chain_1]
            elif set(resnames_chain_1).issubset(set(resnames_chain_2)) or set(resnames_chain_2).issubset(set(resnames_chain_1)):
                pass
            else:
                align = Get_Pairwise(model[chain_1], model[chain_2])
                try:
                    sim_percent = (align[2]/align[4]) * 100
                    print(sim_percent)
                except TypeError as e:
                    print("Comparing aa with nucleic acids")
                    sim_percent = percent
                if sim_percent <= percent:
                    pass
                else:
                    equivalence[chain_2] = equivalence[chain_1]

        for number in range(len(chain_names)):
            if number % 2 == 1 and number != 0:
                interaction = chain_names[number - 1] + chain_names[number]
                new_model = cp.deepcopy(empty_model)
                new_model.add(model[chain_names[number - 1]])
                new_model.add(model[chain_names[number]])
                io = pdb.PDBIO()
                io.set_structure(new_model)
                io.save(os.path.join(folder, "%s_%s%s.pdb" %(file_name, equivalence[interaction[0]], equivalence[interaction[1]])))
                try:
                    interaction = chain_names[number] + chain_names[number + 1]
                    new_model = cp.deepcopy(empty_model)
                    new_model.add(model[chain_names[number]])
                    new_model.add(model[chain_names[number + 1]])
                    io = pdb.PDBIO()
                    io.set_structure(new_model)
                    io.save(os.path.join(folder, "%s_%s%s.pdb" %(file_name, equivalence[interaction[0]], equivalence[interaction[1]])))
                except:
                    print("error")
                    pass

Create_Example("../example_generator/nucleosome_3kuy.pdb","../example_generator/nucleosome_3kuy_gen2")

