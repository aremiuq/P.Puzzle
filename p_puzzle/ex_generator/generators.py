
#generators.py imports

from p_puzzle.constructor.checker import *
import p_puzzle.support.settings as s
import Bio.PDB as pdb
import copy as cp
import os

#generators.py functions

def Create_Example(pdb_file,folder):
    """Parse pdb files into files containing individual pairwise interactions
    Input:
    -pdb_file = target pdb file
    -folder = target folder housing the output
    Output:
    -folder containing the parsed files containing individual pairwise interactions
    """

    pdb_parser = pdb.PDBParser(PERMISSIVE=s.options.strict, QUIET=True)
    file_name = os.path.basename(pdb_file)
    s.Check_folder(folder)
    output_folder = os.path.join(folder, file_name[:-4])
    s.Check_folder(output_folder)
    pdb_structure = pdb_parser.get_structure(file_name, pdb_file)
    s.argprint("File have been parsed.", s.options.verbose, s.options.quiet, 1)
    chain_names = []
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
                io.save(os.path.join(output_folder, "%s_%s.pdb" %(file_name, interaction)))
                s.argprint("%s_%s.pdb created" %(file_name, interaction), s.options.verbose, s.options.quiet, 2)
                try:
                    interaction = chain_names[number] + chain_names[number + 1]
                    new_model = cp.deepcopy(empty_model)
                    new_model.add(model[chain_names[number]])
                    new_model.add(model[chain_names[number + 1]])
                    io = pdb.PDBIO()
                    io.set_structure(new_model)
                    io.save(os.path.join(output_folder, "%s_%s.pdb" %(file_name, interaction)))
                    s.argprint("%s_%s.pdb created" %(file_name, interaction), s.options.verbose, s.options.quiet, 2)
                except:
                    s.argprint("Example complete!", s.options.verbose, s.options.quiet, 0)

def Create_Example2(pdb_file, folder, percent = 95):
    """Parse pdb files into files containing individual pairwise interactions with the similar chains having the same name

    Be carrefoul, the program is nop optimized, big files can take so much time.
    Input:
    -pdb_file = target pdb file
    -folder = target folder housing the output
    Output:
    -folder containing the parsed files containing individual pairwise interactions
    """

    equivalence = {}

    pdb_parser = pdb.PDBParser(PERMISSIVE=s.options.strict, QUIET=True)
    file_name = os.path.basename(pdb_file)
    s.Check_folder(folder)
    output_folder = os.path.join(folder, file_name[:-4])
    s.Check_folder(output_folder)
    pdb_structure = pdb_parser.get_structure(file_name, pdb_file)
    s.argprint("File have been parsed.", s.options.verbose, s.options.quiet, 1)
    chain_names = []
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
                    s.argprint("Similarity between chains:", s.options.verbose, s.options.quiet, 3)
                    s.argprint(sim_percent, s.options.verbose, s.options.quiet, 3)
                except TypeError as e:
                    s.argprint("Comparing aa with nucleic acids", s.options.verbose, s.options.quiet, 2)
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
                io.save(os.path.join(output_folder, "%s_%s%s.pdb" %(file_name, equivalence[interaction[0]], equivalence[interaction[1]])))
                s.argprint("%s_%s%s.pdb created" %(file_name, equivalence[interaction[0]], equivalence[interaction[1]]), s.options.verbose, s.options.quiet, 2)
                try:
                    interaction = chain_names[number] + chain_names[number + 1]
                    new_model = cp.deepcopy(empty_model)
                    new_model.add(model[chain_names[number]])
                    new_model.add(model[chain_names[number + 1]])
                    io = pdb.PDBIO()
                    io.set_structure(new_model)
                    io.save(os.path.join(output_folder, "%s_%s%s.pdb" %(file_name, equivalence[interaction[0]], equivalence[interaction[1]])))
                    s.argprint("%s_%s%s.pdb created" %(file_name, equivalence[interaction[0]], equivalence[interaction[1]]), s.options.verbose, s.options.quiet, 2)
                except:
                    s.argprint("Example complete", s.options.verbose, s.options.quiet, 0)
