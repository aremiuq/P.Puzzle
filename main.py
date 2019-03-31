
#main.py imports

from p_puzzle.ex_generator.generators import *
from p_puzzle.support import settings as s
from p_puzzle.constructor.checker import *
from p_puzzle.constructor.start import *
from p_puzzle.constructor.core import *
from p_puzzle.parser.pdb_files import *
import pickle
import string
import os


#main.py script

if __name__ == "__main__":



    unable_pickle = s.options.pickle
    folder_arg = s.options.folder
    verbose = s.options.verbose
    quiet = s.options.quiet
    initial_chain_arg = s.options.initial
    collisions_accepted_arg = s.options.collisions
    radius_arg = s.options.radius
    percent_arg = s.options.similarity
    pdb_list_arg = s.options.pdb_list
    output_arg = s.options.output
    s.similarity_dict()

    s.argprint("P.Puzzle start", verbose, quiet, 0)

    try:
        level = s.options.exgen_level
    except AttributeError:
        level = False
    if not level == False:
        if level == 1:
            Create_Example(s.options.pdb_file,output_arg)
        elif level == 2:
            Create_Example2(s.options.pdb_file,output_arg)
        exit()

    available_chain_names = list(string.ascii_uppercase + string.ascii_lowercase + string.digits)
    available_chain_names.reverse()

    if unable_pickle:
        relation_pickle = os.path.join(folder_arg, "relationships.p")
        pairs_pickle = os.path.join(folder_arg, "pairs.p")
        try:
            s.argprint("Pickle used!", verbose, quiet, 2)
            relationships = pickle.load(open(relation_pickle,"rb"))
            pairs = pickle.load(open(pairs_pickle,"rb"))
        except:
            s.argprint("No pickle available! Generating a new one...", verbose, quiet, 2)
            onlyfiles = s.Get_PdbList(folder_arg, pdb_list_arg)
            relationships, pairs = Parse_List(onlyfiles)
            out_fd = open(relation_pickle,"wb")
            pickle.dump(relationships ,out_fd)
            out_fd.close()
            out_fd = open(pairs_pickle,"wb")
            pickle.dump(pairs,out_fd)
            out_fd.close()
    else:
        onlyfiles = s.Get_PdbList(folder_arg, pdb_list_arg)
        relationships, pairs = Parse_List(onlyfiles)


    s.argprint("files_parsed", verbose, quiet, 0)
    s.argprint("Number of chains detected: %s" %(len(relationships.keys())), verbose, quiet, 2)

    initial_model = Chose_Start(relationships, pairs, available_chain_names, initial_chain_arg)
    models = Build_model(initial_model, relationships, pairs, collisions_accepted_arg, radius_arg, percent_arg)

    s.argprint("Writing %s pdb file(s)..."%(len(models)), s.options.verbose, s.options.quiet, 0)

    s.argprint("Models created:", s.options.verbose, s.options.quiet, 3)
    s.argprint(models, s.options.verbose, s.options.quiet, 3)


    io = pdb.PDBIO()
    i = 1
    for model in models:
        print(folder_arg)
        filename = "%s_model_%s.pdb" %(i, os.path.basename(folder_arg))
        io.set_structure(model[0])
        io.save(os.path.join(output_arg, filename))
        i += 1

    s.argprint("P.Puzzle finished", s.options.verbose, s.options.quiet, 0)
