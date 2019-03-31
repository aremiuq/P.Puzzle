
#settings.py imports
import argparse
import os
import sys

#argparse

main_parser = argparse.ArgumentParser(prog="P.Puzzle",description="Join your protein interaction pairs into a complex(s)",epilog="If you have any doubt, contact us.")
input_grup = main_parser.add_mutually_exclusive_group(required=False)
input_grup.add_argument("-f","--folder", default = ".", type = str, metavar = "FOLDER", help = "Folder path with the input PDB files. Current folder is used as a default.)")
input_grup.add_argument("-l","--list", type = str, dest = "pdb_list", metavar = "PDB", nargs = "*", help = "Input a list of PDB files instead of a folder.")
info_grup = main_parser.add_mutually_exclusive_group(required=False)
info_grup.add_argument("-v","--verbose", action = "count", default = 0, help = "Activate the verbose mode, add multiple times to increase the level to a maximum of 3.")
info_grup.add_argument("-q","--quiet", action = "store_true", default = False, help = "Deactivate all the STDERR, STDOUT progress status.")
main_parser.add_argument("-s","--strict", action = "store_false", default = True, help = "Disable the Permissive mode in the PDB parser.")
main_parser.add_argument("-p","--pickle", action = "store_false", default = True, help = "Disable the pickle creation and usage.")
main_parser.add_argument("-r","--radius", metavar = "RADIUS", type=int, default = 2, help = "Empty radius in Amstrongs around the atom that aren't considered in collision with another. Default value is 2.)" )
main_parser.add_argument("-c","--collisions", type=int, metavar = "COLLISIONS_ACCEPTED", default = 30, help = "Maximum number of collisions accepted in each check. Default value is 30.")
main_parser.add_argument("-%","--similarity", type=int, metavar = "SIMILARITY_%", choices=range(101), default = 95, help = "Percentage of similarity between sequences with the same name. If less than this percentage, the sequences are considered different. Default value is 95.")
main_parser.add_argument("-o","--output", metavar = "OUTPUT", default = ".", type=str, help = "Introduce a path where the program will store the output model(s). Current folder is used as default.")
main_parser.add_argument("-i","--initial", type=str, metavar = "INITIAL_CHAIN", help = "Insert manually an initial chain from where the program starts.")

sub_parser = main_parser.add_subparsers(help = "sub-command help")
exgen_parser = sub_parser.add_parser("exgen", help = "Example generator mode")
exgen_parser.add_argument("pdb_file", type=str, metavar = "PDB_FILE", help = "PDB file for split in related pairs" )
exgen_parser.add_argument("exgen_level", type=int, default = 1, choices=range(1,3), help = "1 = Split a pdb file into related pairs giving to each chain a diferent name, 2 = Split a pdb file into related pairs giving the same name to similar chains" )

options = main_parser.parse_args()

#Global variable

def similarity_dict():
    global similarity
    similarity = {}

#settings.py functions

def argprint(printable, verbose_level, quiet, print_level):
    """Prints in the STDERR, if the verbose_level is higer or equal to the print_level"""
    if not quiet:
        if verbose_level >= print_level:
            print(printable,file=sys.stderr)

def Get_PdbList(folder, pdb_list):
    """Check if pdb_list is empty, and return a filled pdb_list"""

    if pdb_list is None:
        onlyfiles = list(os.path.join(folder,f) for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f)) and os.path.join(folder,f)[-4:] == ".pdb")
    else:
        onlyfiles = pdb_list
    return onlyfiles

def Check_folder(folder):
    """If folder doesn't exist, create it in the actual path

    Input:
    folder = folder name
    """
    path = os.getcwd()

    new_path = os.path.join(path,folder)

    if not os.path.isdir(new_path):

        try:
            os.mkdir(new_path)
        except OSError:
            argprint("Creation of the directory %s failed" % path, options.verbose, options.quiet, 1)
            exit()
        else:
            argprint("Successfully created the directory %s " % path, options.verbose, options.quiet, 1)

    return True

#Custom error classes

class IncorrectName(Exception):
    """Wrong input naming format"""

    def __init__(self, interaction):
        self.interaction = interaction

    def __str__(self):
        return "Error: %s is captured. This is an incorrect file name format. Please put the interaction chain names (one character maximum for name) between '_' and the .pdb" %(self.interaction)

class RepeatedChain(Exception):
    """Repeated chain file"""

    def __init__(self, interaction):
        self.interaction = interaction

    def __str__(self):
        return "Error: %s file is repeated. If the copies are not equal, problems can occur in the final result" %(self.interaction)

class CollisionAppears(Exception):
    """Notifies a collision between chains"""

    def __init__(self, fixed_chain, mobile_chain, collisions):
        self.fixed_chain = fixed_chain
        self.mobile_chain = mobile_chain
        self.collisions = collisions

    def __str__(self):
        return "The number of collisions was greater than the threshold entered when %s is added to %s" %(self.mobile_chain, self.fixed_chain)

    def Get_Collisions(self):
        return self.collisions
