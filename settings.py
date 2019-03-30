import argparse
import os
import sys

#argparse

main_parser = argparse.ArgumentParser(prog="P.Puzzle",description="Join your protein interaction pairs into a complex(s)",epilog="If you have any doubt, contact us")
input_grup = main_parser.add_mutually_exclusive_group(required=False)
input_grup.add_argument("-f","--folder", default = ".", type = str, metavar = "FOLDER", help = "Folder path with the input PDB files.")
input_grup.add_argument("-l","--list", type = str, dest = "pdb_list", metavar = "PDB", nargs = "*", help = "Input a list of PDB files instead of a folder.")
main_parser.add_argument("-s","--strict", action = "store_false", default = True, help = "Disable the Permisive mode in the PDB parser.")
main_parser.add_argument("-v","--verbose", action = "count", default = 0, help = "Activate the verbose mode, add multiple times for increase the level to a maximum of 3.")
main_parser.add_argument("-q","--quiet", action = "store_true", default = False, help = "Desactivate all the STDERR, STDOUT progress status.")
main_parser.add_argument("-p","--pickle", action = "store_false", default = True, help = "Disable the pickle creation and usage.")
main_parser.add_argument("-r","--radius", metavar = "RADIUS", type=int, default = 2, help = "Empty radious in Amstrongs arround the atom for don't be considered in collision with another, default value is 2.)" )
main_parser.add_argument("-c","--collisions", type=int, metavar = "COLLISIONS_ACCEPTED", default = 30, help = "Maximum number of collisions accepted in each check, default value is 30.")
main_parser.add_argument("-%","--similarity", type=int, metavar = "SIMILARITY_%", choices=range(101), default = 95, help = "Percetage of similarity between sequences with the same name under that, the sequences are considered different, default value is 95.")
main_parser.add_argument("-o","--output", metavar = "OUTPUT", default = ".", type=str, help = "introduce a path where store the output model(s).")
main_parser.add_argument("-i","--initial", type=str, metavar = "INITIAL_CHAIN", help = "Insert manually an initial chain from where the program starts.")

options = main_parser.parse_args()

print(options)
print(options.initial)


#Global variables

def init():
    global similarity
    similarity = {}

#functions

def argprint(printable, verbose_level, quiet, print_level):
    """Prints in the STDERR, if the verbose_level is higer or equal to the print_level"""
    if not quiet:
        if verbose_level >= print_level:
            print(printable,file=sys.stderr)

def Get_PdbList(folder, pdb_list):
    """Check if pdb_list is empty, and return a filled pdb_list"""

    if pdb_list is None:
        onlyfiles = list(os.path.join(folder,f) for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f)))
    else:
        onlyfiles = pdb_list
    return onlyfiles




#Custom error classes

class IncorrectName(Exception):
    """Wrong input naming format"""

    def __init__(self, interaction):
        self.interaction = interaction

    def __str__(self):
        return "Error: %s is captured, this is an incorrect file name format, put interaction chain names (one character maximum for name) between an '_' and the .pdb" %(self.interaction)

class RepeatedChain(Exception):
    """Repeated chain file"""

    def __init__(self, interaction):
        self.interaction = interaction

    def __str__(self):
        return "Error: %s file is repeated, if the copys are not equal, this can generate problems in the final result" %(self.interaction)

class CollisionAppears(Exception):
    """Notifies a collision between chains"""

    def __init__(self, fixed_chain, mobile_chain, collisions):
        self.fixed_chain = fixed_chain
        self.mobile_chain = mobile_chain
        self.collisions = collisions

    def __str__(self):
        return "The number of collisions were greater than the threshold entered when %s is added to %s" %(self.mobile_chain, self.fixed_chain)

    def Get_Collisions(self):
        return self.collisions
