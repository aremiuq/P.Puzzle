import Bio.PDB

pdb_parser = Bio.PDB.PDBParser(PERMISSIVE=True)

pdb_file = pdb_parser.get_structure("pdb_file", "./example1/pair_his3_sc_XA.pdb")

class chain1(Bio.PDB.Select):
    def accept_chain(self, chain):
        if chain.get_id()=="A":
            return True
        else:
            return False

class chain2(Bio.PDB.Select):
    def accept_chain(self, chain):
        if chain.get_id()=="B":
            return True
        else:
            return False

io = Bio.PDB.PDBIO()

io.set_structure(pdb_file)

io.save("chain1.pdb", chain1())

io.save("chain2.pdb", chain2())

