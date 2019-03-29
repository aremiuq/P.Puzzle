
def init():
    global similarity
    similarity = {}

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

class RepeatedChain(Exception):
    """Notifies a collision between chains"""

    def __init__(self, fixed_chain, mobile_chain, collisions):
        self.fixed_chain = fixed_chain
        self.mobile_chain = mobile_chain
        self.collisions = collisions

    def __str__(self):
        return "The number of collisions were greater than the threshold entered when %s is added to %s" %(self.mobile_chain, self.fixed_chain)

    def Get_Collisions(self):
        return self.collisions
