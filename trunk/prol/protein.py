class Protein(object):
    """ A protein. """

    def __init__ (self, sequence):
        self.sequence = sequence

    def __str__ (self):
        """ Builds a string representation of this protein """
        return str(self.sequence)
