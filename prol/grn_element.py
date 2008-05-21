class GRNElement(object):
    """
    Base class for all GRN graph elements
    """

    def __init__(self, parent):
        """
        Creates a new GRN element with the given parent

        Arguments:
        - `parent`: The parent node of this element
        """
        self.parent = parent
        
