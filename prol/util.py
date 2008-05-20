def complement (n):
    """ Convert 0 <-> 1 and 2 <-> 3 """ 
    return abs(1 - n * ((n % 2) * 2 - 1))

def complement_string (sequence):
    """ Complements an entire sequence """
    
    return "".join([str(complement(int(n))) for n in sequence])
