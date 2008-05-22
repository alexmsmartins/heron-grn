from pydot import *
import networkx as NX


def small_worlds_from_file(dot_file):
    """recebe um nome de ficeiro como input e 
    retorna se esta rede é small world ou não de
    acordo com os critérios definidos.
    Um grafo é considerado small-world se Ci for maior
    do que um random graph com o mesmo set de vértices
    e se tem uma média das distências mais curtas baixa."""
    small_worlds(read(path))
    
    
    
    
def small_worlds(graph):
    """recebe um objecto Dot como input e 
    retorna o """
    
    

def clustering_coefficient_from_file():
    """recebe um nome de ficeiro como input e 
    retorna o queficiente de clustering do
    respectivo grafo"""
    average_clustering(read(path))
    
def clustering_coefficient(graph):
    """recebe um obejcto Dot como input e 
    retorna o queficiente de clustering do
    respectivo grafo"""
    average_clustering(graph)




def random_graph_generator(nNodes):
    graph = Graph("random graph")