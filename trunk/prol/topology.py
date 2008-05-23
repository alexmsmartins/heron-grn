# coding=UTF-8

import pydot
import networkx as NX

def small_worlds_from_file(path):
    """
    recebe um nome de ficeiro como input e 
    retorna se esta rede é small world ou não de
    acordo com os critérios definidos.
    Um grafo é considerado small-world se Ci for maior
    do que um random graph com o mesmo set de vértices
    e se tem uma média das distâncias mais curtas baixa.
    """
    return small_worlds(dot_to_NXGraph(pydot.graph_from_dot_file(path)))
    
def small_worlds(graph):
    """
    recebe um objecto NX.Graph como input e 
    retorna o ...
    """
    avgcoef = calc_random_graph_clust_coef(graph)
    above_avg = 0
    
    for coef in NX.clustering(graph):
        if coef > avgcoef:
            above_avg += 1
    
    return (above_avg > graph.size() * 0.5) and (avgcoef > calc_random_graph_clust_coef(graph)) and (calc_avg_graph_shortest_path(graph) < calc_random_graph_avg_smallest_path(graph))

def calc_random_graph_clust_coef(graph):
    """
    Calcula o coeficiente de clustering de um random graph com
    o mesmo número de nós que o grafo fornecido
    """
    return calc_avg_edge_count(graph) / graph.size()

def calc_random_graph_avg_smallest_path(graph):
    """
    Calcula a média de distências mais curtas de um random graph com
    o mesmo número de nós que o grafo fornecido
    """
    return log(graph.size()) / log(calc_avg_edge_count(graph))

def calc_avg_graph_shortest_path(graph):
    """
    Calcula a média das distâncias mais curtas entre todos os nós
    """
    listshortestpaths = NX.path.all_pairs_shortest_path_length(graph)
    avg = 0
    i = 0
    
    for i in listshortestpaths:
        for i in listshortestpaths[i]:
            avg = listshortestpaths[i][j]
            i += 1
    
    return avg / i

def dot_to_NXGraph(dotgraph):
    """
    Converte um grafo pydot.Dot em networkx.Graph
    """
    graph = NX.Graph()
    
    for edge in dotgraph.get_edges():
        graph.add_edge((edge.get_source(), edge.get_destination()))
        
    return graph

def calc_avg_edge_count(graph):
    """
    Calcula o número médio de ligações por nó
    """
    avg = 0;
    
    for i in graph.degree():
        avg += i
    
    return avg / graph.size() 

def scale_free(graph):
    """
    Recebe o grafo e verifica se é small world
    em caso positivo verifica se existem clusters no grafo
    """
    if not small_worlds(graph):
        return False
    
    numHigerClust = 0
    media = NX.average_clustering(graph) * 1.1
    
    for x in NX.clustering(graph):
        if x > media :
            numHigerClust += 1
    
    return (numHigerClust / graph.size()) < 0.1
                
def scale_free_from_file(path):
    """
    recebe um ficheiro para verificar se o grafo é scale free
    """
    return scale_free(dot_to_NXGraph(pydot.graph_from_dot_file(path)))