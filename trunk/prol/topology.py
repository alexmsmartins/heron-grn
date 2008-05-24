# coding=UTF-8
"""
prol-evolution GRN, a model for regulatory gene networks

Copyright (C) 2008 Luis Pureza, Oseias Santos, Pedro Martins and
Ricardo Pereira.
 
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import pydot
import networkx as NX
import math
from math import sqrt
from numpy import *

def small_worlds_from_file(path):
    """
    recebe um nome de ficeiro como input e 
    retorna se esta rede Ã© small world ou nÃ£o de
    acordo com os critÃ©rios definidos.
    Um grafo Ã© considerado small-world se Ci for maior
    do que um random graph com o mesmo set de vÃ©rtices
    e se tem uma mÃ©dia das distÃ¢ncias mais curtas baixa.
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
    o mesmo nÃºmero de nÃ³s que o grafo fornecido
    """
    return calc_avg_edge_count(graph) / graph.size()

def calc_random_graph_avg_smallest_path(graph):
    """
    Calcula a mÃ©dia de distÃªncias mais curtas de um random graph com
    o mesmo nÃºmero de nÃ³s que o grafo fornecido
    """
    return log(graph.size()) / log(calc_avg_edge_count(graph))

def calc_avg_graph_shortest_path(graph):
    """
    Calcula a mÃ©dia das distÃ¢ncias mais curtas entre todos os nÃ³s
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
    Calcula o nÃºmero mÃ©dio de ligaÃ§Ãµes por nÃ³
    """
    avg = 0;
    
    for i in graph.degree():
        avg += i
    
    return avg / graph.size() 

def scale_free(graph):
    """
    Recebe o grafo e verifica se Ã© small world
    em caso positivo verifica se existem clusters no grafo
    """
    
    #Array de graus
    array = graph.degree()

    #lista de pontos com o mesmo grau
    list = {}
    for i in array:
        try:
            list[i] = list[i] + 1
        except:
            list[i] = 1

    list_log_x = []
    list_log_y = []
    #list_log = {}
    #avg_listlog_x = 0
    #avg_listlog_y = 0
    for k, v in list.iteritems():
        list_log_x = list_log_x + [math.log(k)]
        list_log_x = list_log_y + [math.log(v)]
        #avg_listlog_x += math.log(k)
        #avg_listlog_y += math.log(v)










    #avg_listlog_x = avg_listlog_x/len(list_log)
    #avg_listlog_y = avg_listlog_y/len(list_log)





    """
    numerator = 0
    denominator = 0
    for k, v in list_log.iteritems():
        numerator = k * v
        denominator = k*k
    
    numerator = numerator - len(list_log)*avg_listlog_x*avg_listlog_y
    denominator = denominator - len(list_log)*avg_listlog_x*avg_listlog_x
    
    
    b = numerator/denominator
    a = avg_listlog_y - b*avg_listlog_x
    """
    coefficients = polyfit(x,y,1)
    
    print quadratic_error(coefficients[0], coefficients[1], list_log)
    
    
def quadratic_error(a, b, list):
    d = 0;
    for x, y in list.iteritems():
        d += pow(y - a*x+b, 2)
    return sqrt(d/len(list))

      
def scale_free_from_file(path):
    """
    recebe um ficheiro para verificar se o grafo Ã© scale free
    """

    return scale_free(dot_to_NXGraph(pydot.graph_from_dot_file(path)))