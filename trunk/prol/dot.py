import networkx as NX
import re

def read_graph(filename):
    """
    Reads a networkx graph from a .dot file.
    """
    graph = NX.xdigraph.XDiGraph()
    fx = open(filename, "r")
    edge_regex = re.compile("(.+) -> (.+) \[label=(-?\d+)\]")
    for line in fx.readlines():
        match = edge_regex.search(line)
        if match != None:
            a = match.group(1).strip()
            b = match.group(2).strip()
            w = int(match.group(3).strip())

            graph.add_edge(a, b, w)

    fx.close()
    return graph


def write_graph(graph, filename):
    """
    Writes a graph to a .dot file.
    """
    fx = open(filename, "w")
    fx.write("digraph grn\n{\n")
    for node, out_edges in graph.nodes.items():
        for neighbor in out_edges:
            fx.write("    %s -> %s [label=%d]\n" % \
                         (node.id, neighbor.id, graph.weights[(node, neighbor)]))
            
    fx.write("}")
    fx.close()


def write_nx_graph(graph, filename):
    """
    Writes a networkx graph to a .dot file.
    """
    fx = open(filename, "w")
    fx.write("digraph grn\n{\n")
    for edge in graph.edges():
        fx.write("    %s -> %s [label=%d]\n" % edge)
            
    fx.write("}")
    fx.close()
