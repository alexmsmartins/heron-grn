import os
import topology as PRG
import networkx as NX
import pydot
import csv
import sys
import time
import simulator as SIM

argv = sys.argv

if len(argv) != 5:
    print "error - args.."
    print "Script.py <file name> <number genes initial> <number genes final> <increment>"
    exit()
else:
    name = argv[1]
    init = int(argv[2])
    end = int(argv[3])
    increment = int(argv[4])
    
os.system("mkdir results\\"+name)
writer = csv.writer(open("results\\"+name +".csv", "w"), delimiter=";")
writer.writerow(["Name",
                            "Number Genes",
                            "AVG Edge",
                            "AVG Shortest Path",
                            "RG AVG Smallest Path",
                            "RG Clustering Coefficient",
                            "Is Scale-Free",
                            "Is Small World",
                            "Genome Creation Time",
                            "Analysis Time"])
    
for x in range(init, end+1, increment):
    inittime = time.time()
    os.system("python heron.py " + str(x) + " results\\" + name + "\\" + str(x) + ".txt -d " + " results\\" + name + "\\" + str(x) + ".dot")
    
    
    creationfenom = time.time() - inittime
    graph = PRG.dot_to_NXGraph(pydot.graph_from_dot_file("results\\"+ name + "\\" + str(x) + ".dot"))
    
    inittime = time.time()
    tablerow = [ name + str(x),
                        str(x),
                        str(PRG.calc_avg_edge_count(graph)),
                        str(PRG.calc_avg_graph_shortest_path(graph)),
                        str(PRG.calc_random_graph_avg_smallest_path(graph)),
                        str(PRG.calc_random_graph_clust_coef(graph)),
                        str(PRG.scale_free(graph)),
                        str(PRG.small_worlds(graph)),
                        str(creationfenom),
                        str(time.time() - inittime)]
    writer.writerow(tablerow)
    
