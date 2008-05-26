import os
import topology as PRG
import networkx as NX
import pydot
import csv
import sys
import time
import yaml
import random
import simulator as SIM

def read_configuration(filename):
    """
    Read the configuration from some file
    """
    try:
        fx = open(filename, "r")
        config = yaml.load(fx)
        return config
    except IOError, (errno, strerror):
        print "Error while loading the configuration: %s" % strerror
    finally:
        fx.close()

def dump_configuration(config, filename):
    fx = file(filename, 'w')
    yaml.dump(config, fx, default_flow_style=False)
    fx.close()


argv = sys.argv

if len(argv) != 14:
    print "error - args.."
    print "Script.py <file name> <number genes initial> <number genes final> <genes increment> <initial activation probability> <final activation probability> <probability increment> <initial treshold> <final treshold> <treshold increment> <initial inibition> <final inibition> <inibition increment>"
    exit()
else:
    name = argv[1]
    init = int(argv[2])
    end = int(argv[3])
    increment = int(argv[4])
    initProb = float(argv[5])
    endProb = float(argv[6])
    incrementProb = float(argv[7])
    initTreshold = float(argv[8])
    endTreshold = float(argv[9])
    incrementTreshold = float(argv[10])
    initInibition = float(argv[11])
    endInibition = float(argv[12])
    incrementInibition = float(argv[13])
    config_file = "config%d.yaml" % random.randint(0, 10000000)
  
try:        
    os.makedirs("%s" % os.path.join("results", name))
except:
    print "Error!! - Impossible to create folders"

try:
    writer = csv.writer(open(os.path.join("results", name+".csv"), "w"), delimiter=";")
except:
    print "Error!! - Impossible to create cvs file"

try:
    writer.writerow(["Name",
                            "Number Genes",
                            "AVG Edge",
                            "AVG Shortest Path",
                            "RG AVG Smallest Path",
                            "RG Clustering Coefficient",
                            "Genome Creation Time",
                            "Analysis Time",
                            "Genes",
                            "mRNAs",
                            "ncRNAs",
                            "miRNAs"])
except:
    print "Error!! - Impossible to write table headers in cvs"
    
for x in range(init, end+1, increment):
    for t in [initTreshold + ttt*incrementTreshold for ttt in range(int((endTreshold-initTreshold)/incrementTreshold) + 1)]:
        for inib in [initInibition + iiii*incrementInibition for iiii in range(int((endInibition-initInibition)/incrementInibition) + 1)]:
            
            config = read_configuration("config.yaml")
            config["protein/gene binding threshold"] = t
            config["protein inhibition rate"] = inib
            dump_configuration(config, config_file)


            
            inittime = time.time()
            #cmd = popen2.popen4("python heron.py " + str(x) + " " + os.path.join(os.path.join("results", name), str(x) + ".txt") + " -d " + os.path.join(os.path.join("results", name), str(x) + ".dot"))
            #sys.stdout = open('out.log', 'w')
            print "python heron.py " + str(x) + " " + os.path.join(os.path.join("results", name), str(x) + 't' + str(t) + 'i' + str(inib) + ".txt") + " -d " + os.path.join(os.path.join("results", name), (str(x) + 't' + str(t) + 'i' + str(inib)) + ".dot -v -c %s" % config_file)
            bu = os.popen("python heron.py " + str(x) + " " + os.path.join(os.path.join("results", name), str(x) + 't' + str(t) + 'i' + str(inib) + ".txt") + " -d " + os.path.join(os.path.join("results", name), (str(x) + 't' + str(t) + 'i' + str(inib)) + ".dot -v -c %s" % config_file))
        
            aux = bu.readline().split(' ')
            while len(aux) < 5 or aux[4] != "genes\n":
                aux = bu.readline().split(' ')
            genes =  aux[3]
            
            aux = bu.readline().split(' ')
            while len(aux) < 5  or aux[4] != "mRNAs\n":
                aux = bu.readline().split(' ')
            mrnas =  aux[3]
            
            aux = bu.readline().split(' ')
            while len(aux) < 7 or aux[6] != "ncRNAs\n":
                aux = bu.readline().split(' ')
            ncrnas =  aux[5]
            
            aux = bu.readline().split(' ')
            while len(aux) < 5 or aux[4] != "miRNAs\n":
                aux = bu.readline().split(' ')
            mirnas =  aux[3]
            
            
            
            creationfenom = time.time() - inittime
            graph = PRG.dot_to_NXGraph(os.path.join(os.path.join("results", name), (str(x) + 't' + str(t) + 'i' + str(inib)) + ".dot"))
            
            inittime = time.time()
            tablerow = [ name + str(x),
                                str(x),
                                str(PRG.average_degree(graph)),
                                str(PRG.average_shortest_path(graph)),
                                str(PRG.average_shortest_path_random_graph(graph.number_of_nodes(), graph.number_of_edges())),
                                str(PRG.average_clustering_random_graph(graph.number_of_nodes(), graph.number_of_edges())),
                                str(creationfenom),
                                str(time.time() - inittime),
                                genes,
                                mrnas,
                                ncrnas,
                                mirnas]
            
            try:
                writer.writerow(tablerow)
            except:
                print "Error!! - Impossible to write row to cvs"
                
            for y in [initProb + j*incrementProb for j in range(int((endProb-initProb)/incrementProb) + 1)]:
                os.system("python simulator.py " + os.path.join(os.path.join("results", name), str(x) + 't' + str(t) + 'i' + str(inib) + ".txt") + " -p " + str(y) + " -o "+  os.path.join(os.path.join("results", name), str(x) + 't' + str(t) + 'i' + str(inib) + "p" + str(y) + ".png") )

    
