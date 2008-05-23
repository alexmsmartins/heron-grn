#! /usr/bin/python

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

from gene import *
import pickle
import random
from optparse import OptionParser
import Gnuplot, Gnuplot.funcutils

def parse_args():
    """
    Parse command-line arguments
    """
    usage = "usage: %prog [options] <input>"
    parser = OptionParser(usage=usage)
    parser.add_option("-s", "--steps", dest="steps", type="int", default=100,
                      help="Number of steps")
    parser.add_option("-v", action="store_true", dest="verbose",
                      help="Verbose mode")
    
    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.error("Incorrent number of arguments")
        exit()

    input_file = args[0]

    return (options, input_file)


def initialize_network(grn, probability):
    """
    Activates some genes choosen randomly.
    """
    for node in grn.get_nodes():
        if isinstance(node, Gene) and random.random() <= probability:
            node.enabled = True
    

def simulate_network(grn, steps, options):
    """
    Simulate the execution of the GRN for a number of steps
    """

    data = []
    for i in range(steps):
        if options.verbose:
            print " * Step %3d of %d" % (i + 1, steps)
        on = []
        off = []
        j = 0
        for node in grn.get_nodes():
            if isinstance(node, Gene):
                if node.enabled:
                    data.append((j, i))
                j += 1

            # For each neighbor...
            for neighbor in grn.get_node(node):
                weight = grn.weights[(node, neighbor)]

                if node.enabled == True:
                    if weight > 0:
                        # Activate the neighbor
                        on.append(neighbor)
                    else:
                        # Repress the neighbor
                        off.append(neighbor)
                else:
                    # Any active element turns inactive if its activator
                    # is not active
                    if weight > 0 and neighbor.enabled:
                        off.append(neighbor)

        for element in on:
            element.enabled = True

        # Repression takes precedence and so is applied after activation
        for element in off:
            element.enabled = False

    return data


if __name__ == '__main__':
    (options, input_file) = parse_args()

    try:
        fx = open(input_file, "r")
        grn = pickle.load(fx)
    except IOError, (errno, strerror):
        print "Error while reading the graph file: %s" % strerror
        exit()

    chart = Gnuplot.Gnuplot()
    chart.title("Network")
    chart.xlabel("Gene")
    chart.ylabel("Step")
    chart("set pointsize 0.1")
    chart("set yrange[0:%d]" % options.steps)

    initialize_network(grn, 0.5)
    data = simulate_network(grn, options.steps, options)

    chart.plot(data)

    raw_input("Press enter to continue...")
