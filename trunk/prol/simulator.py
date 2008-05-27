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
from optparse import OptionParser, OptionValueError
import Gnuplot, Gnuplot.funcutils

class Simulator(object):
    """
    Simulator object
    """

    def __init__(self, grn, probability):
        """
        Initializes the simulator
        """
        self.grn = grn
        self.__initialize_network(probability)
  

    def simulate_network(self, steps, types, img_output=None, verbose=False):
        """
        Simulate the network and show the behavior of the desired nodes
        on gnuplot
        """
        chart = Gnuplot.Gnuplot(persist=1)
        chart.title("Network dynamics")
        chart.xlabel(",".join(types))
        chart.ylabel("Steps")
        chart("set pointsize 0.3")

        if img_output != None:
            chart("set terminal png")
            chart("set output '%s'" % img_output)
        data = []
        for step in range(steps):
            if verbose:
                print " * Step %3d of %d" % (step + 1, steps)

            chart("set yrange[0:%d]" % (step + 1))
            data += [(node, step) for node in self.__simulate_one_step(types)]

            if img_output == None and len(data) > 0:
                chart.plot(data)

        if len(data) > 0:
            chart.plot(data)


    def __initialize_network(self, probability):
        """
        Activates some genes choosen randomly.
        """
        for node in self.grn.get_nodes():
            if node.id[0:node.id.find("_")] == "Gene" and random.random() <= probability:
                node.enabled = True


    def __simulate_one_step(self, types):
        """
        Simulate the execution of the GRN for one step
        """

        data = []
        on = []
        off = []
        j = 0
        for node in self.grn.get_nodes():
            if node.id[0:node.id.find("_")] in types:
                if node.enabled:
                    data.append(j)
                j += 1

            # For each neighbor...
            for neighbor in self.grn.get_node(node):
                weight = self.grn.weights[(node, neighbor)]

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


def parse_args():
    """
    Parse command-line arguments
    """

    def parse_types(option, opt, value, parser):
	"""
	Parses the -t command-line option
	"""
	type_str = parser.rargs[0]
	try:
	    types = type_str.split(",")
	    setattr(parser.values, option.dest, types)
	except:
	    raise OptionValueError("Unknown type")
	del parser.rargs[0]


    usage = "usage: %prog [options] <input>"
    parser = OptionParser(usage=usage)
    parser.add_option("-p", "--probability", dest="probability", type="float",
		      default=0.5, help="Initial gene activation probability")
    parser.add_option("-o", "--output", dest="output",
		      help="Save the resulting graphic as an image")
    parser.add_option("-s", "--steps", dest="steps", type="int", default=100,
                      help="Number of steps")
    parser.add_option("-t", "--types", dest="types", action="callback",
		      callback=parse_types, default=["Gene"],
		      help="Types to plot in the graph, separated by commas.")
    parser.add_option("-v", action="store_true", dest="verbose",
                      help="Verbose mode")
    
    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.error("Incorrent number of arguments")
        exit()

    input_file = args[0]

    return (options, input_file)


if __name__ == '__main__':
    (options, input_file) = parse_args()

    try:
        fx = open(input_file, "r")
        grn = pickle.load(fx)
    except IOError, (errno, strerror):
        print "Error while reading the graph file: %s" % strerror
        exit()

    simulator = Simulator(grn, options.probability)
    simulator.simulate_network(options.steps, options.types, img_output=options.output, verbose=options.verbose)
