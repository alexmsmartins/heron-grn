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

import wx
import os
import networkx as NX
import topology
import math
import re
import dot
import random
import graph
import pickle
from simulator import Simulator
from grn_element import GRNElement
import Gnuplot, Gnuplot.funcutils


class StatisticsWindow(wx.Dialog):
    """ Statistics window """
    
    def __init__(self, graph):
        """
        Initializes the window
        """
        self.graph = graph
        self.height = 270
        self.width = 370
        wx.Dialog.__init__(self, None, size=(self.width, self.height), \
                               title='Statistics', \
                               style = wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)

        sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(sizer)

        panel = wx.Panel(self, -1)
        sizer.Add(panel, 1, wx.EXPAND)
        sizer = wx.BoxSizer(wx.VERTICAL)
        panel.SetSizer(sizer)
        self.CenterOnScreen()

        self.create_list(panel)      
        self.fill_list()

    def create_list(self, panel):
        self.list = wx.ListCtrl(panel, -1, style=wx.LC_REPORT)
        self.list.InsertColumn(0, "Property")
        self.list.InsertColumn(1, "Value", wx.LIST_FORMAT_RIGHT)
        panel.GetSizer().Add(self.list, 1, wx.EXPAND)

    def fill_list(self):
        def check_small_worlds_conditions(graph):
            """
            Checks that n > k > ln(n) > 1

            See the paper "Collective dynamics of 'small-world' networks", 
            by Duncan J. Watts & Steven H. Strogatz
            """
            n = graph.number_of_nodes()
            k = topology.average_degree(graph)
            return n > k > math.log(n) > 1

        statistics = [("Number of nodes (n)", lambda graph: graph.number_of_nodes()),
                      ("Number of edges", lambda graph: graph.number_of_edges()),
                      ("Average degree (k)", topology.average_degree),
                      ("n > k > ln(n) > 1", check_small_worlds_conditions),
                      ("Directed?", lambda graph: graph.is_directed()),
                      ("Clustering coefficient", NX.average_clustering),
                      ("Clustering coefficient (random)", lambda graph: \
                           topology.average_clustering_random_graph(graph.number_of_nodes(), graph.number_of_edges())),
                      ("Average shortest path", topology.average_shortest_path),
                      ("Average shortest path (random)", lambda graph: \
                           topology.average_shortest_path_random_graph(graph.number_of_nodes(), graph.number_of_edges()))
                      ]  

        for label, func in statistics:
            self.list.Append([label, "%s" % func(self.graph)])

        self.list.SetColumnWidth(0, 230)
        self.list.SetColumnWidth(1, 120)


class PlotOptionsWindow(wx.Dialog):
    """ Plot options window """
    
    def __init__(self, graph):
        """
        Initializes the window
        """
        self.graph = graph
        self.height = 145
        self.width = 350
        wx.Dialog.__init__(self, None, size=(self.width, self.height), \
                               title='Plot options')

        sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(sizer)

        panel = wx.Panel(self, -1)
        sizer.Add(panel, 1, wx.EXPAND)
        sizer = wx.BoxSizer(wx.VERTICAL)
        panel.SetSizer(sizer)
        self.CenterOnScreen()

        self.fill_panel(panel)


    def fill_panel(self, panel):
        wx.StaticText(panel, -1, "Plot only nodes that match", pos=(8, 14))
        self.txt_filter = wx.TextCtrl(panel, -1, ".*", pos=(190, 8), size=(150, 27))
        self.chk_log_log = wx.CheckBox(panel, -1, "log/log plot", pos=(8, 40))

        wx.StaticText(panel, -1, "Plot distribution of node", pos=(8, 72))
        self.rb_edge_direction = wx.RadioBox(panel, -1, "", (165, 52), wx.DefaultSize, ["inputs", "outputs"], 1, style=wx.NO_BORDER)

        
        self.btn_ok = wx.Button(panel, wx.ID_OK, "OK", pos=(165, 102))
        self.Bind(wx.EVT_BUTTON, self.on_ok, self.btn_ok)
        btn_cancel = wx.Button(panel, wx.ID_CANCEL, "Cancel", pos=(255, 102))


    def on_ok(self, event):
        def logaritmize(data):
            def surreal_log(x):
                return x and math.log(x) or 0
            
            return [(surreal_log(x), surreal_log(y)) for x, y in data]
            
        regex = re.compile(self.txt_filter.GetValue(), re.I)
        plot_log = self.chk_log_log.IsChecked()
        plot_outputs = self.rb_edge_direction.GetSelection() == 1

        connectivity = {}
        data = {}
        for edge in self.graph.edges():
            node = plot_outputs and edge[0] or edge[1]
            if not regex.search(node):
                continue

            if connectivity.has_key(node):
                data[connectivity[node]] -= 1
            connectivity[node] = connectivity.setdefault(node, 0) + 1
            data[connectivity[node]] = data.setdefault(connectivity[node], 0) + 1

        tuples = data.items()

        if plot_log:
            tuples = logaritmize(tuples)

        chart = Gnuplot.Gnuplot(persist=1)
        chart.title("Node connectivity")
        chart.xlabel("Number of %s" % ((plot_outputs and "outputs" or "inputs") + (plot_log and " (log)" or "")))
        chart.ylabel("Number of nodes%s" % (plot_log and " (log)" or ""))
        chart("set xrange [0:]")
        chart("set yrange [0:]")
        chart("set style fill solid border -1")
        chart("plot '-' with %s" % (plot_log and "points" or "boxes"))
        chart("%s" % "\n".join(["%s %s" % (a, b) for (a, b) in tuples]))
        chart("e")

        self.Close()


class SimulationOptionsWindow(wx.Dialog):
    """ Simulation options window """
    
    def __init__(self, graph, graph_path):
        """
        Initializes the window
        """
        self.graph = graph
        self.graph_path = graph_path
        self.height = 245
        self.width = 350
        wx.Dialog.__init__(self, None, size=(self.width, self.height), \
                               title='Simulation options')

        sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(sizer)

        panel = wx.Panel(self, -1)
        sizer.Add(panel, 1, wx.EXPAND)
        sizer = wx.BoxSizer(wx.VERTICAL)
        panel.SetSizer(sizer)
        self.CenterOnScreen()

        self.fill_panel(panel)


    def fill_panel(self, panel):
        wx.StaticText(panel, -1, "Number of steps:", pos=(8, 14))
        self.spin_steps = wx.SpinCtrl(panel, -1, "100", pos=(273, 8), size=(50, 27))
        self.spin_steps.SetRange(1, 100000)

        wx.StaticText(panel, -1, "Initial node activation probability:", pos=(8, 46))
        self.spin_act_prob = wx.SpinCtrl(panel, -1, "50", pos=(273, 40), size=(50, 27))
        self.spin_act_prob.SetRange(1, 100)
        wx.StaticText(panel, -1, "%", pos=(328, 46))

        wx.StaticText(panel, -1, "Plot nodes of type:", pos=(8, 78))

        self.node_types = {
            "Genes" : "Gene",
            "mRNA" : "MessengerRNA",
            "Proteins" : "Protein",
            "miRNA" : "MicroRNA"
        }

        self.chk_types = []
        for (i, item) in enumerate(self.node_types.items()):
            self.chk_types.append(wx.CheckBox(panel, -1, item[0], pos=(12, 100 + i * 25)))
            self.chk_types[-1].SetValue(True)

        self.btn_ok = wx.Button(panel, wx.ID_OK, "OK", pos=(165, 202))
        self.Bind(wx.EVT_BUTTON, self.on_ok, self.btn_ok)
        btn_cancel = wx.Button(panel, wx.ID_CANCEL, "Cancel", pos=(255, 202))

    def on_ok(self, event):
        steps = self.spin_steps.GetValue()
        probability = self.spin_act_prob.GetValue() / 100.0
        types = [self.node_types[chk.GetLabel()] for chk in self.chk_types if chk.IsChecked()]

        grn = self.open_grn_graph()
        simulator = Simulator(grn, probability)
        simulator.simulate_network(steps, types)

        self.Close()

    def open_grn_graph(self):
        """
        Asks the user for the path to the GRN graph
        """
        test_path = self.graph_path[0:-4]
        if os.path.exists(test_path):
            path = test_path
        else:        
            dialog = wx.FileDialog(self, message="Select the graph",
                                   wildcard="HeRoN files (*.her)|*.her|" \
                                       "All files (*.*)|*.*",
                                   style=wx.OPEN | wx.CHANGE_DIR)

            if dialog.ShowModal() == wx.ID_OK:
                path = dialog.GetPath()

        try:
            fx = open(path, "r")
            return pickle.load(fx)
        except IOError, (errno, strerror):
            print "Error while reading the graph file: %s" % strerror


    def build_grn_graph(self):
        """
        Converts from a networkx graph to a graph.graph
        """
        grn = graph.graph()
        for a, b, w in self.graph.edges():
            first = GRNElement(None, a)
            second = GRNElement(None, b)
            if not grn.nodes.has_key(first):
                grn.add_nodes([first])

            if not grn.nodes.has_key(second):
                grn.add_nodes([second])

            grn.add_arrow(first, second, w)
        
        assert len(grn.nodes) == self.graph.number_of_nodes()
        return grn


class MainWindow(wx.Frame):
    """ Program's main window """
    
    def __init__(self):
        """
        Initializes the window
        """
        self.height = 600
        self.width = 700
        wx.Frame.__init__(self, None, size=(self.width, self.height), title='Connect the dots', style=wx.CAPTION | wx.CLOSE_BOX)
        self.CenterOnScreen()

        self.create_menu()
        self.create_toolbar()
        self.create_status_bar()

        self.open_graph = None
        self.on_file_open(None)

    def create_menu(self):
        """
        Initializes the menu bar
        """
       
        def create_submenu(items):
            """
            Creates a single submenu
            """
            menu = wx.Menu()
            for label, status, image, handler in items:
                if not label:
                    menu.AppendSeparator()
                    continue

                menu_item = menu.Append(-1, label, status)
                if image != None:
                    menu_item.SetBitmap(self.get_icon(image))

                self.Bind(wx.EVT_MENU, handler, menu_item)
                    
            return menu

        menu_data = (("&File",
                    #    ('&New...\tCtrl+N', 'Create a new graph...', 'new', self.on_file_new),
                      ('&Open...\tCtrl+O', 'Open an existing project', 'open', self.on_file_open),
                      ('','','',''),
                      ('&Quit\tCtrl+Q', 'Terminate the application', None, self.on_file_exit)),
                     ('&View',
                      ('&Statistics', 'View statistics', None, self.on_view_statistics),
                      ('&Visualize graph...', 'Visualize the graph', None, self.on_view_visualize_graph)),                     
                     ('&Plot',
                      ('&Node connectivity...', 'Node connectivity', None, self.on_plot_node_connectivity)),
                     
                     ('&Tools',
                      ('&Run simulation...', 'Simulates the graph with the HeRoN algorithm', None, self.on_tools_run)),
#                      ('&Config', 'Changes the settings of the current GA', None, self.OnConfig)
#                      ), 
#                     ('&About',
#                      ('&About', 'Info about the program', None, self.OnAbout)
                     )
        
        menu_bar = wx.MenuBar()
        for submenu_data in menu_data:
            label = submenu_data[0]
            items = submenu_data[1:]
            menu_bar.Append(create_submenu(items), label)

        self.SetMenuBar(menu_bar)

    
    def create_toolbar(self):
        """
        Creates the toolbar
        """
        
        def create_toolbar_item(window, toolbar, label, icon, help, handler):
            if not label:
                toolbar.AddSeparator()
            else:
                button = toolbar.AddLabelTool(-1, label, icon, wx.NullBitmap)
                window.Bind(wx.EVT_MENU, handler, button)
        
        toolbar_data = (("Open", self.get_icon('open'), "Opens a new project", self.on_file_open),
                        ("", "", "", ""),
                        ("View statistics", self.get_icon('stats'), "View statistics", self.on_view_statistics),
                        ("Plot connectivity", self.get_icon('plot'), "Plot graph connectivity", self.on_plot_node_connectivity),
                        ("", "", "", ""),
                        ("Run", self.get_icon('run'), "Simulates the network", self.on_tools_run))
        
        self.toolbar = self.CreateToolBar()
        self.toolbar.SetToolBitmapSize((24,24))
        for button_data in toolbar_data:
            create_toolbar_item(self, self.toolbar, *button_data)

        self.toolbar.Realize()


    def create_status_bar(self):
        statusBar = self.CreateStatusBar()
        statusBar.SetFieldsCount(2)
        statusBar.SetStatusWidths([-3, -1])



    def get_icon(self, name):
        """
        Returns the icon corresponding to the given id
        """
        icon_path = 'imgs/'
        icon_dict = {
                      'open':           'document-open.png',
                      'run':            'emblem-system.png',
                      'plot':           'utilities-system-monitor.png',
                      'stats':          'document-properties.png'
                    }
    
        return wx.Bitmap(os.path.join(icon_path + icon_dict[name]), wx.BITMAP_TYPE_PNG)


    def on_file_new(self, event):
        """
        Event handler for the "New" action
        """
        random_graph2 = NX.generators.random_graphs.watts_strogatz_graph(800, 70, 0.8)

        print NX.average_clustering(random_graph2)

        #random_graph = NX.generators.random_graphs.gnm_random_graph(800, 30 * 1000)
        random_graph = self.my_small_world(800, 70, 0.8)
        graph = NX.xdigraph.XDiGraph()

        for edge in random_graph.edges():
            first, second, w = edge
            if random.random() < 0.5:
                second, first, w = edge

            graph.add_edge("Gene_" + str(first), "Gene_" + str(second), random.choice([-1, 1]))
          
#        assert graph.number_of_nodes() == random_graph.number_of_nodes()
#        assert graph.number_of_edges() == random_graph.number_of_edges()
#        self.open_graph = graph


    def my_small_world(self, n, k, p):
        """ 
        Return a Watts-Strogatz small world graph. 
  
        First create a ring over n nodes.  Then each node in the ring is 
        connected with its k nearest neighbors.  Then shortcuts are 
        created by rewiring existing edges as follows: for each edge u-v 
        in the underlying "n-ring with k nearest neighbors"; with 
        probability p replace u-v with a new edge u-w with 
        randomly-chosen existing node w. In contrast with 
        newman_watts_strogatz_graph(), the random rewiring does not 
        increase the number of edges. 
        
        
        :Parameters: 
        - `n`: the number of nodes 
        - `k`: each node is connected to k neighbors in the ring topology 
        - `p`: the probability of rewiring an edge 
        - `seed`: seed for random number generator (default=None) 
        
        """ 
        G = NX.generators.classic.empty_graph(n,create_using=NX.xdigraph.XDiGraph())
        

        G.name="watts_strogatz_graph(%s,%s,%s)"%(n,k,p) 
        nlist = G.nodes() 
        fromv = nlist 
        # connect the k/2 neighbors 
        for n in range(1, k/2+1):
            tov = fromv[n:] + fromv[0:n] # the first n are now last 
            for i in range(len(fromv)): 
                G.add_edge(fromv[i], tov[i]) 
#                G.add_edge(tov[i], fromv[i]) 
         # for each edge u-v, with probability p, randomly replace with 
         # edge u-w 
        e = G.edges() 
        for (u, v, w) in e: 
            if random.random() < p: 
                newv = random.choice(nlist) 
                # avoid self-loops and reject if edge u-newv exists 
                # s that the correct WS model? 
                while newv == u or G.has_edge(u, newv):  
                    newv = random.choice(nlist) 
                G.delete_edge(u,v)  # conserve number of edges  
                G.add_edge(u, newv, w) 


        return G             

    def on_file_open(self, event):
        """
        Event handler for the "Open" action
        """
        dialog = wx.FileDialog(self, message="Select the graph",
                               wildcard="DOT files (*.dot)|*.dot|" \
                                   "All files (*.*)|*.*",
                               style=wx.OPEN | wx.CHANGE_DIR)

        if dialog.ShowModal() == wx.ID_OK:
            path = dialog.GetPath()
            self.open_graph = dot.read_graph(path) 
            self.open_graph_path = path
            self.SetTitle(path)
            for i in range(self.toolbar.GetToolsCount()):
                self.toolbar.EnableTool(i, True)
        else:
            exit()

        dialog.Destroy()    

   


    def on_file_exit(self, event):
        """
        Event handler for when the user wants to quit the application
        """
        exit(0)


    def on_view_statistics(self, event):
        """
        Show the statistics window
        """
        dialog = StatisticsWindow(self.open_graph)
        dialog.Show()


    def on_view_visualize_graph(self, event):
        """
        Visualize the graph
        """
        while True:
            filename = "/tmp/output%d" % random.randint(1, 10000)
            if not os.path.exists(filename):
                break;

        dot.write_nx_graph(self.open_graph, filename)
        os.system("dotty %s" % filename)
        os.remove(filename)


    def on_plot_node_connectivity(self, event):
        """
        Opens the dialog to let the user choose plot options
        """
        dialog = PlotOptionsWindow(self.open_graph)
        dialog.Show()


    def on_tools_run(self, event):
        """
        Show the simulation options window
        """
        dialog = SimulationOptionsWindow(self.open_graph, self.open_graph_path)
        dialog.Show()


class Application(wx.App):
    """
    ConnectDots application
    """

    def OnInit(self):
        """
        Creates the main window
        """
        window = MainWindow()
        window.Show()

        return True
        
if __name__ == '__main__':
    app = Application(False)
    app.MainLoop()

