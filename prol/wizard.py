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
import wx.lib.newevent as newevent
import wx.wizard as wiz
import threading
import os
from heron import Heron

# Create an "Advance" event - the event that HeRoN generates when it advances
# one step in the network creation algorithm. By handling this event,
# the GUI is updated.
(HeronAdvanceEvent, EVT_HERON_ADVANCE) = newevent.NewEvent()

class HeronRunner(threading.Thread):
    """
    Runs HeRoN inside a thread
    """

    def __init__(self, mainWindow, genome_size, config):
        """
        Initializes the thread
        """
        threading.Thread.__init__(self)
        self.mainWindow = mainWindow
        self.heron = Heron(genome_size, config)
            
    def run(self):
        """
        Runs HeRoN's graph creation algorithm
        """
        self.heron.create_network(window=self.mainWindow, event=HeronAdvanceEvent)


class ParametersPage(wiz.WizardPageSimple):
    """
    Wizard page where the user chooses the options to generate the network
    """

    def __init__(self, parent):
        """
        Initializes the page
        """
        wiz.WizardPageSimple.__init__(self, parent)

        sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(sizer)
        panel = wx.Panel(self, -1)
        sizer.Add(panel, 0, wx.EXPAND)
       
        self.fill_panel(panel)

        self.Bind(wx.wizard.EVT_WIZARD_PAGE_CHANGING, self.post_processing)

    def fill_panel(self, panel):
        """
        Fill the panel with the necessary widgets
        """

        wx.StaticText(panel, -1, "Genome size:", pos=(8, 10))
        self.txt_genome_size = wx.TextCtrl(self, -1, "10000", pos=(290, 4), size=(90, 27))
        
        wx.StaticText(panel, -1, "Promoter sequence:", pos=(8, 42))
        self.txt_promoter = wx.TextCtrl(self, -1, "0101", pos=(290, 36), size=(90, 27))

        wx.StaticText(panel, -1, "Termination sequence:", pos=(8, 74))
        self.txt_termination = wx.TextCtrl(self, -1, "1111", pos=(290, 68), size=(90, 27))

        wx.StaticText(panel, -1, "U1 left:", pos=(8, 106))
        self.txt_u1_left = wx.TextCtrl(self, -1, "30", pos=(330, 100), size=(50, 27))

        wx.StaticText(panel, -1, "U1 right:", pos=(8, 138))
        self.txt_u1_right = wx.TextCtrl(self, -1, "13", pos=(330, 132), size=(50, 27))

        wx.StaticText(panel, -1, "miRNA/mRNA binding site size:", pos=(8, 170))
        self.spin_miRNA_size = wx.SpinCtrl(self, -1, "6", pos=(330, 164), size=(50, 27))
        self.spin_miRNA_size.SetRange(1, 100)

        wx.StaticText(panel, -1, "Protein inhibition rate:", pos=(8, 202))
        self.spin_protein_inhibition = wx.SpinCtrl(self, -1, "50", pos=(330, 196), size=(50, 27))
        self.spin_protein_inhibition.SetRange(1, 100)
        wx.StaticText(panel, -1, "%", pos=(385, 202))

        wx.StaticText(panel, -1, "Protein/gene binding site size:", pos=(8, 234))
        self.spin_pg_binding_size =  wx.SpinCtrl(self, -1, "6", pos=(330, 228), size=(50, 27))
        self.spin_pg_binding_size.SetRange(1, 100)

        wx.StaticText(panel, -1, "Protein/gene binding threshold:", pos=(8, 266))
        self.spin_pg_threshold =  wx.SpinCtrl(self, -1, "20", pos=(330, 260), size=(50, 27))
        self.spin_pg_threshold.SetRange(1, 50)
        wx.StaticText(panel, -1, "%", pos=(385, 266))

        self.functions = { 
            "max" : "max", 
            "min" : "min", 
            "average" : "lambda list: sum(list) / len(list)",
            "random" : "random.choice"
        }
        wx.StaticText(panel, -1, "Binding function:", pos=(8, 298))
        self.rb_functions = wx.RadioBox(self, -1, "", (10, 302), wx.DefaultSize, self.functions.keys(), \
                                            -1, wx.NO_BORDER)

    def post_processing(self, event):
        """
        Called when the user tries to move to the next page
        """
        errors = self.validate()
        if len(errors) > 0:
            dlg = wx.MessageDialog(self, "The following errors were found:\n\n - " \
                                       + "\n - ".join(errors), "Erros found", wx.OK)
            dlg.ShowModal()
            dlg.Destroy()
            event.Veto()
        else:
            self.set_configuration()
            event.Allow()
 
    def validate(self):
        """
        Validate user's input
        """

        def is_valid_sequence(sequence):
            """
            Checks if some sequence is valid, i.e, every element is
            0, 1, 2 or 3
            """
            try:
                numbers = map(int, sequence)
                if len(numbers) == 0:
                    return False
                return all(x >= 0 and x <= 3 for x in numbers)
            except:
                return False

        errors = []
        try:
            genome_size = int(self.txt_genome_size.GetValue())
            if genome_size <= 0:
                errors.append("The genome's size must be a positive integer")
        except:
                errors.append("The genome's size must be a positive integer")

        if not is_valid_sequence(self.txt_promoter.GetValue()):
            errors.append("Invalid promoter")

        if not is_valid_sequence(self.txt_termination.GetValue()):
            errors.append("Invalid termination sequence")

        if not is_valid_sequence(self.txt_u1_left.GetValue()):
            errors.append("Invalid U1 left sequence")

        if not is_valid_sequence(self.txt_u1_right.GetValue()):
            errors.append("Invalid U1 right sequence")

        return errors

    def set_configuration(self):
        """
        Save the configuration in a hashtable
        """
        self.GetParent().genome_size = int(self.txt_genome_size.GetValue())
        self.GetParent().config = {
                "miRNA/mRNA binding site size" : self.spin_miRNA_size.GetValue(),
                "protein/gene binding site size" : self.spin_pg_binding_size.GetValue(),
                "protein/gene binding threshold" : self.spin_pg_threshold.GetValue() / 100.0,
                "U1 left" : self.txt_u1_left.GetValue(),
                "U1 right" : self.txt_u1_right.GetValue(),
                "promoter" : self.txt_promoter.GetValue(),
                "termination" : self.txt_termination.GetValue(),
                "binding function" : self.functions[self.rb_functions.GetStringSelection()],
                "protein inhibition rate" : self.spin_protein_inhibition.GetValue() / 100.0
        }


class ProgressPage(wiz.WizardPageSimple):
    """
    Wizard page that shows progress while the network is created
    """

    def __init__(self, parent):
        """
        Initializes the page
        """
        wiz.WizardPageSimple.__init__(self, parent)

        sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(sizer)
        panel = wx.Panel(self, -1)
        sizer.Add(panel, 0, wx.EXPAND)
       
        self.fill_panel(panel)
        self.Bind(wx.wizard.EVT_WIZARD_PAGE_CHANGED, self.run_heron)
        self.Bind(EVT_HERON_ADVANCE, self.on_heron_advance)

    def fill_panel(self, panel):
        """
        Fill the panel with the necessary widgets
        """

        self.steps = [
            wx.StaticText(panel, -1, "1. Create the genome", (8, 10)),
            wx.StaticText(panel, -1, "2. Splice the genes", (8, 32)),
            wx.StaticText(panel, -1, "3. Translate the mRNAs into proteins", (8, 54)),
            wx.StaticText(panel, -1, "4. Create miRNAs", (8, 76)),
            wx.StaticText(panel, -1, "5. Create connections between proteins and genes", (8, 98)),
            wx.StaticText(panel, -1, "6. Create connections between miRNAs and mRNAs", (8, 120))
        ]

        bold = wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD)
        self.steps[0].SetFont(bold)

        wx.StaticText(panel, -1, "Progress:", pos=(8, 160)).SetFont(bold)
        self.progress = wx.Gauge(self, -1, 50, (8, 180), (385, 25))
        self.progress.SetRange(100)

        wx.StaticText(panel, -1, "Statistics:", (8, 230)).SetFont(bold)
        wx.StaticText(panel, -1, "- Number of genes:", (8, 250))
        wx.StaticText(panel, -1, "- Number of mRNAs:", (8, 270))
        wx.StaticText(panel, -1, "- Number of proteins:", (8, 290))
        wx.StaticText(panel, -1, "- Number of ncRNAs:", (8, 310))
        wx.StaticText(panel, -1, "- Number of miRNAs:", (8, 330))
        wx.StaticText(panel, -1, "- Number of protein/gene bindings:", (8, 350))
        wx.StaticText(panel, -1, "- Number of miRNA/mRNA bindings:", (8, 370))

        self.lbl_statistics = {
            "genes" : wx.StaticText(panel, -1, "", (340, 250), (50, 20), wx.ALIGN_RIGHT),
            "mRNAs" : wx.StaticText(panel, -1, "", (340, 270), (50, 20), wx.ALIGN_RIGHT),
            "proteins" : wx.StaticText(panel, -1, "", (340, 290), (50, 20), wx.ALIGN_RIGHT),
            "ncRNAs" : wx.StaticText(panel, -1, "", (340, 310), (50, 20), wx.ALIGN_RIGHT),
            "miRNAs" : wx.StaticText(panel, -1, "", (340, 330), (50, 20), wx.ALIGN_RIGHT),
            "protein/gene bindings" : wx.StaticText(panel, -1, "", (340, 350), (50, 20), wx.ALIGN_RIGHT),
            "miRNA/mRNA bindings" : wx.StaticText(panel, -1, "", (340, 370), (50, 20), wx.ALIGN_RIGHT)
        }

    def run_heron(self, event):
        """
        Creates a regulatory network on HeRoN with the parameters defined
        by the user
        """
        self.progress.SetValue(0)
        wx.FindWindowById(wx.ID_FORWARD).Disable()

        [label.SetLabel("") for label in self.lbl_statistics.values()]
        thread = HeronRunner(self, self.GetParent().genome_size, \
                                 self.GetParent().config)
        self.GetParent().heron = thread.heron
        thread.start()


    def on_heron_advance(self, event):
        bold = wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD)
        regular = wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL)

        if event.step > 0:
            self.steps[event.step - 1].SetFont(regular)

        if event.step < len(self.steps):
            self.steps[event.step].SetFont(bold)
        else:
            wx.FindWindowById(wx.ID_FORWARD).Enable()
            

        if hasattr(event, "genes"):
            self.lbl_statistics["genes"].SetLabel(str(event.genes))

        if hasattr(event, "mRNAs"):
            self.lbl_statistics["mRNAs"].SetLabel(str(event.mRNAs))

        if hasattr(event, "ncRNAs"):
            self.lbl_statistics["ncRNAs"].SetLabel(str(event.ncRNAs))

        if hasattr(event, "proteins"):
            self.lbl_statistics["proteins"].SetLabel(str(event.proteins))

        if hasattr(event, "miRNAs"):
            self.lbl_statistics["miRNAs"].SetLabel(str(event.miRNAs))

        if hasattr(event, "protein_gene_bindings"):
            self.lbl_statistics["protein/gene bindings"].SetLabel(str(event.protein_gene_bindings))

        if hasattr(event, "miRNA_mRNA_bindings"):
            self.lbl_statistics["miRNA/mRNA bindings"].SetLabel(str(event.miRNA_mRNA_bindings))

        if hasattr(event, "done"):
            self.progress.SetValue(event.done * 100)

            
class SavePage(wiz.WizardPageSimple):
    """
    Wizard page that lets the user save the network
    """

    def __init__(self, parent):
        """
        Initializes the page
        """
        wiz.WizardPageSimple.__init__(self, parent)

        sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(sizer)
        panel = wx.Panel(self, -1)
        sizer.Add(panel, 0, wx.EXPAND)
       
        self.fill_panel(panel)

        self.Bind(wx.wizard.EVT_WIZARD_PAGE_CHANGED, self.on_changed)
        self.Bind(wx.wizard.EVT_WIZARD_PAGE_CHANGING, self.post_processing)

    def fill_panel(self, panel):
        """
        Fill the panel with the necessary widgets
        """

        wx.StaticText(panel, -1, "Choose the location where you want to save the network:", (8, 8))
        
        self.path_button = wx.Button(panel, -1, "...", pos=(8, 32), size=(385, 34))
        self.Bind(wx.EVT_BUTTON, self.on_choose_path, self.path_button)

        self.chk_save_as_dot = wx.CheckBox(self, -1, "Also save as .dot file", pos=(8, 74))

    def on_choose_path(self, event):
        """
        Opens a "save file" dialog so that the user can choose the location
        where he wants to save the network.
        """
        dlg = wx.FileDialog(self, message="Save file as ...", \
                                defaultDir=os.getcwd(), defaultFile="", \
                                wildcard="HeRoN graphs (*.her)|*.her", \
                                style=wx.SAVE)

        if dlg.ShowModal() == wx.ID_OK:
            dlg.Destroy()
            path = dlg.GetPath()
            if path[-4:] != ".her":
                path += ".her"

            self.path_button.SetLabel(path)
            wx.FindWindowById(wx.ID_FORWARD).Enable()

    def on_changed(self, event):
        """
        Called when the page is displayed
        """
        if self.path_button.GetLabel() == "...":
            wx.FindWindowById(wx.ID_FORWARD).Disable()

    def post_processing(self, event):
        """
        Saves the network
        """
        path = self.path_button.GetLabel()
        save_as_dot = self.chk_save_as_dot.IsChecked()
        heron = self.GetParent().heron
        heron.dump(path)

        if save_as_dot:
            heron.save_as_dot(path + ".dot")
        

        
class Wizard(object):
    """
    The wizard
    """
   
    def __init__(self):
        """
        Initializes the wizard
        """
        self.wizard = wiz.Wizard(None, -1, "HeRoN network creation wizard")
        self.wizard.SetPageSize(wx.Size(402, 416))

        self.page1 = ParametersPage(self.wizard)
        self.page2 = ProgressPage(self.wizard)
        self.page3 = SavePage(self.wizard)

        wiz.WizardPageSimple_Chain(self.page1, self.page2)
        wiz.WizardPageSimple_Chain(self.page2, self.page3)

    def show(self):
        """
        Shows the wizard
        """
        self.wizard.RunWizard(self.page1)
 
        self.wizard.Destroy()
        self.wizard.DestroyChildren()


class Application(wx.App):
    """
    Wizard application
    """

    def OnInit(self):
        """
        Creates a new wizard
        """
        wizard = Wizard()
        wizard.show()

        return True
        
if __name__ == '__main__':
    app = Application(False)
    app.MainLoop()

