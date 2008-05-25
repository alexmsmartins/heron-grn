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

class MainWindow(wx.Frame):
    """ Program's main window """
    
    def __init__(self):
        """
        Initializes the window
        """
        self.height = 60
        self.width = 700
        wx.Frame.__init__(self, None, size=(self.width, self.height), title='Connect the dots', style=wx.CAPTION | wx.CLOSE_BOX)
        self.CenterOnScreen()

        self.create_menu()
        self.create_toolbar()

        self.open_graph = None

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
                      #                      ('&New...\tCtrl+N', 'Create a new project', 'new', self.OnNewProject),
                      ('&Open...\tCtrl+O', 'Open an existing project', 'open', self.on_file_open),
                      ('','','',''),
                      ('&Quit\tCtrl+Q', 'Terminate the application', None, self.on_file_exit)),
                     ('&View',
                      ('&Statistics', 'View statistics', None, self.on_view_statistics)),
                     
#                     ('&Tools',
#                      ('&Run', 'Runs the current GA', None, self.OnRun),
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
                button = toolbar.AddSimpleTool(-1, icon, label, help)
                window.Bind(wx.EVT_MENU, handler, button)
        
        toolbar_data = (#("New", self.get_icon('new'), "Makes a new project", self.OnNewProject),
                        ("Open", self.get_icon('open'), "Opens a new project", self.on_file_open)),
                        #("", "", "", ""),
                        #("Run", self.get_icon('run'), "Runs the current GA", self.OnRun),
                        #("Config", self.get_icon('config'), "Changes the settings of the current GA", self.OnConfig),
                        #("", "", "", ""),
                        #("Exit", self.get_icon('exit'), "Closes the stuf", self.OnExit))
        
        toolbar = self.CreateToolBar()
        toolbar.SetToolBitmapSize((24,24))
        for button_data in toolbar_data:
            create_toolbar_item(self, toolbar, *button_data)
        toolbar.Realize()


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
            self.open_graph = NX.nx_pydot.read_dot(path)
            
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
        pass

    def get_icon(self, name):
        """
        Returns the icon corresponding to the given id
        """
        icon_path = 'imgs/'
        icon_dict = { 'simulator':      'preferences-desktop-sound.png',
                      ### File 
                      'new':            'document-new.png',
                      'open':           'document-open.png',
                      'save':           'document-save.png',
                      'save_as':                'document-save-as.png',
                      ### Config
                      'run':            'media-playback-start.png',
                      'config':         'stock_chart-edit-type.png',
                      ### Other
                      'exit':           'system-log-out.png'
                    }
    
        return wx.Bitmap(os.path.join(icon_path + icon_dict[name]), wx.BITMAP_TYPE_PNG)

    




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

