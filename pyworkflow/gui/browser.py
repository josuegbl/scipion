# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Jose Gutierrez (jose.gutierrez@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
In this module a simple ObjectBrowser is implemented.
This class can be subclasses to extend its functionality.
A concrete use of ObjectBrowser is FileBrowser, where the
elements to inspect and preview are files.
"""

import os
import os.path
import stat

import Tkinter as tk
import ttk

import xmipp
import gui
from pyworkflow.utils import dirname, getHomePath, prettySize, getExt, dateStr
from tree import BoundTree, TreeProvider
from text import TaggedText
from pyworkflow.utils.properties import Icon


class ObjectBrowser(tk.Frame):
    """ This class will implement a simple object browser.
    Basically, it will display a list of elements at the left
    panel and can display a preview and description on the
    right panel for the selected element.
    An ObjectView will be used to grab information for
    each element such as: icon, preview and description.
    A TreeProvider will be used to populate the list (Tree).
    """
    def __init__(self, parent, treeProvider, showPreview=True, **args):
        tk.Frame.__init__(self, parent, **args)
        self.treeProvider = treeProvider
        self._lastSelected = None
        gui.configureWeigths(self)
        # The main layout will be two panes, 
        # At the left containing the elements list
        # and the right containing the preview and description
        p = tk.PanedWindow(self, orient=tk.HORIZONTAL)
        p.grid(row=0, column=0, sticky='news')
        
        leftPanel = tk.Frame(p)
        self._fillLeftPanel(leftPanel)
        p.add(leftPanel, padx=5, pady=5)
        p.paneconfig(leftPanel, minsize=300)
        
        if showPreview:
            rightPanel = tk.Frame(p)
            gui.configureWeigths(rightPanel)
            self._fillRightPanel(rightPanel)
            p.add(rightPanel, padx=5, pady=5)    
            p.paneconfig(rightPanel, minsize=200)    
        
            # Register a callback when the item is clicked
            self.tree.itemClick = self._itemClicked
        
    def _fillLeftPanel(self, frame):
        gui.configureWeigths(frame)
        self.tree = BoundTree(frame, self.treeProvider)
        self.tree.grid(row=0, column=0, sticky='news')
        self.itemConfig = self.tree.itemConfig
        self.getImage = self.tree.getImage
    
    def _fillRightPanel(self, frame):
        top = tk.Frame(frame)
        top.grid(row=0, column=0, sticky='news')
        gui.configureWeigths(top)
        top.rowconfigure(0, minsize=200)
        self._fillRightTop(top)
        
        bottom = tk.Frame(frame)
        bottom.grid(row=1, column=0, sticky='news')
        gui.configureWeigths(bottom)
        bottom.rowconfigure(1, weight=1)
        self._fillRightBottom(bottom)
        
    def _fillRightTop(self, top):
        self.noImage = self.getImage('no-image128.png')
        self.label = tk.Label(top, image=self.noImage)
        self.label.grid(row=0, column=0, sticky='news')
        
    def _fillRightBottom(self, bottom):
        self.text = TaggedText(bottom, width=40, height=15, bg='white')
        self.text.grid(row=0, column=0, sticky='news')
        
    def _itemClicked(self, obj):
        self._lastSelected = obj
        img, desc = self.treeProvider.getObjectPreview(obj)
        self.text.clear()
        if isinstance(img, str):
            img = self.getImage(img)
        if img is None:
            img = self.noImage
        self.label.config(image=img)
        if desc is not None:
            self.text.addText(desc)
            
    def getSelected(self):
        """ Return the selected object. """
        return self._lastSelected
      

#------------- Classes and Functions related to File browsing --------------

class FileInfo(object):
    """ This class will store some information about a file.
    It will serve to display files items in the Tree.
    """
    def __init__(self, path, filename):
        self._fullpath = os.path.join(path, filename)
        self._filename = filename
        self._stat = os.stat(self._fullpath)
        
    def isDir(self):
        return stat.S_ISDIR(self._stat.st_mode)
    
    def getFileName(self):
        return self._filename
    
    def getPath(self):
        return self._fullpath
    
    def getSize(self):
        """ Return a human readable string of the file size."""
        return prettySize(self._stat.st_size)
    
    def getDate(self):
        return dateStr(self._stat.st_mtime)
    
    
class FileHandler(object):
    """ This class will be used to get the icon, preview and info
    from the different types of objects.
    It should be used with FileTreeProvider, where different
    types of handlers can be registered.
    """
    def getFileIcon(self, objFile):
        """ Return the icon name for a given file. """
        if objFile.isDir():
            icon = 'file_folder.gif'
        else:
            icon = 'file_generic.gif'
        
        return icon
    
    def getFilePreview(self, objFile):
        """ Return the preview image and description for the specific object. """
        if objFile.isDir():
            return 'fa-folder-open.png', None
        return None, None
    
    
class TextFileHandler(FileHandler):   
    def __init__(self, textIcon):
        FileHandler.__init__(self)
        self._icon = textIcon
         
    def getFileIcon(self, objFile):
        return self._icon
    
class MdFileHandler(FileHandler):
    def getFileIcon(self, objFile):
        return 'file_md.gif'
    
class SqlFileHandler(FileHandler):
    def getFileIcon(self, objFile):
        return 'file_sqlite.gif'    
    
class ImageFileHandler(FileHandler):
    _image = xmipp.Image()
    _index = ''
    
    def getFilePreview(self, objFile):
        fn = self._index + objFile.getPath()
        dim = 128
        self.tkImg = gui.getTkImage(self._image, fn, dim)
        
        return self.tkImg, None 
    
class ParticleFileHandler(ImageFileHandler):
    def getFileIcon(self, objFile):
        return 'file_image.gif'
    
class VolFileHandler(ImageFileHandler):
    def getFileIcon(self, objFile):
        return 'file_vol.gif'
    
class StackHandler(ImageFileHandler):
    _index = '1@'
    
    def getFileIcon(self, objFile):
        return 'file_stack.gif'
    
    
class FileTreeProvider(TreeProvider):
    """ Populate a tree with files and folders of a given path """
    
    _FILE_HANDLERS = {}
    _DEFAULT_HANDLER = FileHandler()
    
    @classmethod
    def registerFileHandler(cls, fileHandler, *extensions):
        """ Register a FileHandler for a given file extention. 
        Params:
            fileHandler: the FileHandler that will take care of extensions.
            *extensions: the extensions list that will be associated to this FileHandler.
        """
        for fileExt in extensions:
            cls._FILE_HANDLERS[fileExt] = fileHandler
        
    def __init__(self, currentDir=None, showHidden=False):
        self._currentDir = os.path.abspath(currentDir)
        self._showHidden = showHidden
        self.getColumns = lambda: [('File', 300), ('Size', 70), ('Time', 150)]
    
    def getFileHandler(self, obj):
        filename = obj.getFileName()
        fileExt = getExt(filename)
        return self._FILE_HANDLERS.get(fileExt, self._DEFAULT_HANDLER)
        
    def getObjectInfo(self, obj):
        filename = obj.getFileName()
        fileHandler = self.getFileHandler(obj)
        icon = fileHandler.getFileIcon(obj)
        
        info = {'key': filename, 'text': filename, 
                'values': (obj.getSize(), obj.getDate()), 'image': icon
                }
            
        return info
    
    def getObjectPreview(self, obj):
        fileHandler = self.getFileHandler(obj)
        return fileHandler.getFilePreview(obj)
    
    def getObjectActions(self, obj):
        return []
    
    def getObjects(self):
        files = os.listdir(self._currentDir)
        files.sort()
        for f in files:
            if self._showHidden or not f.startswith('.'):
                yield FileInfo(self._currentDir, f)

    def getDir(self):
        return self._currentDir
    
    def setDir(self, newPath):
        self._currentDir = newPath
        
# Some constants for the type of selection
# when the file browser is opened
SELECT_NONE = 0 # No selection, just browse files                  
SELECT_FILE = 1
SELECT_FOLDER = 2
SELECT_PATH = 3 # Can be either file or folder


class FileBrowser(ObjectBrowser):
    """ The FileBrowser is a particular class of ObjectBrowser
    where the "objects" are just files and directories.
    """
    def __init__(self, parent, initialDir='.', 
                 selectionType=SELECT_FILE, selectionSingle=True, 
                 allowFilter=True, filterFunction=None, previewDim=144,
                 showHidden=False):
        """ 
        """
        tp = FileTreeProvider(initialDir, showHidden)
        ObjectBrowser.__init__(self, parent, tp)
        
        if selectionType != SELECT_NONE:
            buttonsFrame = tk.Frame(self)
            self._fillButtonsFrame(buttonsFrame)
            buttonsFrame.grid(row=1, column=0)

    def _fillLeftPanel(self, frame):
        """ Redefine this method to include a buttons toolbar and
        also include a filter bar at the bottom of the Tree.
        """
        # Tree with files
        frame.columnconfigure(0, weight=1)
        
        treeFrame = tk.Frame(frame)
        ObjectBrowser._fillLeftPanel(self, treeFrame)
        # Register the double-click event
        self.tree.itemDoubleClick = self._itemDoubleClick
        treeFrame.grid(row=1, column=0, sticky='news')
        # Toolbar frame
        toolbarFrame = tk.Frame(frame)
        self._fillToolbar(toolbarFrame)
        toolbarFrame.grid(row=0, column=0, sticky='new')
        # Filter frame
        tk.Label(frame, text="Filter").grid(row=2, column=0)
        
        frame.rowconfigure(1, weight=1)

    def _fillToolbar(self, frame):
        """ Fill the toolbar frame with some buttons. """
        self._col = 0
        
        def addButton(text, image, command):
            btn = tk.Label(frame, text=text, image=self.getImage(image), 
                       compound=tk.LEFT, cursor='hand2')
            btn.bind('<Button-1>', command)
            btn.grid(row=0, column=self._col, sticky='nw',
                     padx=(0, 5), pady=5)
            self._col += 1
            
        addButton('Refresh', Icon.ACTION_REFRESH, self._actionRefresh)
        addButton('Home', Icon.HOME, self._actionHome)
        addButton('Back', Icon.ARROW_LEFT, self._actionUp)
        addButton('Up', Icon.ARROW_UP, self._actionUp)
        
    def _fillButtonsFrame(self, frame):
        """ Add button to the bottom frame if the selectMode
        is distinct from SELECT_NONE.
        """
        tk.Button(frame, text="Close", image=self.getImage(Icon.BUTTON_CLOSE),
                  command=self._close,
                  compound=tk.LEFT).grid(row=0, column=0, padx=(0,5))                        
        tk.Button(frame, text="Select", image=self.getImage(Icon.BUTTON_SELECT),
                  command=self._select, 
                  compound=tk.LEFT).grid(row=0, column=1)
                
    def _actionRefresh(self, e=None):
        self.tree.update()
        
    def _goDir(self, newDir):
        self.treeProvider.setDir(newDir)
        self.tree.update()
        
    def _actionUp(self, e=None):
        self._goDir(dirname(self.treeProvider.getDir()))
        
    def _actionHome(self, e=None):
        self._goDir(getHomePath())
        
    def _itemDoubleClick(self, obj):
        if obj.isDir():
            self._goDir(obj.getPath())
            
    def onClose(self):
        pass
    
    def onSelect(self, obj):
        pass
    
    def _close(self, e=None):
        self.onClose()
        
    def _select(self, e=None):
        self.onSelect(self.getSelected())
        
        
class BrowserWindow(gui.Window):
    """ Windows to hold a browser frame inside. """
    def __init__(self, title, master=None, **args):
        if 'minsize' not in args:
            args['minsize'] = (800, 400)
        gui.Window.__init__(self, title, master, **args)
        
    def setBrowser(self, browser):
        browser.grid(row=0, column=0, sticky='news')
        self.itemConfig = browser.tree.itemConfig
        
class FileBrowserWindow(BrowserWindow):
    """ Windows to hold a browser frame inside. """
    def __init__(self, title, master=None, path=None, **args):
        BrowserWindow.__init__(self, title, master, **args)
        self.registerHandlers()
        browser = FileBrowser(self.root, path)
        def selected(obj):
            print obj.getPath()
        browser.onClose = self.close
        browser.onSelect = selected
        self.setBrowser(browser) 
        
    def registerHandlers(self):
        FileTreeProvider.registerFileHandler(TextFileHandler('file_text.gif'), '.txt', '.log', '.out', '.err')
        FileTreeProvider.registerFileHandler(TextFileHandler('file_python.gif'), '.py')
        FileTreeProvider.registerFileHandler(TextFileHandler('file_java.gif'), '.java')
        FileTreeProvider.registerFileHandler(MdFileHandler(), '.xmd', '.star')
        FileTreeProvider.registerFileHandler(SqlFileHandler(), '.sqlite')
        FileTreeProvider.registerFileHandler(ParticleFileHandler(), '.mrc', '.spi')
        FileTreeProvider.registerFileHandler(VolFileHandler(), '.vol')
        FileTreeProvider.registerFileHandler(StackHandler(), '.stk', '.spi')
    
