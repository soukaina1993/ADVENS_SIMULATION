import os
import sys
import html.parser
import html.entities
import unittest.mock
from shiny._main import main
path = os.path.dirname(os.path.abspath(__file__))
apath = os.path.join(path, "app.py")

# these next two lines are only if you are using Windows OS
drive, apath = os.path.splitdrive(apath)
apath = apath.replace("\\","/")
#

sys.argv = ['shiny', 'run', apath]
main()
