import os
import sys
from shiny._main import main
path = os.path.dirname(os.path.abspath(__file__))
apath = os.path.join(path, "file1.py")

# these next two lines are only if you are using Windows OS
drive, apath = os.path.splitdrive(apath)
apath = apath.replace("\\","/")
#

sys.argv = ['shiny', 'run', apath]
main()