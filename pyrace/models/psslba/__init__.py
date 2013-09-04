import os, glob

mods = [ os.path.basename(f)[:-3] for f in glob.glob(os.path.dirname(__file__)+"/*.py") if not os.path.basename(f).startswith('_')]

#mods=[ os.path.basename(f)[:-3] for f in glob.glob(os.path.dirname(__file__)+"/*.py")]
for mod in mods:
    cmd='from .%s import *'%mod
#    print cmd
    exec(cmd)
    

