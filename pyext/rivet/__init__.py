"Python interface to the Rivet MC validation system"

# ## Try to bootstrap the Python path for rivet loading
# import os, sys
# PROGPATH = sys.argv[0]
# PROGNAME = os.path.basename(PROGPATH)
# import commands
# try:
#     modname = sys.modules[__name__].__file__
#     binpath = os.path.dirname(modname)
#     rivetconfigpath = os.path.join(binpath, "rivet-config")
#     rivetpypath = commands.getoutput(rivetconfigpath + " --pythonpath")
#     sys.path.append(rivetpypath)
# except:
#     pass

## Change dlopen status to GLOBAL for Rivet lib
import sys
try:
    import ctypes
    sys.setdlopenflags(sys.getdlopenflags() | ctypes.RTLD_GLOBAL)
    del ctypes
except:
    import dl
    sys.setdlopenflags(sys.getdlopenflags() | dl.RTLD_GLOBAL)
    del dl
del sys

## Import SWIG-generated wrapper core
from rivet.rivetwrap import *

## Import plot info helper
from rivet.plotinfo import *

## Import submodules into the visible namespace
import spiresbib, util

## Try to use Psyco optimiser if available
# TODO: or PyPy, Cython, Weave?
try:
    import psyco
    psyco.full()
except Exception, e:
    pass
