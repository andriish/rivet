"Python utility functions for use by Rivet scripts (and anyone else who wants to)"

def check_python_version(req_version=(2,4,0)):
    "Enforce the Rivet scripts' minimal Python version requirement"
    import sys
    if sys.version_info[:3] < req_version:
        sys.stderr.write( "Python version >= %s is required... exiting\n" % ".".join(req_version) )
        sys.exit(1)

def set_process_name(name):
    "Try to rename the process on Linux so it doesn't appear as 'python <scriptpath>'"
    try:
        import ctypes
        libc = ctypes.cdll.LoadLibrary("libc.so.6")
        libc.prctl(15, name, 0, 0, 0)
    except Exception:
        pass

def import_ET():
    "Try to import the ElementTree XML parser, which has many historical import signatures"
    ET = None
    try:
        import xml.etree.cElementTree as ET
    except ImportError:
        try:
            import cElementTree as ET
        except ImportError:
            try:
                import xml.etree.ElementTree as ET
            except:
                raise ImportError("Can't load the ElementTree XML parser (any of three historical ways)")
    return ET
