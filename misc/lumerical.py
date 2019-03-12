import platform
import imp

# Location of Python API depending of the OS in use
if platform.system() == "Linux":
    lumapi = imp.load_source("lumapi", "/opt/lumerical/FDTD/api/python/lumapi.py")
elif platform.system() == "Windows":
    lumapi = imp.load_source("lumapi", "C:\\Program Files\\Lumerical\\FDTD\\api\\python\\lumapi.py")
elif platform.system() == "Darwin":
    lumapi = imp.load_source("lumapi", "/Applications/Lumerical/DEVICE/DEVICE.app/Contents/API/Python/lumapi.py")

class lumerical(object):

    def close_session(self):
        """ Close a session of Lumerical MODE. """
        self.inst.close()

class mode(lumerical):

    def __init__(self):
        self.inst = lumapi.MODE()

class fdtd(object):

    def __init__(self):
        self.inst = lumapi.FDTD()

class device(object):

    def __init__(self):
        self.inst = lumapi.DEVICE()

class interconnect(object):

    def __init__(self):
        self.inst = lumapi.INTERCONNECT()