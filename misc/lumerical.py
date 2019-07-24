import platform
import imp

# Location of Python API depending of the OS in use
if platform.system() == "Linux":
    lumapi = imp.load_source("lumapi", "/opt/lumerical/2019b/api/python/lumapi.py")
elif platform.system() == "Windows":
    lumapi = imp.load_source("lumapi", "C:\\Program Files\\Lumerical\\2019b\\api\\python\\lumapi.py")
elif platform.system() == "Darwin":
    lumapi = imp.load_source("lumapi", "/Applications/Lumerical 2019b.app/Contents/API/Python/lumapi.py")

if __name__ == "__main__":
    
    #from lumerical import lumapi
    mode1 = lumapi.FDTD()