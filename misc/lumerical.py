import platform, imp

# Location of Python API depending of the OS in use
if platform.system() == "Linux":
    lumapi = imp.load_source("lumapi", str("/opt/lumerical/fdtd/api/python/lumapi.py"))
elif platform.system() == "Windows":
    lumapi = imp.load_source("lumapi", "C:\\Program Files\\Lumerical\\FDTD\\api\\python\\lumapi.py")
elif platform.system() == "Darwin":
    lumapi = imp.load_source("lumapi", "/Applications/Lumerical 2020a.app/Contents/API/Python/lumapi.py")

# Methods used with the Lumerical Python API
def putOptions(options, lumApp=None):
    " Reads an option structure and for each field it puts the name in the Lumerical environment. "
    if lumApp is None:
        print('Error: No active lumerical software to pass the options.')
        return
    for name, value in options.items():
        lumApp.putv(str(name), float(value))
        print(type(name))
        print(type(value))


if __name__ == "__main__":
    
    #from lumerical import lumapi
    mode1 = lumapi.FDTD()