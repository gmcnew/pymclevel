from distutils.core import setup
import py2exe

setup(
    console = ['bestofboth.py'],
    options = {
        "py2exe": {
            "unbuffered": True,
            "optimize": 2,
            "dll_excludes": ["mswsock.dll", "powrprof.dll", "w9xpopen.exe"],
            "excludes": ["Tkconstants", "Tkinter", "tcl", "socket", "ssl", "bz2", "ctypes", "unicodedata", "pyexpat"],
        }
    }
)