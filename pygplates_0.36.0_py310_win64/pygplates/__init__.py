# Import the pygplates shared library (C++).
from .pygplates import *
# Import any private symbols (with leading underscore).
from .pygplates import __version__
from .pygplates import __doc__

# Let the pygplates shared library (C++) know of its imported location.
import os.path
pygplates._post_import(os.path.dirname(__file__))

# Rename '.pygplates' to '._impl' so we don't have both pygplates.<symbol> and pygplates.pygplates.<symbol>.
_impl = pygplates
del pygplates
