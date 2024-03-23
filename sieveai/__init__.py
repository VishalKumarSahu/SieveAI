__program__ = "SieveAI"
__description__ = "An automated drug discovery pipeline"
__build__ = "20231229"
__author__ = "Vishal Kumar Sahu"
__url__ = "https://github.com/VishalKumarSahu/SieveAI"
__email__ = "mail@vishalkumarsahu.in"

from .exe import Vina
from .process import Docking
from .cli import dock, rescore
from .lib import ChimeraX, MolConverter, OpenBabel

__all__ = ["__program__", "__description__", "__build__", "__author__", "__url__", "__email__"]


def rock():
    print("Rocking")
