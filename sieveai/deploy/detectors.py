"""Deploy helpers - detectors and installers."""

import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional


def which(exe: str) -> Optional[Path]:
    """Find executable in PATH."""
    path = shutil.which(exe)
    return Path(path) if path else None


def check_exe(exe: str) -> bool:
    """Check if executable exists."""
    return which(exe) is not None


def check_pkg(pkg: str) -> bool:
    """Check if Python package is importable."""
    try:
        __import__(pkg.replace("-", "_"))
        return True
    except ImportError:
        return False


def install_package(package: str) -> bool:
    """
    Install Python package using pip.
    
    Returns:
        True if successful, False otherwise.
    """
    try:
        subprocess.check_call([
            sys.executable, "-m", "pip", "install", package
        ])
        return True
    except subprocess.CalledProcessError:
        return False


def get_install_guide(exe: str) -> str:
    """Get installation instructions for executable."""
    guides = {
        "vina": """
AutoDock VINA:
  - Ubuntu: Download from https://github.com/ccsb-scripps/AutoDock-Vina/releases
  - macOS: brew install autodock-vina
  - conda: conda install -c conda-forge autodock-vina
""",
        "chimerax": """
ChimeraX:
  - Download from https://www.rbvi.ucsf.edu/chimerax/download.html
  - Add to PATH after installation
""",
        "prepare_receptor": """
MGLTools (ADFR Suite):
  - Download from http://ccsb.scripps.edu/mgltools/downloads/
  - Install ADFR Suite
  - Add to PATH
""",
        "obabel": """
OpenBabel:
  - Ubuntu: sudo apt install openbabel
  - macOS: brew install open-babel
  - conda: conda install -c conda-forge openbabel
""",
    }
    
    return guides.get(exe, f"Install {exe} from official source")
