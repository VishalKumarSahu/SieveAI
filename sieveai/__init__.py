from .__metadata__ import __version__, __description__, __build__, __name__

def dock():
  from .managers import Manager
  _m = Manager()
  _m.cli_dock()

def server():
  import subprocess, os, sys
  _entrypoint = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gui", "web.py")
  subprocess.run([sys.executable, "-m", "streamlit", "run", _entrypoint])
