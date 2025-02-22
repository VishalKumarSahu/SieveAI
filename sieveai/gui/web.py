import streamlit as Stream
import time

from sieveai.managers.manager import Manager

class SieveAIDashboard:
  def __init__(self):
    Stream.session_state['sieveai'] = Manager(
      log_level='debug',
      Report_Interval=(6,),
      Report_Check_Interval=(30, 'seconds')
    )
    self.log_placeholder = None

  def layout(self):
    Stream.title("SieveAI Dashboard")
    self.setup_sidebar()
    self.setup_footer()

  def setup_footer(self):
    log_ph = self.log_placeholder = (f"🕑 {time.strftime('%Y-%m-%d %H:%M:%S')}")
    footer_html = f"""<div style='text-align: center; background-color: lightblue;'>
      <p>Developed with ❤️ by TheBiomics | {log_ph}</p>
    </div>"""
    Stream.markdown(footer_html, unsafe_allow_html=True)

  def setup_sidebar(self):
    _entrypoints = {
        "home.py": {"title": "SieveAI", "icon": "😄", "default": True},
        "project.py": {"title": "Project", "icon": "😄", "default": False},
        "projects.py": {"title": "Projects", "icon": ":material/dashboard:", "default": False},
        "settings.py": {"title": "Settings", "icon": ":material/settings:", "default": False},
      }

    _eps = []
    for _py, _props in _entrypoints.items():
      _epr = Stream.Page(_py, **_props)
      _eps.append(_epr)

    _nv = Stream.navigation(_eps)
    # _nv.run()

    with Stream.sidebar:
      Stream.title("SieveAI")
      # Stream.sidebar.page_link('./home.py', label="SieveAI", icon="😄")

      if Stream.button("version " + Stream.session_state['sieveai'].__version__):
        self.version_details()
        # Stream.session_state['sieveai'].__version__

      # "[Visit GitHub](https://github.com/VishalKumarSahu/SieveAI)"

      with Stream.spinner("Loading..."):
        time.sleep(2)
      Stream.toast("Loaded...")

  @Stream.dialog("SieveAI")
  def version_details(self, *args, **kwargs):
    Stream.write(f"SieveAI version: {Stream.session_state['sieveai'].__version__}")

  def run(self):
    # Stream.set_page_config(page_title="Data manager", page_icon=":material/edit:")
    self.layout()

if __name__ == "__main__":
  dashboard = SieveAIDashboard()
  dashboard.run()
