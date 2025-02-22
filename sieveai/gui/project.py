import streamlit as Stream

Stream.header("Project", divider=True)

_Manager = Stream.session_state['sieveai']

from io import StringIO
s = StringIO()
print(_Manager.SETTINGS.user, file=s)
result = s.getvalue()
Stream.success(result)

with Stream.expander(f"Project Settings", expanded=False):
    _Manager.path_base = Stream.text_input(f"Project Path: {_Manager.path_base}", _Manager.path_base)
    _Manager.SETTINGS.dir_receptors = Stream.text_input(f"Receptors: {_Manager.SETTINGS.dir_receptors}", _Manager.SETTINGS.dir_receptors)
    _Manager.SETTINGS.dir_receptors = Stream.text_input(f"Ligands: {_Manager.SETTINGS.dir_receptors}", _Manager.SETTINGS.dir_receptors)

    for _key_val in _Manager.SETTINGS.user.items():
        Stream.write(_key_val)

_n_receptors = 50
with Stream.expander(f"Select Receptors (n={_n_receptors})", expanded=True):
  "Not discovered yet"
  if Stream.button('Add 50 Receptors'):
    _n_receptors += 50
    f"rec {_n_receptors}"

with Stream.expander("Select Ligands", expanded=True):
  "Not discovered yet"

with Stream.expander("Select Complexes", expanded=True):
  "Choose the receptor and ligand combination"
  df = _Manager.DF(
      [
        {"Receptor": "st.selectbox", "Ligand": 4, "Dock": True},
        {"Receptor": "st.balloons", "Ligand": 5, "Dock": False},
        {"Receptor": "st.time_input", "Ligand": 3, "Dock": True},
    ]
  )
  edited_df = Stream.data_editor(df)
  Stream.write(edited_df)

if Stream.button('Start Docking'):
    # _Manager.handle_process()
    _info = Stream.info('Starting docking...')
    _Manager.time_sleep(3) # Wait for 3 seconds
    _info.empty() # Clear the alert
