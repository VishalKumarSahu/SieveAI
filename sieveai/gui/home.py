import streamlit as Stream

tab1, tab2, tab3 = Stream.tabs(["Summary", "Create New", "Important Links"])

with tab1:
  _n_projects = 0
  if _n_projects:
    Stream.header("Browse Projects", divider=True)
    f":blue[You have {_n_projects} projects.]"
  else:
    f":red[You have no projects.]"

with tab2:
  Stream.header("Create a new project", divider=True)
  Stream.write("Create a new project by providing the project name and description.")
  name = Stream.text_input("Project Name")
  description = Stream.text_area("Project Description")
  if Stream.button("Create Project"):
    Stream.success(f"Project '{name}' created successfully.")
    _n_projects += 1

with tab3:
  with Stream.popover("Open popover"):
    Stream.markdown("Hello World 👋")
    name = Stream.text_input("What's your name?")

  Stream.write("Your name:", name)

  popover = Stream.popover("Filter items")
  red = popover.checkbox("Show red items.", True)
  blue = popover.checkbox("Show blue items.", True)

  if red:
    Stream.write(":red[This is a red item.]")
  if blue:
    Stream.write(":blue[This is a blue item.]")
