import os as OS
import math as MATH

import pandas as PD

import warnings as WARNINGS

from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from .ExecutableBase import ExecutableBase

class VinaBase(ExecutableBase):
  def __init__(self, *args, **kwargs):
    WARNINGS.simplefilter('ignore', PDBConstructionWarning)
    super(VinaBase, self).__init__(**kwargs)
    self.__defaults = {
        "path_base": None,
        "path_config": None,

        # Allowed Extensions
        "ext_ligand": ["*.pdb", "*.sdf", "*.mol2"],
        "ext_receptor": ["*.pdb"],

        # Setting Extensions
        "ext__contacts": "contacts.txt",
        "ext__hbonds": "hbonds.txt",
        "ext__chimera_cxc": "cxc",
        "ext__dock_result": "result.pdbqt",
        "ext__dock_log": "result.log",

        "dir_ligand": "ligand",
        "dir_receptor": "receptor",
        "dir_docking": "docking",

        "dir_receptor_pdbqt": "receptor-pdbqt",
        "dir_receptor_clean": "receptor-clean",
        "dir_receptor_summary": "receptor-summary",
        "dir_receptor_config": "receptor-config",
        "dir_ligand_pdbqt": "ligand-pdbqt",
        "dir_analysis": "vina-analysis",

        "path_receptor_pdbqt": None,
        "path_receptor_clean": None,
        "path_receptor_summary": None,
        "path_receptor_config": None,
        "path_ligand_pdbqt": None,
        "path_analysis": None,

        "path_docking": None,
        "path_receptor": None,
        "path_ligand": None,

        "path_mgltools": None,

        "processes": ["prepare_ligand", "prepare_receptor", "prepare_grid", "perform_docking", "analyse_docking", "filter_results", "cleanup_files"],
        "skipped_receptor": [],
        "skipped_ligand": [],
        "file__result_score": "result.score.csv",
        "file__config": "sieveai.cfg",
        "vina_config_keys": ["flex", # "receptor", "ligand",
              "center_x", "center_y", "center_z",
              "size_x", "size_y", "size_z",
              "out", "log",
              "cpu", "seed",
              "exhaustiveness", "num_modes", "energy_range"],

        "multiprocess": False,
      }

    self.utility.update_attributes(self, kwargs, self.__defaults)

    # To run single process as string rather a list of existing methods
    if isinstance(self.processes, str) and len(self.processes) > 3:
      self.processes = [self.processes]

    self.path_config = f"{self.path_base}/{self.file__config}"

  def __write_vina_config(self, *args, **kwargs):
    """Writes VINA config"""
    __structure = []
    __coordinates = []
    self.properties = []

    _vina_config_settings = {
      "spacing": 1,
      "residues": [],
      "energy_range": 3,
      "seed": 41103333,
      "num_modes": 10,
      "exhaustiveness": 16,
    }

    _vina_config_settings.update(kwargs)

    _receptor_path = _vina_config_settings.get("receptor")
    parser = PDBParser(get_header=False)
    __structure = parser.get_structure("PDB", _receptor_path)

    """Filter coordinates to calculate the grid of selected res only"""

    if _vina_config_settings.get("residues") and len(_vina_config_settings.get("residues")):
      # If residues are set, calculate around specified residues (residue ids)
      for chains in __structure.get_chains():
        for chain in chains:
          chain_vars = vars(chain)
          if chain_vars.get("resname") in _vina_config_settings.get("residues"):
            # print([[k for k in res.get_coord()] for res in chain])
            __coordinates.extend([[k for k in res.get_coord()] for res in chain])
    else:
      # Else draw around whole the structure
      __coordinates = []
      for atom in __structure.get_atoms():
        __coordinates.append(atom.get_coord())

    x_total, y_total, z_total = (0, 0, 0)
    x_min, x_max = (0, 0)
    y_min, y_max = (0, 0)
    z_min, z_max = (0, 0)

    d_x = 0
    d_y = 0
    d_z = 0

    center = (0, 0, 0)
    size = (0, 0, 0)
    distance = 0

    if len(__coordinates):
      for c in __coordinates:
        at_x = float(c[0])
        at_y = float(c[1])
        at_z = float(c[2])
        x_total += at_x
        y_total += at_y
        z_total += at_z
        x_min = min(x_min, at_x)
        x_max = max(x_max, at_x)
        y_min = min(y_min, at_y)
        y_max = max(y_max, at_y)
        z_min = min(z_min, at_z)
        z_max = max(z_max, at_z)

      center = [x_total/len(__coordinates), y_total/len(__coordinates), z_total/len(__coordinates)]
      size = [(x_max - x_min), (y_max - y_min), (z_max - z_min)]

      for c in __coordinates:
        d_x += (c[0]-center[0])**2
        d_y += (c[1]-center[1])**2
        d_z += (c[2]-center[2])**2

      distance = MATH.sqrt((d_x+d_y+d_z)/len(__coordinates))

    self.properties = {
      "center": center,
      "size": size,
      "distance": distance,
    }

    center_x, center_y, center_z = self.properties['center']
    center_x = round(center_x, 4)
    center_y = round(center_y, 4)
    center_z = round(center_z, 4)

    """
    VINA has maximum box size of 12.6 nm
    """
    size_x, size_y, size_z = self.properties['size']
    size_x = min(int(size_x), 126)
    size_y = min(int(size_y), 126)
    size_z = min(int(size_z), 126)

    _spacing = float(_vina_config_settings.get("spacing"))

    self.vina_config = {
      "center_x": center_x,
      "center_y": center_y,
      "center_z": center_z,
      "size_x": size_x + _spacing,
      "size_y": size_y + _spacing,
      "size_z": size_z + _spacing,
    }

    _vina_config_settings.update(self.vina_config)

    __vina_allowed_config = {key: _vina_config_settings[key] for key in _vina_config_settings if key in self.vina_config_keys}
    __vina_config_lines = [f"{config_key} = {__vina_allowed_config[config_key]}" for config_key in __vina_allowed_config.keys() if __vina_allowed_config[config_key] is not None]
    __vina_config_file_location = f"{_vina_config_settings.get('destination')}{OS.sep}{_vina_config_settings.get('config_file')}"

    if _vina_config_settings.get('config_file'):
      self.utility.write(__vina_config_file_location, __vina_config_lines)

  def __parse_res_score_line(self, line):
    records = None
    if line.startswith(" "):
      records = str(line).split()
      records = [r.strip() for r in records]

    return records

  def __parse_score_log(self, file_path):
    __v_scores = []
    with open(file_path, "r") as __vsf:
      __v_scores = __vsf.readlines()

    _score_table_dict = []
    _score_headers = ["mode", "affinity", "rmsd_lb", "rmsd_ub"]
    _score_flag = False
    for _res_line in __v_scores:
        if _score_flag:
            score_recs = self.__parse_res_score_line(_res_line)
            if score_recs:
              _score_table_dict.append(score_recs)

        if not _score_flag and _res_line.startswith("-----+------------+----------+----------"):
            _score_flag = True

    _score_df = PD.DataFrame.from_records(_score_table_dict, columns=_score_headers)
    return _score_df

  def prepare_grid(self, *args, **kwargs):
    # return # to by pass
    self.utility.log_info("Docking preparing grid...")
    _rec_pdbqt_path_list = self.utility.find_files(self.path_receptor_pdbqt, '.pdbqt')
    _rec_pre_paths = self.utility.find_files(self.path_receptor_pdbqt, ".config")

    if len(_rec_pdbqt_path_list) == len(_rec_pre_paths):
      self.utility.log_warning(f"{len(_rec_pdbqt_path_list)} receptor config files already exist. Continuing...")
      return

    for _rec in _rec_pdbqt_path_list:
      file_name = self.utility.file_name(_rec)
      params = {
        "receptor": f"{_rec}",
        "config_file": f"{file_name}.config",
        "destination": self.path_receptor_config
      }
      self.__write_vina_config(**params)

  def gen_chimerax_scripts(self, *args, **kwargs):
    _rec_pdbqt_path_list = self.utility.find_files(self.path_receptor_pdbqt, '.pdbqt')
    _lig_pdbqt_path_list = self.utility.find_files(self.path_ligand_pdbqt, '.pdbqt')

    _result_matrix = []

    _complexes = self.utility.product([_rec_pdbqt_path_list, _lig_pdbqt_path_list]) # combination
    _complexes_to_process = []

    for _rec, _lig in _complexes:
      _rec_fn = self.utility.filename(_rec)
      _lig_fn = self.utility.filename(_lig)
      _comp = _rec_fn.strip() + "--" + _lig_fn.strip()

      _res_file = f"{self.path_docking}{OS.sep}{_comp}.result.pdbqt"
      _log_file = self.utility.change_ext(_res_file, "log")
      _cxc_file = f"{self.path_analysis}{OS.sep}{_comp}.{self.ext__chimera_cxc}"

      if not self.utility.check_path(_cxc_file):
        _complexes_to_process.append((_rec_fn, _lig_fn, _comp, _rec, _lig, _res_file, _log_file, _cxc_file))

    for _rec_fn, _lig_fn, _comp, _rec, _lig, _res_file, _log_file, _cxc_file in _complexes_to_process:
        # If log file is not present, or
        if not self.utility.check_path(_log_file):
          continue

        _scores = self.__parse_score_log(_log_file)

        # Continue if out file doesn't exist
        if not self.utility.check_path(_res_file):
          self.utility.log_error(f"{_comp}.log file doesn't exist.")
          _result_matrix.append({
            "receptor": _rec_fn,
            "ligand": _lig_fn,
            "conformer_id": 0,
            "conformer_score": 0,
          })
          continue

        # ToDo: Skip chimerax execution if hbonds or contacts file

        _complex_commads = [
          f"close;"
          f"set bgColor white; open {_rec}; wait; hide surfaces; hide atoms; show cartoons; wait;"
          "addh;",
          "~sel;",
          "wait;",
          f"open {_res_file}; wait;",
        ]

        _total_conformers = len(_scores['mode'].tolist())

        for _idx, _conformer_details in _scores.iterrows():
          _model_id = _conformer_details["mode"]
          _complex_commads.extend([
            f"# MODEL-NO-{_model_id}",
            "hide #!2.1-%s target m;" % (_total_conformers),
            f"show #!2.{_model_id} models;",
            f"view;",
            f"sel #!2.{_model_id};", # Select the model
            f"contacts (#1 & ~hbonds) restrict sel radius 0.05 log t saveFile {self.path_analysis}/{_comp}--{_model_id}.{self.ext__contacts};",
            f"wait;",
            f"hb #1 restrict sel reveal t show t select t radius 0.05 log t saveFile {self.path_analysis}/{_comp}--{_model_id}.{self.ext__hbonds};",
            f"wait;",
            "label sel residues text {0.name}-{0.number} height 1.5 offset -2,0.25,0.25 bgColor #00000099 color white;",
            f"~sel;",
            # f"save {self.path_analysis}/{_comp}--{_model_id}.complex.png width 1200 height 838 supersample 4 transparentBackground true;",
            f"sel #!2.{_model_id};", # Select the model
            f"view sel;",
            f"~sel;",
            # f"save {self.path_analysis}/{_comp}--{_model_id}.ligand.png width 1200 height 838 supersample 4 transparentBackground true;",
            f"turn x 45;",
            # f"save {self.path_analysis}/{_comp}--{_model_id}.T45.ligand.png width 1200 height 838 supersample 4 transparentBackground true;",
            f"\n",
          ])

          _result_matrix.append({
            "receptor": _rec_fn,
            "ligand": _lig_fn,
            "conformer_id": _model_id,
            "conformer_score": _conformer_details["affinity"],
          })

        _complex_commads.extend([
          f"exit;",
        ])

        self.utility.write(_cxc_file, _complex_commads)
        # self.chimera_files.append(f"open {_cxc_file};")

    _score_file_path = f"{self.path_base}{OS.sep}{self.file__result_score}"
    if _result_matrix and not self.utility.check_path(_score_file_path):
      _f2 = PD.DataFrame(_result_matrix)
      _f2.to_csv(_score_file_path, index=False)
      self.utility.log_info(f"AutoDock result score table updated to {_score_file_path}.")
    else:
      self.utility.log_warning("No score file was saved .")

  def run_chimerax_scripts(self, *args, **kwargs):
    _results_score_file = f"{self.path_base}/{self.file__result_score}"

    self.utility.time_start()
    if not self.utility.check_path(_results_score_file):
      self.utility.log_warning(f"Score file {_results_score_file} does not exist to run chimerax scripts. Skipping...")
      return None

    _results = PD.read_csv(_results_score_file)
    _results["complex_id"] = _results['receptor'] + "--" + _results['ligand']
    _results = _results.drop_duplicates(subset='complex_id', keep="last")
    _total_to_process = _results.shape[0]
    self.utility.log_info(f"{_total_to_process} records to be processed.")

    _processed_count = 0

    for _idx, _complex in self.utility.ProgressBar(_results.iterrows()):
      _processed_count = _processed_count + 1
      _rec = _complex['receptor']
      _lig = _complex['ligand']
      _last_conformer = _complex['conformer_id']
      _complex_name = f"{_rec}--{_lig}"
      _cxc_file = f"{self.path_analysis}{OS.sep}{_complex_name}.cxc"

      # Skip if hbonds.txt file is already generated
      self.utility.log_info(f"Processing {_complex_name} complex with ChimeraX %i/%i having {_last_conformer} conformers." % (_processed_count, _total_to_process))
      if self.utility.check_path(_cxc_file):
        if self.utility.check_path(f"{self.path_analysis}{OS.sep}{_complex_name}--{_last_conformer}.{self.ext__hbonds}"):
          self.utility.log_info(f"{_complex_name}--{_last_conformer} is already processed by ChimeraX.")
          continue
        # @TODO: Multiprocess
        self.chimerax_run_file(_cxc_file)
      else:
        self.utility.log_info(f"File {_cxc_file} not found. Continuing...")

    self.utility.log_info(f"ChimeraX scripts execution is completed. Took {self.utility.time_elapsed()}.")

  def process_chimerax_results(self, *args, **kwargs):
    self.utility.time_start()
    self.utility.log_info("Starting processing ChimeraX results.")
    _results_score_file = f"{self.path_base}{OS.sep}{self.file__result_score}"

    if not self.utility.check_path(_results_score_file):
      self.utility.log_warning(f"Score file {_results_score_file} does not exist for processing chimerax results. Skipping...")
      return None

    _results = PD.read_csv(_results_score_file)
    _result_matrix = []

    for _idx, _complex in self.utility.ProgressBar(_results.iterrows()):
      _rec = _complex['receptor']
      _lig = _complex['ligand']
      __conf_id = _complex['conformer_id']
      __conf_file = f"{_rec}--{_lig}--{__conf_id}"
      __contacts = self.chimerax_get_contacts_residues(f"{self.path_analysis}{OS.sep}{__conf_file}.{self.ext__contacts}")
      __hbonds = self.chimerax_get_hbonds_residues(f"{self.path_analysis}{OS.sep}{__conf_file}.{self.ext__hbonds}")

      __contacts_count, __hbonds_count = 0, 0
      if isinstance(__contacts, PD.DataFrame):
        __contacts["residues"] =  __contacts["atom1__resname"].astype(str) + ":" + __contacts["atom1__resid"].astype(str)
        __contacts = __contacts["residues"].tolist()
        __contacts_count = len(__contacts)
        __contacts = ",".join(__contacts)
      else:
        __contacts = ""

      if isinstance(__hbonds, PD.DataFrame) :
        __hbonds["residues"] =  __hbonds["donor__resname"].astype(str) + ":" + __hbonds["donor__resid"].astype(str)
        __hbonds = __hbonds["residues"].tolist()
        __hbonds_count = len(__hbonds)
        __hbonds = ",".join(__hbonds)
      else:
        __hbonds = ""

      _result_matrix.append({
          "receptor": _rec,
          "ligand": _lig,
          "conformer_id": _complex['conformer_id'],
          "conformer_score": _complex['conformer_score'],
          "contacts_count": __contacts_count,
          "hbonds_count": __hbonds_count,
          "contacts": __contacts,
          "hbonds": __hbonds,
        })

      self.utility.log_info(f"Processed {_rec}--{_lig}--{_complex['conformer_id']}.")

    _f2 = PD.DataFrame(_result_matrix)
    _f2.to_csv(_results_score_file, index=False)
    self.utility.log_info(f"Processed ChimeraX results and updated to {_results_score_file}.")

  def cleanup_files(self, *args, **kwargs):
    self.utility.update_attributes(**kwargs)

    # Check if tgz file for the directory name is present then skip the process
    self.utility.log_info("Cleaning up the files.")

    self.path_analysis = f"{self.path_base}/{self.dir_analysis}"
    self.path_analysis_tgz = f"{self.path_analysis}.tgz"

    _analysis_files = self.utility.find_files(self.path_analysis, [self.ext__contacts, self.ext__hbonds, self.ext__chimera_cxc])

    if len(_analysis_files) > 0:
      self.utility.log_info(f"{len(_analysis_files)} files to be compressed.")
      self.utility.add_tgz_files(self.path_analysis_tgz, _analysis_files)

    if self.utility.check_path(self.path_analysis_tgz):
      self.utility.read_tgz_file(self.path_analysis_tgz)
      _archived_files = self.utility.tgz_files

      self.utility.log_info(f"{len(_analysis_files)} files in {self.path_analysis_tgz} to be deleted.")
      for _l in self.utility.ProgressBar(_analysis_files):
        _fname = self.utility.file_name(_l, with_ext=True)
        if _fname in _archived_files:
          # Delete the files
          self.utility.delete_path(_l)

      # As directory is empty now, delete the directory
      self.utility.delete_path(self.path_analysis)
    else:
      self.utility.log_error(f"{self.path_analysis_tgz} file is not present.")

    # ext__dock_log
    self.path_docking = f"{self.path_base}/{self.dir_docking}"
    self.path_docking_log_tgz = f"{self.path_docking}.log.tgz"
    _dock_files = self.utility.find_files(self.path_docking, [self.ext__dock_log])
    if len(_dock_files) > 0:
      self.utility.log_info(f"{len(_dock_files)} files in {self.path_docking_log_tgz} to be compressed.")
      self.utility.add_tgz_files(self.path_docking_log_tgz, _dock_files)

    if self.utility.check_path(self.path_docking_log_tgz):
      self.utility.read_tgz_file(self.path_docking_log_tgz)
      _arch_dock_files = self.utility.tgz_files
      for _l in self.utility.ProgressBar(_dock_files):
        _fname = self.utility.file_name(_l, with_ext=True)
        if _fname in _arch_dock_files:
          # Delete the log files only
          self.utility.delete_path(_l)

      # ToDo: Update using Project Manager
    else:
      self.utility.log_error(f"{self.path_docking_log_tgz} file is not present.")

  def analyse_docking(self, *args, **kwargs):
    # return # to by pass
    if self.path_analysis is None and self.dir_analysis:
      self.path_analysis = self.utility.validate_dir(f"{self.path_base}/{self.dir_analysis}")

    self.gen_chimerax_scripts(**kwargs)
    self.run_chimerax_scripts(**kwargs)
    self.process_chimerax_results(**kwargs)

  """
  Filter results based on DL Classifiers
  """
  def filter_results(self, *args, **kwargs):
    ...

  def dock(self, *args, **kwargs):
    _command = args[0] if len(args) > 0 else kwargs.get("command")
    _log_file = args[1] if len(args) > 1 else kwargs.get("log_file")

    _output = self.utility.cmd_run(_command)
    try:
      # Write single log
      self.utility.write(_log_file, str(_output))
    except:
      self.utility.log_info(f"Problem occurred in processing command {_command}. File '{_log_file}' not generated.")

    return True

  def perform_docking(self, *args, **kwargs):
    self.utility.time_start()

    # _docked_results = []

    if self.path_docking is None and self.dir_docking:
      self.path_docking = self.utility.validate_dir(f"{self.path_base}/{self.dir_docking}")

    # @TODO: Filter rec for skipped list
    _rec_pdbqt_path_list = self.utility.find_files(self.path_receptor_pdbqt, '.pdbqt')

    # @TODO: Filter lig for skipped list
    _lig_pdbqt_path_list = self.utility.find_files(self.path_ligand_pdbqt, '.pdbqt')

    _complexes = self.utility.product([_rec_pdbqt_path_list, _lig_pdbqt_path_list]) # combination
    _complexes_to_process = []

    total_complexes = 0
    for _rec, _lig in _complexes:
      total_complexes += 1
      _rec_fn = self.utility.filename(_rec)
      _lig_fn = self.utility.filename(_lig)
      _comp = _rec_fn.strip() + "--" + _lig_fn.strip()
      _res_file = f"{self.path_docking}{OS.sep}{_comp}.{self.ext__dock_result}"
      _log_file = f"{self.path_docking}{OS.sep}{_comp}.{self.ext__dock_log}"

      if not all((self.utility.check_path(_res_file), self.utility.check_path(_log_file))):
        _complexes_to_process.append((_rec, _lig, _comp, _res_file))

    if self.multiprocess:
      self.utility.log_info(f"Adding remaining {len(_complexes_to_process)}/{total_complexes} complex(s) to multiprocess queue.")
    else:
      self.utility.log_info(f"Processing remaining {len(_complexes_to_process)}/{total_complexes} complex(s).")

    for _rec, _lig, _comp, _res_file in self.utility.ProgressBar(_complexes_to_process):
      _log_file = self.utility.change_ext(_res_file, "log")
      _config_file = f"{self.path_receptor_config}{OS.sep}{self.utility.filename(_rec)}.config"
      _rec_pdbqt = f"{self.path_receptor_pdbqt}{OS.sep}{self.utility.filename(_rec)}.pdbqt"
      if not self.utility.check_path(_config_file):
        self.utility.log_info(f"Config file '{_config_file}' for the receptor '{self.utility.filename(_rec)}' is missing.")
        continue

      _command = ['vina',
            "--cpu", "4",
            "--receptor", _rec_pdbqt,
            "--config", _config_file,
            "--ligand", _lig,
            "--exhaustiveness", "32",
            "--out", _res_file,
            "--verbosity", "2",
            # "--log", f"{self.path_docking}{OS.sep}{result_file}.result.log"
          ]
      if self.multiprocess:
        _process = self.utility.multiprocess_add(self.dock, args=(_command, _log_file))
      else:
        self.dock(_command, _log_file)

    if self.multiprocess:
      self.utility.multiprocess_start()

    self.utility.log_info(f"Completed the docking process. Took {self.utility.time_elapsed()}.")

    del _rec_pdbqt_path_list
    del _lig_pdbqt_path_list
    del _complexes
    del _complexes_to_process

  def prepare_dynamics(self, *args, **kwargs):
    pass

  def add(self, *args, **kwargs):
    self.utility.update_attributes(self, kwargs)

  def run(self, *args, **kwargs):
    """Run all options from the given settings."""
    _start_time = self.utility.time_get()
    self.utility.log_info("Starting %s." % " -> ".join(self.processes))
    _result = {}
    for _process in self.processes:
      _method = getattr(self, _process, None)
      if hasattr(_method, '__call__'):
        _result[_process] = _method(**kwargs)
      else:
        self.utility.log_error(f"{_method} is not callable.")
    self.utility.log_info(f"Completed all the processes.")
    return _result
