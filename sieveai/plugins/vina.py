from __future__ import annotations

from ..sieveaibase import StepManager, DictConfig
from ..managers import Structures
from ..process.docking import PluginDockingBase

from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings as WARNINGS

"""
ToDo:
  Clean PDB files
  Fix alternate occupancy along with residue name for correct PDBQT coordinates
  * https://www.cgl.ucsf.edu/chimerax/docs/user/commands/altlocs.html
  * alt list; alt clean;

"""

class Vina(PluginDockingBase):
  is_ready = False
  plugin_name = "AutoDock VINA"
  plugin_uid = "AutoDockVINA"
  plugin_version = "1.0"
  assignments = ['docking']
  current_assignment = None
  dependencies = ['chimerax', 'cir']
  url = "https://autodock-vina.readthedocs.io/en/latest/docking_python.html"

  Vina_config = None

  _num_modes = 5

  path_vina_exe = None

  def __init__(self, *args, **kwargs):
    WARNINGS.simplefilter('ignore', PDBConstructionWarning)
    super().__init__(**kwargs)

  def _init_molecules(self, *args, **kwargs) -> None:
    self.TASKS.start_step('init_molecules', self.plugin_uid, self.plugin_uid)
    self.update_attributes(self, kwargs)

    self.path_excel_results = self._get_result_excel_path()

    self.path_vina_exe = 'vina' # self.which('vina')

    self.require('re', 'RegEx')
    self.re_contacts = self.RegEx.compile(f"(\d+) contacts")

    self._restore_progress()

    self._set_defualt_config()

    self._steps_map_methods = {
      "init": self._prepare_molecules,
      "config": self._write_receptor_vina_config,
      "dock": self._run_cli_docking,
      # "dock": self._run_api_docking,
      "analyse": self._parse_analyse_interactions,
      "final": self._finalise_complex,
    }

    self._step_sequence = tuple(self._steps_map_methods.keys())

    if True: # Compare hash of file changes and reattach selectively
      self.Receptors = Structures(self.SETTINGS.user.path_receptors, ['protein', 'dna', 'rna'])
      self.Ligands = Structures(self.SETTINGS.user.path_ligands, ['compound', 'protein', 'rna'])
      self.log_debug('Storing Molecules')

    self.TASKS.end_step('init_molecules', self.plugin_uid, self.plugin_uid)

  def _set_defualt_config(self):

    self.Vina_config = DictConfig()

    self.Vina_config.default = {
        "exhaustiveness": 16,
        "verbosity": 2,
        "center_x": None,
        "center_y": None,
        "center_z": None,
        "size_x": None,
        "size_y": None,
        "size_z": None,
        "out": None,
        "cpu": None,
        "log": None,
        "seed": 41103333,
        "num_modes": 10,
      }

    self.Vina_config.other.spacing = 1
    self.Vina_config.other.residues = []

    self.Vina_config.energy_range = 3

    self.Vina_config.allowed_keys = ["flex", "receptor", "ligand",
                  "center_x", "center_y", "center_z",
                  "size_x", "size_y", "size_z",
                  "out", "log",
                  "cpu", "verbosity", # "seed",
                  "exhaustiveness", "num_modes", "energy_range"]

  def _run_api_docking(self, cuid):
    """
    @ToDo: Docking using python binding for AutoDock VINA
    """

    from vina import Vina as VinaPy
    _vna = VinaPy(sf_name='vina')
    _cuid_c = self.Complexes[cuid]
    _vna.set_receptor(str(_cuid_c.path_receptor))

    _vna.set_ligand_from_file(str(_cuid_c.path_ligand))
    _vna.compute_vina_maps(center=_cuid_c.center, box_size=_cuid_c.box_size)

    # Score the current pose
    energy = _vna.score()
    print('Score before minimization: %.3f (kcal/mol)' % energy[0])

    # Minimized locally the current pose
    energy_minimized = _vna.optimize()
    print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
    _vna.write_pose(str(_cuid_c.path_ligand).with_suffix('.min.pdbqt'), overwrite=True)

    # Dock the ligand
    _vna.dock(exhaustiveness=32, n_poses=20)
    _vna.write_poses(str(_cuid_c.path_out), n_poses=5, overwrite=True)

  def _run_cli_docking(self, cuid) -> None:
    """Runs vina command."""
    self.TASKS.start_step('run_cli_docking', cuid, self.plugin_uid)
    _cuid_c = self.Complexes[cuid]

    if _cuid_c.path_out.exists():
      self.log_debug(f'Vina result for {cuid} exists. Returning...')
      self.TASKS.end_step('run_cli_docking', cuid, self.plugin_uid)
      return

    _config = {
            # "--cpu": 1,
            "--receptor": _cuid_c.path_receptor.resolve(),
            "--ligand": _cuid_c.path_ligand.resolve(),
            "--config": _cuid_c.path_vina_config.resolve(),
            "--out": _cuid_c.path_out.resolve(),
            "--num_modes": self._num_modes,
            # "--verbosity": self.Vina_config.default.get('verbosity', 2),
            "cwd": _cuid_c.path_docking.resolve(),
            # ">": _cuid_c.path_score.resolve() # from command line to file
          }
    _result = self.cmd_run(str(self.path_vina_exe), **_config)
    _cuid_c.path_score.write_text(_result)
    self.TASKS.end_step('run_cli_docking', cuid, self.plugin_uid)

  _score_headers = ["mode", "affinity", "rmsd_lb", "rmsd_ub"]

  def _parse_analyse_interactions(self, cuid):
    self.TASKS.start_step('parse_analyse_interactions', cuid, self.plugin_uid)
    _cuid_c = self.Complexes[cuid]

    _cuid_c.path_cxc_cmd = (_cuid_c.path_docking / 'analysis.cxc').rel_path()

    if not _cuid_c.path_score.exists():
      self.TASKS.end_step('parse_analyse_interactions', cuid, self.plugin_uid)
      return

    _score = list(_cuid_c.path_score.readlines())
    _score_flag = False
    _score_records = []
    for _res_line in _score:
      if _score_flag and _res_line.startswith(" "):
          _record = str(_res_line).split()
          _record = [r.strip() for r in _record]
          _score_records.append(_record)

      if not _score_flag and _res_line.startswith("-----+------------+----------+----------"):
        _score_flag = True

    _df_score = self.DF(_score_records, columns=self._score_headers)

    _models = _df_score['mode'].tolist()
    _CX =  self.SETTINGS.PLUGIN_REFS.chimerax()
    _file_template_cxc_contacts = 'CXC-Result-Model-%s.contacts.txt'
    _file_template_cxc_hbonds = 'CXC-Result-Model-%s.hbonds.txt'

    _n_models_contacts = len([*_cuid_c.path_docking.search(_file_template_cxc_contacts % '*')])
    _n_models_hb = len([*_cuid_c.path_docking.search(_file_template_cxc_hbonds % '*')])

    if (_n_models_contacts +  _n_models_hb) > 1:
      self.log_debug(f'ChimeraX analysis results (n={_n_models_contacts +  _n_models_hb}) for {cuid} exists, skipping CXC analysis...')
    else:
      _complex_commads = [
        f"close;"
        f"set bgColor white; open {_cuid_c.path_receptor.resolve()}; wait; hide surfaces; hide atoms; show cartoons; wait;"
        "addh;",
        "~sel;",
        "wait;",
        f"open {_cuid_c.path_out.resolve()}; wait;",
      ]

      for _model_id in _models:
        _path_contacts = (_cuid_c.path_docking / (_file_template_cxc_contacts % _model_id)).resolve()

        if _path_contacts.exists(): continue

        _complex_commads.extend([
          f"# MODEL-NO-{_model_id}",
          "hide #!2.1-%s target m;" % (len(_models)),
          f"show #!2.{_model_id} models;",
          f"view;",
          f"sel #!2.{_model_id};", # Select the model
          f"contacts (#1 & ~hbonds) restrict sel radius 0.05 log t saveFile {_path_contacts};",
          f"wait;",
          f"hb #1 restrict sel reveal t show t select t radius 0.05 log t saveFile {(_cuid_c.path_docking / (_file_template_cxc_hbonds % _model_id)).resolve()};",
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

      _complex_commads.extend([
          f"exit;",
        ])

      _cuid_c.path_cxc_cmd.write_text("\n".join(_complex_commads))
      _CX.exe_cxc_file(_cuid_c.path_cxc_cmd.resolve())

    _summary = []
    for _model_id in _models:
      _contacts_df = _CX.parse_contacts((_cuid_c.path_docking / (_file_template_cxc_contacts % _model_id)).resolve())
      _hbonds_df = _CX.parse_hbonds((_cuid_c.path_docking / (_file_template_cxc_hbonds % _model_id)).resolve())

      self.Complexes[cuid][f'Model_{_model_id}'].contacts = _contacts_df
      self.Complexes[cuid][f'Model_{_model_id}'].hbonds = _hbonds_df

      _contacts, _hbonds = [], []
      _NoneType = type(None)

      if not isinstance(_contacts_df, _NoneType):
        _contacts = _contacts_df.copy()
        _contacts["residues"] =  _contacts["atom1__resname"].astype(str) + ":" + _contacts["atom1__resid"].astype(str)
        """Return only unique interacting residues"""
        _contacts = set(_contacts["residues"].tolist())

      _len_contacts = len(_contacts)
      _contacts = ",".join(_contacts)

      if not isinstance(_hbonds_df, _NoneType):
        _hbonds = _hbonds_df.copy()
        _hbonds["residues"] =  _hbonds["donor__resname"].astype(str) + ":" + _hbonds["donor__resid"].astype(str)
        """Return only unique interacting residues"""
        _hbonds = set(_hbonds["residues"].tolist())

      _len_hbonds = len(_hbonds)
      _hbonds = ",".join(_hbonds)

      _summary.append({
        'mode': _model_id,
        'total_contacts': _len_contacts,
        'total_hbonds': _len_hbonds,
        'contacts': _contacts,
        'hbonds': _hbonds,
      })

    if len(_summary) > 0:
      _df_summary = self.DF(_summary)
      _df_score = self.PD.merge(_df_score, _df_summary, on='mode', how='outer')
      # _df_score.columns = ['Conformer ID', 'VINA Score', 'RMSD LB', 'RMSD UB', 'Total Contacts', 'Total Hbonds', 'Contacts', 'Hbonds']
      self.Complexes[cuid].vina_results = _df_score

    self.TASKS.end_step('parse_analyse_interactions', cuid, self.plugin_uid)

  def _prepare_molecules(self, cuid):
    return cuid

  def _finalise_complex(self, cuid):
    return cuid

  def _process_complex(self, cuid) -> None:
    self.TASKS.start_step('process_complex', cuid, self.plugin_uid)
    if not cuid in self.Complexes:
      self.TASKS.end_step('process_complex', cuid, self.plugin_uid)
      return cuid

    for _step in self.Complexes[cuid].step:
      if _step in self.Complexes[cuid].steps_completed:
        continue

      _step = self.Complexes[cuid].step._current

      if not self.Complexes[cuid].step.is_last:
        self.log_debug(f"{cuid}:: Step: {_step}")
      else:
        self.log_debug(f"{cuid}:: Last Step: {_step}")

      self._steps_map_methods[_step](cuid)
      self.Complexes[cuid].steps_completed.append(_step)

    self.TASKS.end_step('process_complex', cuid, self.plugin_uid)
    return cuid

  def _queue_complexes(self) -> None:
    self.TASKS.start_step('queue_complexes', self.plugin_uid, self.plugin_uid)

    self.log_debug(f"Receptors: {self.Receptors}")
    self.log_debug(f"Ligands: {self.Ligands}")

    if self.path_vina_exe is None:
      self.log_error(f'{self.path_vina_exe} executable is not found. Please provide the executable command or path to executable file.')
      return

    # Perform Docking
    self.SETTINGS.user.multiprocessing and self.init_multiprocessing()
    _combs = list(self.product(self.Receptors.keys(), self.Ligands.keys()))

    for _rck, _ligk in _combs:
      _rec = self.Receptors[_rck]
      _lig = self.Ligands[_ligk]

      _complex_uid = self.slug(f"{_rec.mol_id}--{_lig.mol_id}")

      if not _complex_uid in self.Complexes:
        self.Complexes[_complex_uid] = DictConfig()

        _complex_path = (self.path_plugin_docking / _complex_uid).validate()

        # Convert mol_path to PDBQT using meeko/OpenBabel???
        _c_rec = _complex_path / f'REC.pdbqt'
        _c_lig = _complex_path / f'LIG.pdbqt'

        if not _rec.formats['.pdbqt'].mol_path.exists():
          self.log_debug(f'{_rec.mol_id} PDBQT does not exist.')
          continue

        if not _lig.formats['.pdbqt'].mol_path.exists():
          self.log_debug(f'{_lig.mol_id} PDBQT does not exist.')
          continue

        _rec.formats['.pdbqt'].mol_path.copy(_c_rec.resolve())
        _lig.formats['.pdbqt'].mol_path.copy(_c_lig.resolve())

        self.Complexes[_complex_uid].update({
            "step": StepManager(self._step_sequence),
            "steps_completed": [],
            "uid": _complex_uid,
            "rec_uid": _rec.mol_id,
            "lig_uid": _lig.mol_id,
            "path_receptor": _c_rec,
            "path_ligand": _c_lig,
            "path_docking": _complex_path,
            "path_out": _complex_path / f'{_complex_uid}.out.pdbqt',
            "path_score": _complex_path / f'{_complex_uid}.vina.txt',
            "path_vina_config": _complex_path / f"{_complex_uid}.vina.config"
          })

        self.log_debug(f'VINA_01: {_complex_uid} initiated.')

      if self.SETTINGS.user.multiprocessing:
        self.queue_task(self._process_complex, _complex_uid)
        self.log_debug(f'VINA_02: {_complex_uid} queued for multiprocessing.')
      else:
        self.log_debug(f'VINA_03: {_complex_uid} is being processed...')
        self._process_complex(_complex_uid)

    self.TASKS.end_step('queue_complexes', self.plugin_uid, self.plugin_uid)

  def _rank_conformers(self, _df_all_conformers):
    self.TASKS.start_step('rank_conformers', self.plugin_uid, self.plugin_uid)
    _ranking_columns = {
        'affinity': True,
        'total_contacts': False,
        'total_hbonds': False,
        # 'rmsd': True # for site specific docking
      }

    _all_result_order = ['Composite_Rank', 'affinity', 'mode', 'rec_uid', 'lig_uid', 'total_contacts', 'total_hbonds', 'contacts', 'hbonds']
    _top_result_order = ['Composite_Rank', 'affinity', 'mode', 'rec_uid', 'lig_uid', 'total_contacts', 'total_hbonds', 'contacts', 'hbonds']

    _column_rename_map = {
      'conformer_uid': 'Complex UID',
      'complex_uid': 'Complex ID',
      'rec_uid': 'Receptor ID',
      'lig_uid': 'Ligand ID',

      'mode': 'Conformer ID',
      'affinity': 'VINA Score',
      'rmsd_lb': 'RMSD LB',
      'rmsd_ub': 'RMSD UB',
      'rmsd_ub': 'RMSD UB',

      'total_contacts': 'Total Contacts',
      'total_hbonds': 'Total Hbonds',
      'contacts': 'Contact Residues',
      'hbonds': 'HBond Residues',
      'contacts_count': 'Total Contacts', # Deprecated
      'hbonds_count': 'Total Hbonds', # Deprecated

      'Rank': 'Composite Rank',  # Deprecated
      'Composite_Rank': 'Composite Rank',
      'Composite_Score': 'Composite Score',
    }

    _ranker = self.SETTINGS.PLUGIN_REFS.CompositeIndexRanker() # CIR
    _df_all_conformers, _top_res = _ranker.rank_conformers(
          _df_all_conformers, _ranking_columns,
            grouper = 'lig_uid'
        )

    _df_all_conformers = _df_all_conformers[_all_result_order]
    _top_res = _top_res[_top_result_order]

    _df_all_conformers.rename(columns=_column_rename_map, inplace=True)
    _top_res.rename(columns=_column_rename_map, inplace=True)

    self.TASKS.end_step('rank_conformers', self.plugin_uid, self.plugin_uid)
    return _df_all_conformers, _top_res

  _df_all_results = None
  _df_top_ranked = None
  def _tabulate_results(self, *args, **kwargs):
    # Combine all vina_results, concatenate them, and rank conformers
    self.require('pandas', 'PD')

    # Combine all the interactions
    _score_tables = []
    for _idx in self.Complexes._keys:
      _cmplx = self.Complexes[_idx]
      if not all([isinstance(_cmplx, (dict)), 'vina_results' in _cmplx]):
        continue

      _s = _cmplx.vina_results
      _s['rec_uid'] = _cmplx.rec_uid
      _s['lig_uid'] = _cmplx.lig_uid
      _s['complex_uid'] = _cmplx.uid
      _score_tables.append(_s)

    if not len(_score_tables) > 0:
      self.log_debug('No results were found to be concatenated.')
      return

    _score_tables = self.PD.concat(_score_tables)

    # Save conformers and ranks as excel
    self._df_all_results, self._df_top_ranked = self._rank_conformers(_score_tables)

    # self.pd_excel(self.path_excel_results, _score_tables, sheet_name=f"Raw-Results")
    self.pd_excel(self.path_excel_results, self._df_all_results, sheet_name=f"{self.plugin_uid}-All-Ranked")
    self.pd_excel(self.path_excel_results, self._df_top_ranked, sheet_name=f"{self.plugin_uid}-Top-Ranked")

    return self._df_all_results, self._df_top_ranked

  _df_all_results = None
  def _cxc_generate_images(self, *args, **kwargs):
    # Generate CXC script to generate bulk interaction images
    # Generate HTML viewer file (All/Filtered?) to manually update images
    _top_ranked = kwargs.get('top_ranked_df', args[0] if len(args) > 0 else self._df_top_ranked)
    if None is _top_ranked:
      # Read sheet from excel

      # If no top ranked data
      self.log_debug('No data was provided for image generation')

    _CX =  self.SETTINGS.PLUGIN_REFS.chimerax()

    _cx_cxc_data = ""
    _cx_html_data = _CX.get_html_head()

    _cmd_contacts = _CX.get_cmd_contacts()
    _cmd_hbonds = _CX.get_cmd_hbonds()

    # self.write(_cx_html_path, _cx_html_data)
    # self.write(_cx_cxc_path, _cx_cxc_data)

    # """Run on ChimeraX silently"""
    # _CX.exe_cxc_file(_cx_cxc_path.resolve())

    # """Open in ChimeraX"""
    # self.OS.system(f'chimerax --cmd "open {_cx_html_path}" &')

  def _finalise_results(self, *args, **kwargs):
    self.TASKS.start_step('Finalise_Results', 'FINAL_STEP', self.plugin_uid)
    self._tabulate_results()
    self._cxc_generate_images()
    self.TASKS.end_step('Finalise_Results', 'FINAL_STEP', self.plugin_uid)

  def _write_receptor_vina_config(self, cuid):
    self.TASKS.start_step('write_receptor_vina_config', cuid, self.plugin_uid)
    if not cuid in self.Complexes:
      return

    _complx = self.Complexes[cuid]
    _check = ['path_vina_config' in _complx, _complx.path_vina_config.exists()]

    if all(_check) and _complx.path_vina_config.size > 0:
      self.log_debug(f'{cuid} vina config file already exists. Skipping...')
      self.TASKS.end_step('write_receptor_vina_config', cuid, self.plugin_uid)
      return

    _config_path = self.Complexes[cuid].path_vina_config

    _mol_obj = self.Receptors[_complx.rec_uid]

    _Parser = PDBParser(get_header=False)
    _structure = _Parser.get_structure(_mol_obj.mol_id, _mol_obj.mol_path)

    _coordinates = []

    # If config settings has specific residues for site specific docking then prepare grid around specific residues
    if self.Vina_config.other.get("residues") and len(self.Vina_config.other.get("residues")):
      for _chains in _structure.get_chains():
        for _chain in _chains:
          _chain_vars = vars(_chain)
          if _chain_vars.get("resname") in self.Vina_config.other.get("residues"):
            _coordinates.extend([[k for k in res.get_coord()] for res in _chain])
    else:
      _coordinates = [_atom.get_coord() for _atom in _structure.get_atoms()]

    # Calculate center, size, and distance
    _coord_x, _coord_y, _coord_z = zip(*_coordinates) if _coordinates else ([], [], [])

    _center = (sum(_coord_x) / len(_coord_x), sum(_coord_y) / len(_coord_y), sum(_coord_z) / len(_coord_z))

    _size = (max(_coord_x) - min(_coord_x), max(_coord_y) - min(_coord_y), max(_coord_z) - min(_coord_z))

    # Prepare VINA config
    _center_x, _center_y, _center_z = (round(_coord, 4) for _coord in _center)
    _size_x, _size_y, _size_z = (min(int(dim), 126) for dim in _size)

    _spacing = float(self.Vina_config.other.get("spacing", 1))

    _vina_config = {
        **self.Vina_config.default,
        "center_x": _center_x,
        "center_y": _center_y,
        "center_z": _center_z,
        "size_x": _size_x + _spacing,
        "size_y": _size_y + _spacing,
        "size_z": _size_z + _spacing,
      }

    self.Complexes[cuid].center = [_center_x, _center_y, _center_z]
    self.Complexes[cuid].box_size = [_size_x + _spacing, _size_y + _spacing, _size_z + _spacing]

    _vina_config = {_k: _vina_config[_k] for _k in _vina_config if _k in self.Vina_config.allowed_keys}

    _vina_config_lines = [f"{config_key} = {_vina_config[config_key]}" for config_key in _vina_config.keys() if _vina_config[config_key] is not None]

    _vina_config_lines = "\n".join(_vina_config_lines + ["\n"])

    _config_path.write(_vina_config_lines)

    self.TASKS.end_step('write_receptor_vina_config', cuid, self.plugin_uid)

  def _start_preparation(self, *args, **kwargs):
    # Setting molecular formats
    self.TASKS.start_step('start_preparation', self.plugin_uid, self.plugin_uid)
    _mgltools = self.SETTINGS.PLUGIN_REFS.MGLTools()
    self.Receptors.set_format('pdbqt', converter=_mgltools.convert_pqbqt) # convert_pqbqt, prepare_receptor

    _openBabel = self.SETTINGS.PLUGIN_REFS.OpenBabel()
    self.Ligands.set_format('pdbqt', converter=_mgltools.convert_pqbqt) # _openBabel.convert

    self._queue_complexes()

    if self.SETTINGS.user.multiprocessing:
      self.process_queue()
      self.queue_final_callback(self._finalise_results)
    else:
      self._finalise_results()

    self.TASKS.end_step('start_preparation', self.plugin_uid, self.plugin_uid)

  def boot(self, *args, **kwargs) -> Vina:
    # Check configuration from SETTINGS
    self.TASKS.start_step('boot', self.plugin_uid, self.plugin_uid)
    if not self.plugin_uid in self.SETTINGS.plugin_data:
      self.SETTINGS.plugin_data[self.plugin_uid] = DictConfig()

    self._init_molecules(**kwargs)
    self._update_progress()
    self.TASKS.end_step('boot', self.plugin_uid, self.plugin_uid)
    return self

  def run(self, *args, **kwargs) -> Vina:
    self.TASKS.start_step('run', self.plugin_uid, self.plugin_uid)
    self._start_preparation(**kwargs)
    self._update_progress()
    self.TASKS.end_step('run', self.plugin_uid, self.plugin_uid)
    return self

  def shutdown(self, *args, **kwargs) -> Vina:
    self.TASKS.start_step('shutdown', self.plugin_uid, self.plugin_uid)
    self._update_progress()
    self.post_docking()
    self.TASKS.end_step('shutdown', self.plugin_uid, self.plugin_uid)
    return self
