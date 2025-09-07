from ..process.docking import PluginDockingBase
from ..sieveaibase import EntityPath

class FlexX(PluginDockingBase):
  """
  Arranges protein pdbs (like downloading and defining binding pockets in mol2 format)
  Aggregates or rearranges ligands (Chunks of 1000 into single SDF for docking with FlexX)
  After docking with LeadIT-FlexX, analysis is performed.
  """
  is_ready = False
  plugin_name = "FlexX"
  plugin_uid = "FlexX"
  plugin_version = "1.0"
  assignments = ['docking']
  current_assignment = None
  dependencies = ['chimerax', 'cir']
  url = "https://www.biosolveit.de/products/#FlexX"

  def __init__(self, *args, **kwargs):
    super().__init__(**kwargs)
    self.SETTINGS.user.path_base = self.path_base

  def setup(self, *args, **kwargs):
    self.update_attributes(self, kwargs)

  def rank_flexx_results(self, *args, **kwargs):
    _ranking_columns = {
            'Score': True,
            'Clash': False,
            'Match': True,
            'Lipo': True,
            'Ambig': True,
          }

    _ranker = self.SETTINGS.PLUGIN_REFS.CompositeIndexRanker() # CIR
  def _split_sdf_file(self, *args, **kwargs):
    _path_sdf_file = kwargs.get('path_sdf_file', args[0] if len(args) > 0 else None)

    _sdf_dir = _path_sdf_file.parent() / _path_sdf_file.stem
    if _sdf_dir.exists():
      _sdf_dir = _sdf_dir + '-SDFs'

    _path_poses_dir = kwargs.get('path_poses_dir', args[1] if len(args) > 1 else _sdf_dir.validate())
    _sdf_model_delimiter = kwargs.get('sdf_model_delimiter', args[2] if len(args) > 2 else "$$$$")

    _file_sdf_uid = None
    _last_line = ""

    _model_n_content = []
    for _line in _path_sdf_file.readlines():
      if _file_sdf_uid is None or _last_line.startswith(_sdf_model_delimiter):
        _file_sdf_uid = _line.strip().strip("\n")

      if _line.startswith(_sdf_model_delimiter):
        _new_sdf_path = _path_poses_dir / f"{self.slug(_file_sdf_uid)}.sdf"
        if not _new_sdf_path.exists():
          _model_n_content[0] = _file_sdf_uid + "\n"
          _new_sdf_path.write("".join(_model_n_content))
          _model_n_content = []
      else:
          _model_n_content.append(_line)

      _last_line = _line

    return True

  def combine_flexx_csv(self, *args, **kwargs):
    _path_csv_files = kwargs.get('path_csv_files', args[0] if len(args) > 0 else [])
    _flexx_xlsx_path = self._get_result_excel_path()

    _combined_res = None
    for _csv in _path_csv_files:
      _csv_df = self.read_csv(_csv, sep=';')
      _csv_df['File_Name'] = _csv.stem
      _combined_res = self.PD.concat([_combined_res, _csv_df]) if not _combined_res is None else _csv_df

    _flexx_comp_extract = self.re_compile(r'\((\d+)\)(.*)_(.*)')
    _combined_res[['ACNP', 'Compound_ID', 'Pose_Number']] = _combined_res.apply(lambda _x: _flexx_comp_extract.findall(str(_x['Posename']))[0], result_type='expand', axis=1)
    _combined_res = _combined_res.reset_index()
    _combined_res['Result_ID'] = _combined_res.index + 1
    _combined_res['Compound_ID'] = _combined_res.Compound_ID.apply(lambda _x: _x.strip())

    self.to_excel(_flexx_xlsx_path, _combined_res, sheet_name='Raw_Results')

  def run_flexx_plugin(self, *args, **kwargs):
    _fxp = self.path_plugin_docking # <<<<<<

    _fxx_files = [*_fxp.search('*.fxx')] + [*_fxp.search('*/*.fxx')]

    # FXX wise results or just search for CSV and SDF?

    for _file in _fxx_files:
      _sdf_files = [*_file.parent().search(f'{_file.stem}*.sdf')]
      _csv_files = [*_file.parent().search(f'{_file.stem}*.csv')]
      self.combine_flexx_csv(_csv_files)
      for _sdf_file in _sdf_files:
        # split sdf conformers
        self._split_sdf_file(_sdf_file)

  def boot(self, *args, **kwargs):
    return self

  def run(self, *args, **kwargs):
    self.run_flexx_plugin()
    return self

  def shutdown(self, *args, **kwargs):
    return self
