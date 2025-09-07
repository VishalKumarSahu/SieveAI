from ..process.converter import PluginConverterBase

class MGLTools(PluginConverterBase):
  is_ready = False
  plugin_name = "Meeko"
  plugin_uid = "Meeko"
  plugin_version = "1.0"
  assignments = ['conversion']
  current_assignment = None
  url = "https://github.com/forlilab/Meeko"

  def __init__(self, *args, **kwargs):
    super().__init__(**kwargs)

  def setup(self, *args, **kwargs):
    self.update_attributes(self, kwargs)

  def prepare_receptor(self, *args, **kwargs) -> None:
    _path_source = kwargs.get('path_source', args[0] if len(args) > 0 else None)
    _path_target = kwargs.get('path_target', args[1] if len(args) > 1 else None)
    _path_log = kwargs.get('path_log', args[2] if len(args) > 2 else (self.path_target.with_suffix('.summary.log')))

    if not _path_target.exists():
      self.cmd_run(*[
          "prepare_receptor",
          '-r', _path_source,
          '-o', _path_target,
          '-A', "bonds_hydrogens",
          '-U', "waters",
          '-v', '-d', _path_log
        ])

  def mk_prepare_receptor(self, *args, **kwargs):
    _path_source = kwargs.get('path_source', args[0] if len(args) > 0 else None)
    _path_target = kwargs.get('path_target', args[1] if len(args) > 1 else None)
    self.cmd_run(*[
      "mk_prepare_receptor.py", '--skip_gpf',
      ], **{
        '--pdb'     : _path_source,
        '-o'        : _path_target,
        '--box_size': "70 70 70",
        'text'      : False,
        'check'     : False,
        'shellx'    : False,
      })

  def convert_pqbqt(self, *args, **kwargs):
    _path_source = kwargs.get('path_source', args[0] if len(args) > 0 else None)
    _path_target = kwargs.get('path_target', args[1] if len(args) > 1 else None)

    _format_to = kwargs.get('format_to', args[2] if len(args) > 2 else None)
    _format_from = kwargs.get('format_from', args[3] if len(args) > 3 else None)
    _addH = kwargs.get('addH', args[4] if len(args) > 4 else True)
    _removeH = kwargs.get('removeH', args[4] if len(args) > 4 else True)
    try:
      from meeko import MoleculePreparation, PDBQTWriterLegacy
      from rdkit import Chem
      _path_source = str(_path_source)
      _molecules = Chem.SDMolSupplier(_path_source, removeHs=_removeH, ) if 'sdf' in _format_from else [Chem.MolFromPDBFile(_path_source)]

      # there is one molecule in this SD file, this loop iterates just once
      for _mol in _molecules:
        _preparator = MoleculePreparation()
        _mol_setups = _preparator.prepare(_mol)
        for _setup in _mol_setups:
          # _setup.show() # optional
          _pdbqt_string, _, _ = PDBQTWriterLegacy.write_string(_setup)
          self.write(_path_target, _pdbqt_string)
    except Exception as _e:
      self.log_error(f'Error in MGLTools: {_e}')
