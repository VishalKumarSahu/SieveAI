from .ChimeraX import ChimeraX

class MolConverter(ChimeraX):
  def __init__(self, *args, **kwargs):
    super(MolConverter, self).__init__(**kwargs)
    self.__update_attr(*args, **kwargs)

  def __update_attr(self, *args, **kwargs):
    if not hasattr(self, "__defaults"): self.__defaults =  {
      "utility": None,
    }

    # Set all defaults
    [setattr(self, _k, self.__defaults[_k]) for _k in self.__defaults.keys() if not hasattr(self, _k)]
    self.__defaults = dict() # Unset defaults to prevent running for second time
    [setattr(self, _k, kwargs[_k]) for _k in kwargs.keys()]

  def mol_split_sdf(self, *args, **kwargs):
    _file_path = args[0] if len(args) > 0 else kwargs.get("file_path")
    _ext_to = args[2] if len(args) > 2 else kwargs.get("ext_to", "sdf")
    _sdf_id_delimiter = args[3] if len(args) > 3 else kwargs.get("id_delimiter", "> <SYBYL.NAME>")
    _sdf_model_delimiter = args[4] if len(args) > 4 else kwargs.get("model_delimiter", "$$$$")

    _new_fn = None
    _last_line = ""
    _splitted_sdf_files = []
    with open(_file_path, 'r', encoding='UTF8') as _fh:
      _model_n_content = []
      for _line in _fh:
        _line = _line.strip('\n')
        if _last_line.startswith(_sdf_id_delimiter):
          _new_fn = _line
        elif _new_fn is None:
          _new_fn = _line

        if _line.startswith(_sdf_model_delimiter):
          _file_dir = self.utility.file_dir(_file_path)
          _new_fp = f"{_file_dir}/{_new_fn}.{_ext_to}"
          if not self.utility.check_path(_new_fp):
            self.utility.write(_new_fp, _model_n_content)
            self.utility.log_info(f"{_new_fp} splitted from {_file_path}.")
          else:
            self.utility.log_warning(f"{_new_fp} already exists. Skipping writing.")

          _splitted_sdf_files.append(_new_fp)
          _model_n_content = []
        else:
          _model_n_content.append(_line)

        _last_line = _line

    return _splitted_sdf_files

  def mol_validate_sdf(self, *args, **kwargs):
    _file_path = args[0] if len(args) > 0 else kwargs.get("file_path")
    _sdf_list = []
    _model_count = 0
    with open(_file_path, 'r', encoding='UTF8') as _fh:
      for _line in _fh:
        if _line.startswith("$$$$"):
          _model_count += 1

    if _model_count > 1:
      _splitted_files = self.mol_split_sdf(_file_path, **kwargs)
      _sdf_list.extend(_splitted_files)
    else:
      _sdf_list.append(_file_path)

    return _sdf_list

  def mol_split_multimodel(self, *args, **kwargs):
    _file_path = args[0] if len(args) > 0 else kwargs.get("path")
    _converted_paths = []

    _file_ext = self.utility.ext(_file_path)

    ## SDF
    if _file_ext == "sdf":
      _sdf_list = self.mol_validate_sdf(_file_path, **kwargs)
      _converted_paths.extend(_sdf_list)

    ## PDB
    if _file_ext == "pdb":
      # self.validate_pdb(_mol_path, **kwargs)
      pass

    return _converted_paths

  def mol_convert_molecules(self, *args, **kwargs):
    _dir_path = args[0] if len(args) > 0 else kwargs.get("path")
    _from = args[1] if len(args) > 1 else kwargs.get("ext_from")
    _ext_to = args[2] if len(args) > 2 else kwargs.get("ext_to")

    _file_paths = self.utility.find_files(_dir_path, _from)
    _files_to_convert = []
    _files_converted = []

    self.utility.log_info(f"Processing {len(_file_paths)} molecules for conversion to {_ext_to}.")

    # Split and get paths if there are multimodel
    for _mol_path in _file_paths:
      _files_to_convert.extend(self.mol_split_multimodel(_mol_path))

    # Convert the single file
    for _mol_path in _files_to_convert:
      _file_ext = self.utility.ext(_mol_path)
      _converted_path = self.utility.change_ext(_mol_path, _ext_to)
      _files_converted.append(_converted_path)

      if not self.utility.check_path(_converted_path):
        # Convert to the extension using OpenBabel
        self.utility.log_info(f"Writing conversion to {_converted_path}")
        _command = f"""obabel -i{_file_ext} {_mol_path} -o{_ext_to} -O {_converted_path} """
        _output = self.utility.cmd_run(_command)

    return _files_converted
