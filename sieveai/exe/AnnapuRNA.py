from .VinaBase import VinaBase
import pandas as PD
import os as OS

class AnnapuRNA(VinaBase):
  def __init__(self, *args, **kwargs):
    super(AnnapuRNA, self).__init__(*args, **kwargs)
    self.__update_attr(*args, **kwargs)

  def __update_attr(self, *args, **kwargs):
    if not hasattr(self, "__defaults"): self.__defaults =  {
      "annapurna_path": "~/annapurna/annapurna.py",
      "annapurna_models": [
            # 'DL_basic',
            # 'DL_modern',
            # 'kNN_basic',
            'kNN_modern',
          ],
    }

    # Set all defaults
    [setattr(self, _k, self.__defaults[_k]) for _k in self.__defaults.keys() if not hasattr(self, _k)]
    self.__defaults = dict() # Unset defaults to prevent running for second time
    [setattr(self, _k, kwargs[_k]) for _k in kwargs.keys()]

  def __read_results_df(self, *args, **kwargs):
    _score_file_path = f"{self.path_base}{OS.sep}{self.file__result_score}"

    if not self.utility.check_path(_score_file_path):
      _error = "Results Files Does Not Exist for Processing ChimeraX Results. Skipping..."
      return _error

    _results = PD.read_csv(_score_file_path)
    _results["complex_id"] = _results['receptor'] + "--" + _results['ligand']
    _results["complex_uid"] = _results['complex_id'] + "--" + _results['conformer_id'].astype(str)

    if kwargs.get("_drop_conformers", False):
      _results = _results.drop_duplicates(subset='complex_id', keep="last")

    return _results

  def annapurna_merge_rescoring(self, *args, **kwargs):
    kwargs["_drop_conformers"] = False
    _results = self.__read_results_df(*args, **kwargs)
    _all_results = PD.DataFrame()

    _score_file_path = f"{self.path_base}{OS.sep}{self.file__result_score}"
    self.utility.log_info(f"Processing {_results.shape}.")

    for _complex_id in self.utility.ProgressBar(_results['complex_id'].unique()):
      _ap_complex_res_merged = PD.DataFrame()
      for _am in self.annapurna_models:
        _model_res_df = PD.DataFrame()
        _updated_cols = ["conformer_id", "complex_id", f"{_am}__im", f"{_am}__iu", f"{_am}__score_RNA-Ligand",f"{_am}__E_ligand", f"{_am}__score_ligand", f"{_am}__score"]
        _m_res_f = f"{self.dir_rescoring}/{_complex_id}/{_complex_id}.{_am}.csv"

        if self.utility.check_path(_m_res_f):
          _model_res_df = PD.read_csv(_m_res_f, sep="\t")
          _model_res_df.columns = _updated_cols
          _model_res_df['complex_id'] = _model_res_df.apply(lambda row: self.utility.file_name(row['complex_id'], num_ext = 2), axis=1)

          _model_res_df["complex_uid"] = _model_res_df['complex_id'] + "--" + _model_res_df['conformer_id'].astype(str)
          _model_res_df = _model_res_df.drop(["conformer_id", "complex_id"], axis=1)

          if _ap_complex_res_merged.empty:
            _ap_complex_res_merged = _model_res_df.copy()
          else:
            _ap_complex_res_merged = PD.merge(_ap_complex_res_merged.set_index("complex_uid"), _model_res_df, how="outer", on="complex_uid")
        else:
          self.utility.log_error(f"{_m_res_f} file doesn't exist.")

      if _all_results.empty:
        _all_results = _ap_complex_res_merged.copy()
      else:
        _all_results = PD.concat([_all_results, _ap_complex_res_merged])

      self.utility.log_info(f"Merge processed {_complex_id}.")

    _drop_columns = set(_results.columns).intersection(set(_all_results.columns))

    if "complex_uid" in _all_results.columns:
      _drop_columns.remove("complex_uid")

      if len(_drop_columns):
        _results = _results.drop(list(_drop_columns), axis = 1)

      _results = PD.merge(_results.set_index("complex_uid"), _all_results, how="outer", on="complex_uid")
      _results.to_csv(_score_file_path, index=False)

    self.utility.log_info(f"Processed and Merged {_results.shape} to {_score_file_path}. Took {self.utility.time_elapsed()}")

  def rescoring_AnnapuRNA(self, *args, **kwargs):
    self.__update_attr(**kwargs)
    self.utility.time_start()

    kwargs["_drop_conformers"] = True
    _results = self.__read_results_df(*args, **kwargs)

    if not isinstance(_results, PD.DataFrame):
      self.utility.log_info(f"{_results}. Skipping rescoring_AnnapuRNA...")
      return

    _models = " ".join(["-m " + _m for _m in self.annapurna_models])
    self.dir_rescoring = self.utility.validate_dir(f"{self.path_base}/rescoring-dl-modern-2")
    self.path_receptor_pdbqt = f"{self.path_base}/{self.dir_receptor_pdbqt}"
    self.path_docking = f"{self.path_base}/{self.dir_docking}"

    for _idx, _complex in self.utility.ProgressBar(_results.iterrows()):
      _rec = _complex['receptor']
      _complex_uid = _complex['complex_id']
      _rec_pdbqt = f"{self.path_receptor_pdbqt}/{_rec}.pdbqt"
      _out_pdbqt = f"{self.path_docking}/{_complex_uid}.result.pdbqt"
      _complex_output_dir = f"{self.dir_rescoring}/{_complex_uid}"
      if self.utility.check_path(_complex_output_dir):
        self.utility.log_info(f"{_complex_uid} already analysed.")
        continue
      _command = f"conda run -n annapurna python {self.annapurna_path} -r {_rec_pdbqt} -l {_out_pdbqt} {_models} -o {_complex_output_dir}/{_complex_uid} -s --overwrite"
      # print(_command)
      _result = self.utility.cmd_run(_command)
      self.utility.log_info(f"Rescoring processed {_complex_uid}.")

    self.utility.log_info(f"Total Time Taken %s." % self.utility.time_elapsed())
    self.annapurna_merge_rescoring(*args, **kwargs)
