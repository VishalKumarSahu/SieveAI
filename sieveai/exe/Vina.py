import os as OS
from pdbtools.pdb_selaltloc import select_by_occupancy as selAltLoc

from .VinaBase import VinaBase

class Vina(VinaBase):
  def __init__(self, *args, **kwargs):
    super(Vina, self).__init__(*args, **kwargs)

  def prepare_receptor(self, *args, **kwargs):
    # return # to by pass
    self.utility.log_info("Preparing Receptors")

    if self.path_receptor is None and self.dir_receptor:
      self.path_receptor = f"{self.path_base}/{self.dir_receptor}"

    if not self.path_receptor:
      raise Exception("Receptor directory is not provided.")

    if self.path_receptor_pdbqt is None and self.dir_receptor_pdbqt:
      self.path_receptor_pdbqt = self.utility.validate_dir(f"{self.path_base}/{self.dir_receptor_pdbqt}")

    if self.path_receptor_clean is None and self.dir_receptor_clean:
      self.path_receptor_clean = self.utility.validate_dir(f"{self.path_base}/{self.dir_receptor_clean}")

    if self.path_receptor_summary is None and self.dir_receptor_summary:
      self.path_receptor_summary = self.utility.validate_dir(f"{self.path_base}/{self.dir_receptor_summary}")

    if self.path_receptor_config is None and self.dir_receptor_config:
      self.path_receptor_config = self.utility.validate_dir(f"{self.path_base}/{self.dir_receptor_config}")


    _rec_paths = self.utility.find_files(self.path_receptor, '*.pdb')
    _rec_pre_paths = self.utility.find_files(self.path_receptor_pdbqt, "*.pdbqt")

    _temp_rec_paths = {self.utility.filename(_f) for _f in _rec_paths}
    _temp_rec_pre_paths = {self.utility.filename(_f) for _f in _rec_pre_paths}
    _rec_not_prepared = _temp_rec_paths.difference(_temp_rec_pre_paths)

    if len(_rec_not_prepared) == 0:
      self.utility.log_warning(f"All {len(_rec_paths)} receptors are already prepared. Continuing...")
      return
    else:
      self.utility.log_info(f"{len(_rec_not_prepared)}/{len(_rec_paths)} receptors are not prepared. Preparing...")

    for _rec_file in _rec_not_prepared:
      _rec_path = f"{self.path_receptor}/{_rec_file}.pdb"
      _rec_name = self.utility.file_name(_rec_path, with_ext = True)
      if self.utility.file_name(_rec_path) in self.skipped_receptor:
        self.utility.log_warning("Receptor in skipped list.")
        continue

      if self.utility.check_path(f"{self.path_receptor_pdbqt}{OS.sep}{_rec_name}qt"):
        self.utility.log_warning("Receptor already prepared. Delete the prepared receptor file. @TODO Force or selective override yet to be implemented.")
        continue

      # Clean PDB
      lines = []
      with open(_rec_path, 'r', encoding="utf8") as file_handle:
        lines = list(selAltLoc(file_handle))

      _atom_records = [x for x in lines if not x.startswith("HETATM")]
      _het_records = [x for x in lines if not x.startswith("ATOM") and not x.startswith("TER")]

      _clean_receptor =  f"{self.path_receptor_clean}{OS.sep}{_rec_name}"

      self.utility.write(_clean_receptor, _atom_records)
      self.utility.log_info(f"{_rec_name} contains ATOM {len(_atom_records)} and {len(_het_records)} non-ATOM records.")

      command = f"""prepare_receptor -r "{self.path_receptor_clean}{OS.sep}{_rec_name}" -o "{self.path_receptor_pdbqt}{OS.sep}{_rec_name}qt" -A 'bonds_hydrogens' -U 'waters' -v -d "{self.path_receptor_summary}{OS.sep}{_rec_name}.summary.log" """

      # @TODO: Process Using CommandManager
      s = OS.popen(command)
      output = s.read()
      self.utility.log_info(output)

    self.utility.log_info("Receptor preparation completed.")

  def prepare_ligand(self, *args, **kwargs):
    # Check if path exists or not and create accordingly
    if self.path_ligand is None and self.dir_ligand:
      self.path_ligand = f"{self.path_base}/{self.dir_ligand}"

    if self.path_ligand_pdbqt is None and self.dir_ligand_pdbqt:
      self.path_ligand_pdbqt = self.utility.validate_dir(f"{self.path_base}/{self.dir_ligand_pdbqt}")

    if not self.path_ligand:
       raise Exception("Ligand directory is not provided.")

    self.mol_convert_molecules(self.path_ligand, self.ext_ligand, "pdb")
    _lig_paths = self.utility.find_files(self.path_ligand, "*.pdb")
    _lig_pre_paths = self.utility.find_files(self.path_ligand_pdbqt, "*.pdbqt")

    _temp_lig_paths = {self.utility.filename(_f) for _f in _lig_paths}
    _temp_lig_pre_paths = {self.utility.filename(_f) for _f in _lig_pre_paths}
    _lig_not_prepared = _temp_lig_paths.difference(_temp_lig_pre_paths)

    if len(_lig_not_prepared) == 0:
      self.utility.log_warning(f"All {len(_lig_paths)} ligand(s) are already prepared. Continuing...")
      return
    else:
      self.utility.log_info(f"{len(_lig_not_prepared)}/{len(_lig_paths)} ligand(s) are not prepared. Preparing...")

    for _lig_file in _lig_not_prepared:
      _lig_path = f"{self.path_ligand}{OS.sep}{_lig_file}.pdb"
      _lig_file_name = self.utility.file_name(_lig_path)
      _lig_pdbqt_path = f"{self.path_ligand_pdbqt}{OS.sep}{_lig_file_name}.pdbqt"

      if self.utility.file_name(_lig_path) in self.skipped_ligand:
        self.utility.log_info("Ligand is in skipped list.")
        continue

      # Continue if already created
      if self.utility.check_path(_lig_pdbqt_path):
        continue

      # Directly convert pdb to pdbqt using openbabel
      _res = self.utility.cmd_run([
          "obabel", _lig_path, "-O", _lig_pdbqt_path, "-h",
          "--quiet"
        ])

      self.utility.log_info(f"Converted ligand {_lig_file_name}.\n{_res}", log_wrap = False)

    self.utility.log_info("Ligand preparation completed.")
