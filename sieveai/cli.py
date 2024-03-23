def process_args():
  import os as OS
  from sieveai import __program__, __build__
  from UtilityLib import CommandUtility

  # key: (['arg_k1', 'arg_k2'], nargs, default, help, {})
  _version_info = f"{__program__} (build-{__build__})"
  _cli_settings = {
    "debug": (['--debug'], None, 0, 'silent/verbose/debug mode from 0, 1, 2, and 3.', {}),
    "db_path": (['-db'], None, None, 'Provide path to the database for Sieve project.', {}),
    "path_base": (['-b'], "*", [OS.getcwd()], 'Provide base directory to run the process.', {}),
    "path_receptor": (['-r'], None, None, 'Specify path of directory containing receptors.', {}),
    "dir_receptor": (['-dr'], None, "receptor", 'Specify name of the directory containing receptor PDB structures under base directory.', {}),
    "ext_receptor": (['-er'], "*", ["*.pdb"], 'Pattern(s) for extension of receptor files.', {}),
    "path_ligand": (['-l'], None, None, 'Specify path of directory containing ligands.', {}),
    "dir_ligand": (['-dl'], None, "ligand", 'Specify name of the directory containing ligand PDB structures under base directory.', {}),
    "ext_ligand": (['-el'], "*", ["*.pdb", "*.sdf", "*.mol2"], 'Pattern(s) for extension of ligand files.', {}),
    "multiprocess": (['-m'], None, False, 'Multiprocessing of the complexes.', {}),
    "mode": ([], None, "prod", 'Other options prod|dev|test.', {}),
    # "action": (['-a'], None, 'docking', 'Action to perform among docking|rescoring|web.', {}),
  }

  _cmdu = CommandUtility()
  _params = _cmdu.get_cli_args(_cli_settings, version=_version_info)

  print("{}\n{}\n{}".format("=" * len(_version_info), _version_info, "=" * len(_version_info)))

  return _params

def dock():
  _args = process_args()
  from . import Docking
  _process = Docking(**_args)
  print(f"Initalizing Docking...")
  _process.process()

def rescore():
  _args = process_args()
  from . import Rescoring
  _process = Rescoring(**_args)
  print(f"Initalizing Rescoring...")
  _process.process()
