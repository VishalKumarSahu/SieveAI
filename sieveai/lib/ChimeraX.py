import re as REGEX
import pandas as PD

from .OpenBabel import OpenBabel

class ChimeraX(OpenBabel):
  def __init__(self, *args, **kwargs):
    super(ChimeraX, self).__init__(**kwargs)
    self.__update_attr(*args, **kwargs)

  def __update_attr(self, *args, **kwargs):
    if not hasattr(self, "__defaults"): self.__defaults =  {
        "aa__polar": ["SER", "THR", "CYS", "ASN", "GLN", "TYR"],
      }
    # Set all defaults
    [setattr(self, _k, self.__defaults[_k]) for _k in self.__defaults.keys() if not hasattr(self, _k)]
    self.__defaults = dict() # Unset defaults to prevent running for second time
    [setattr(self, _k, kwargs[_k]) for _k in kwargs.keys()]

  def __parse_hbond_line(self, line):
    """
    @BUG: HBonds parsing from ChimeraX results needs to recheck
    @TODO: Handle case when hydrogen is not available in the atomic model
    """
    result = None
    if "#" in line:
      pass
    elif "/" in line:
      # Since Model ID is 1 and only Chain ID is given so appending to keep the syntax intact
      line = str(line).replace("/", "Structure #1/")
    else:
      raise Exception("Unknown type of col separators.")

    cols = line.split("#")

    if len(cols):
      try:
        # Use Second Column to get Donor
        donor = tuple(cols[1].split()[:4])
        # Use Third Column to get acceptor
        acceptor = tuple(cols[2].split()[:4])
        # Use fourth colun to get Hydrogen and HBond length
        hydrogen = tuple(cols[3].split()[:4])
        length = tuple(cols[3].split()[-2:])
        result = {
          "donor": donor,
          "acceptor": acceptor,
          "hydrogen": hydrogen,
          "length": length
        }
      except:
        length = None
        print(line)

    return result

  def chimerax_parse_hbonds_file(self, file_url):
    hbond_attribs = []

    if self.utility.check_path(file_url):
      key_term = "H-bonds"
      search_regex = REGEX.compile(f"(\d+) {key_term}")
      with open(file_url, 'r', encoding="utf8") as f:
        hbonds = f.readlines()

      residues_line_flag = False

      for line in hbonds:
        # print(residues_line_flag, search_regex.search(line))
        if residues_line_flag and not line.startswith(key_term) and "#" in line:
          _res = self.__parse_hbond_line(line)
          if _res:
            hbond_attribs.append(_res)
        if not residues_line_flag and search_regex.search(line):
          residues_line_flag = True
          if not int(search_regex.search(line).group(1)) > 0:
            return hbond_attribs

    return hbond_attribs

  def __parse_contacts_line(self, line):
    result = None
    if "#" in line:
      pass
    elif "/" in line:
      # Since Model ID is 1 and only Chain ID is given so appending to keep the syntax intact
      line = str(line).replace("/", "Structure #1/")
    else:
      raise Exception("Unknown type of col separators.")

    cols = line.split("#")
    # print(cols); raise Exception("Unknown type of col separators.")

    if len(cols):
      # Use Second Column to get Atom1
      atom1 = tuple(cols[1].split()[:4])
      # Use Third Column to get Atom2
      atom2 = tuple(cols[2].split()[:4])
      # Use fourth colun to get Hydrogen and HBond length
      overlap = cols[2].split()[-2]
      distance = cols[2].split()[-1]
      result = {
        "atom1": atom1,
        "atom2": atom2,
        "overlap": overlap,
        "distance": distance
      }

    return result

  def chimerax_parse_contacts_file(self, file_url):
    contact_attribs = []
    if self.utility.check_path(file_url):
      key_term = "contacts"
      search_regex = REGEX.compile(f"(\d+) {key_term}")
      with open(file_url, 'r', encoding="utf8") as f:
        contacts = f.readlines()

      residues_line_flag = False

      for line in contacts:
        if residues_line_flag and not len({"atom1", "atom2"}.intersection(set(line.split()))):
          _res = self.__parse_contacts_line(line)
          if _res:
            contact_attribs.append(_res)
        if not residues_line_flag and search_regex.search(line):
          residues_line_flag = True
          if not int(search_regex.search(line).group(1)) > 0:
            return contact_attribs

    return contact_attribs

  def __parse_atom_identity(self, atom_details: tuple = ()):
    if atom_details and type(atom_details) is tuple and len(atom_details) == 4:
      atom_details_model = atom_details[0].split("/")[0]
      atom_details_chain = atom_details[0].split("/")[1]
      atom_details_resname = atom_details[1]
      atom_details_resid = atom_details[2]
      atom_details_atom = atom_details[3]

      model_id = atom_details_model
      sub_model_id = None
      if "." in atom_details_model:
        model_id = atom_details_model.split(".")[0]
        sub_model_id = atom_details_model.split(".")[1]

      __result = {
        "model_id": model_id,
        "sub_model_id": sub_model_id,
        "chain": atom_details_chain,
        "resname": atom_details_resname,
        "resid": atom_details_resid,
        "atom": atom_details_atom,
      }

      return __result
    else:
      raise Exception(f"Problem in atom records {atom_details}")

  def chimerax_get_contacts_residues(self, file_url):
    __contacts = PD.DataFrame(self.chimerax_parse_contacts_file(file_url))

    if not __contacts.shape[0]:
        return None

    __contacts["atom1"] = __contacts["atom1"].apply(self.__parse_atom_identity)
    __contacts["atom2"] = __contacts["atom2"].apply(self.__parse_atom_identity)

    __contacts = PD.concat([__contacts.drop(['atom1'], axis=1), __contacts['atom1'].apply(PD.Series)], axis=1)
    __contacts.rename(columns={'model_id': 'atom1__model_id', 'sub_model_id': 'atom1__sub_model_id', 'chain': 'atom1__chain', 'resname': 'atom1__resname', 'resid': 'atom1__resid', 'atom': 'atom1__atom'}, inplace=True)

    __contacts = PD.concat([__contacts.drop(['atom2'], axis=1), __contacts['atom2'].apply(PD.Series)], axis=1)
    __contacts.rename(columns={'model_id': 'atom2__model_id', 'sub_model_id': 'atom2__sub_model_id', 'chain': 'atom2__chain', 'resname': 'atom2__resname', 'resid': 'atom2__resid', 'atom': 'atom2__atom'}, inplace=True)

    return __contacts

  def chimerax_get_hbonds_residues(self, file_url):
    __hbonds = PD.DataFrame(self.chimerax_parse_hbonds_file(file_url))

    if not __hbonds.shape[0]:
        return None

    __hbonds["donor"] = __hbonds["donor"].apply(self.__parse_atom_identity)
    __hbonds["acceptor"] = __hbonds["acceptor"].apply(self.__parse_atom_identity)
    __hbonds["hydrogen"] = __hbonds["hydrogen"].apply(self.__parse_atom_identity)

    __hbonds = PD.concat([__hbonds.drop(['donor'], axis=1), __hbonds['donor'].apply(PD.Series)], axis=1)
    __hbonds.rename(columns={'model_id': 'donor__model_id', 'sub_model_id': 'donor__sub_model_id', 'chain': 'donor__chain', 'resname': 'donor__resname', 'resid': 'donor__resid', 'atom': 'donor__atom'}, inplace=True)

    __hbonds = PD.concat([__hbonds.drop(['acceptor'], axis=1), __hbonds['acceptor'].apply(PD.Series)], axis=1)
    __hbonds.rename(columns={'model_id': 'acceptor__model_id', 'sub_model_id': 'acceptor__sub_model_id', 'chain': 'acceptor__chain', 'resname': 'acceptor__resname', 'resid': 'acceptor__resid', 'atom': 'acceptor__atom'}, inplace=True)

    __hbonds = PD.concat([__hbonds.drop(['hydrogen'], axis=1), __hbonds['hydrogen'].apply(PD.Series)], axis=1)
    __hbonds.rename(columns={'model_id': 'hydrogen__model_id', 'sub_model_id': 'hydrogen__sub_model_id', 'chain': 'hydrogen__chain', 'resname': 'hydrogen__resname', 'resid': 'hydrogen__resid', 'atom': 'hydrogen__atom'}, inplace=True)

    return __hbonds

  def chimerax_run_file(self, *args, **kwargs):
    _cxc_file = args[0] if len(args) > 0 else kwargs.get("file")

    _res = self.utility.cmd_run(["chimerax",
      "--cmd", f"open {_cxc_file}",
      "--silent",
      "--offscreen",
    ])

    return _res

  def chimerax_convert(self, *args, **kwargs):
    _mol_path = args[0] if len(args) > 0 else kwargs.get("mol_path")
    _ext_to = args[1] if len(args) > 1 else kwargs.get("ext_to")
    _destination = args[2] if len(args) > 2 else kwargs.get("destination")

    _new_file_path = self.utility.change_ext(_mol_path, _ext_to)

    if _destination:
      _file_name = self.utility.filename(_mol_path)
      _new_file_name = self.utility.change_ext(_file_name, _ext_to)
      _destination = self.utility.validate_dir(_destination)
      _new_file_path = f"{_destination}/{_new_file_name}"

    if not self.utility.check_path(_new_file_path):
      self.utility.log_info(f"converting {_mol_path} to {_new_file_path}.")
      _file_content = [
          f"close; open {_mol_path}; wait 5;"
          f"save {_new_file_path}; wait 5; close;",
          "\n"
      ]
    else:
      self.utility.log_warning(f"Not converting {_mol_path} as {_new_file_path} already exists.")
