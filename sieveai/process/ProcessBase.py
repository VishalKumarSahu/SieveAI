from UtilityLib import ProjectManager
from ..exe import Vina, AnnapuRNA

class ProcessBase():
  def __init__(self, *args, **kwargs):
    if not kwargs.get("utility"):
      # This will add db in first project dir only
      kwargs["db_path"] = f"{kwargs.get('path_base')[0]}/SieveAI.db"
      _utility = ProjectManager(**kwargs)
      _utility.db_connect()
    else:
      _utility = getattr(self, "utility")

    self.utility = _utility

    self.__defaults = {
      "path_base": None,
      "file_config": "config.yml",
      "vina": Vina(**kwargs),
      "annapurna": AnnapuRNA(**kwargs),
    }

    self.utility.update_attributes(self, kwargs, self.__defaults)
    self.start_time = self.utility.time_start()
