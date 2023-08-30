from UtilityLib import UtilityManager
from sieveai.lib import LibManager

class ExecutableBase(LibManager):
  def __init__(self, *args, **kwargs):
    if not kwargs.get("utility") or not hasattr(self, "utility"):
          # This will add db in first project dir only
      kwargs["db_path"] = f"{kwargs.get('path_base')[0]}/SieveAI.db"
      _utility = UtilityManager(**kwargs)
      _utility.db_connect()
      kwargs['utility'] = _utility
    self.__defaults =  {}
    _utility.update_attributes(self, kwargs)
