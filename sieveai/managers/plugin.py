import importlib as ImportLib

class PluginManager:
  plugin_map = {}  # Global dictionary to hold plugin references

  def __init__(self, *args, **kwargs):
    pass

  @classmethod
  def add_plugin(cls, name, filename):
    cls.plugin_map[name] = filename

  @staticmethod
  def get_plugin_refs(*args, **kwargs):
    _p_map = kwargs.get('plugin_map',  args[0] if len(args) > 0 else {})

    if not isinstance(_p_map, dict):
      _p_map = {}

    PluginManager.plugin_map = _p_map

    _plugin_refs = {}
    for _pc, _pn in _p_map.items():
      _pref = PluginManager.share_plugin(_pc)
      _plugin_refs[_pc] = _pref

      # When plugins are referred by the filename or referred with lower class
      _plugin_refs[_pc.lower()] = _pref
      _plugin_refs[_pn] = _pref
      _plugin_refs[_pn.lower()] = _pref

    return _plugin_refs

  @staticmethod
  def share_plugin(_plugin_name):
    """
    @ToDo: Implement import from file or different path
    """
    _plugin_filename = PluginManager.plugin_map.get(_plugin_name)
    if not _plugin_filename:
      return None
    _plugin_ref = ImportLib.import_module("..plugins.%s" % _plugin_filename, package=__package__)
    return getattr(_plugin_ref, _plugin_name)
