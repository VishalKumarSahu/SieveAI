# Provide external libraries reference from this base

from UtilityLib.lib.path import EntityPath
from UtilityLib.lib.obj import ObjDict as DictConfig
from UtilityLib.lib.task import TaskManager
from UtilityLib.lib.step import StepManager
from UtilityLib.lib.schedule import ScheduleManager
from UtilityLib import ProjectManager

from tqdm.auto import tqdm as _TQDMPB

from .__metadata__ import __version__
from .managers.plugin import PluginManager

_SM_Ref = ScheduleManager()

class SieveAIBase(ProjectManager):
  name = 'SieveAI'
  version = __version__
  EXE_MAP = DictConfig()
  SETTINGS = DictConfig()
  TASKS = TaskManager()

  path_backup_dir = None
  SchReporter = _SM_Ref
  report_callback_fn = None
  log_level = 'debug'

  def __init__(self, *args, **kwargs):
    super().__init__(**kwargs)
    self.path_backup_dir = (self.path_base / 'temp-progress-backup').validate()
    self.path_sieveai_master_config = (EntityPath('~') / '.SieveAI').validate()
    self._first_check()
    self.list_plugins()

  def __exit__(self):
    self.finalise_all()

  def register_reporter(self):
    if len(self.SchReporter.events) == 0 and self.SETTINGS.user.report_flag == True:
      self.SchReporter.add(self.report_progress, self.SETTINGS.user.report_interval, self.SETTINGS.user.report_interval_unit)
      self.SchReporter.add(self.report_check, self.SETTINGS.user.exit_interval, self.SETTINGS.user.exit_interval_unit)

  def report_check(self):
    self.log_debug('Checking if tasks are completed to stop reporting.')
    _lsd = self.TASKS.get_status_df()
    if _lsd is not None:
      _remaining = _lsd.Status.sum() - _lsd.Status.count()
      if _remaining == 0:
        self.finalise_all()
      else:
        self.log_debug(f'{_remaining} tasks remaining. Will check again in {self.SETTINGS.user.report_interval}')

  def finalise_all(self):
    self.log_info('Finalising all event scheduled.')
    if len(self.SchReporter.events) > 0:
      self.SchReporter.stop_all()

    _bkups = self.clean_file_backups(self.path_backup_dir / self.SETTINGS.user.path_sob_progress.name)

    if self.path_backup_dir.exists() and len(self.path_backup_dir.files) == 0:
      self.path_backup_dir.delete(is_protected=False)

  _report_progress_bar_ref = _TQDMPB(desc='SieveAI-Jobs')
  def report_progress(self, *args, **kwargs):
    _status_df = self.TASKS.get_status_df()
    _done, _total = _status_df.Status.sum(), _status_df.Status.count()
    self._report_progress_bar_ref.n = _done
    self._report_progress_bar_ref.total = _total
    self._report_progress_bar_ref.refresh()
    if callable(self.report_callback_fn):
      self.report_callback_fn(self, *args, **kwargs)

  def save_progress(self, *args, **kwargs):
    """
      Copy previous version with timestamp
      Pickle SETTINGS
      Collect all settings and data from different plugins (Memory Management?)
    """
    self.create_file_backup(self.SETTINGS.user.path_sob_progress, self.path_backup_dir)
    # Limit backup files to last 10
    self.pickle(self.SETTINGS.user.path_sob_progress, self.SETTINGS)

  def restore_progress(self, *args, **kwargs):
    """
      Read pickle file
      If there is any issue check SETTINGS backup and restore progress from there
    """
    if self.SETTINGS.user.path_sob_progress.exists() and self.SETTINGS.user.path_sob_progress.size > 0:
      # If it fails grab latest backup using backup pattern
      _s = self.unpickle(self.SETTINGS.user.path_sob_progress)
      self.SETTINGS.update(_s)

  def get_plugin(self, _plugin_name):
    """Share plugin reference"""
    return self.PLUGINS[_plugin_name]

  def list_plugins(self):
    # Plugins from core
    self.SETTINGS.path_plugins = EntityPath(__file__, '..', 'plugins').resolve()
    self.SYS.path.append(str(self.SETTINGS.path_plugins))

    _file_class_map = self.list_py_classes(self.SETTINGS.path_plugins)
    _plugins = {}
    for _file, _classes in _file_class_map.items():
      for _cl in _classes:
        _plugins[_cl] = _file

    self.SETTINGS.PLUGIN_REFS.update(PluginManager.get_plugin_refs(_plugins))

  def _first_check(self):
    _self_config_toml =  self.path_sieveai_master_config / 'sieveai.config.toml'

    if not _self_config_toml.exists():
      from .managers.once import OnceManager
      _om = OnceManager(path_base=self.path_base, version=self.version)
      _om.check(_self_config_toml)
    else:
      _exe_map = self.load_toml(_self_config_toml) or DictConfig()
      self.EXE_MAP = _exe_map.get('EXE_PATHS', DictConfig())
