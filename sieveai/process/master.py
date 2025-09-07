from .base import PluginBase

class Master(PluginBase):
  def __init__(self, *args, **kwargs):
    super().__init__(**kwargs)

  def process(self, *args, **kwargs):
    self.update_attributes(self, kwargs)

    self.log_debug('SieveAI: Master Process Started.')

    self.register_reporter()
    self.TASKS.start_step('MasterProcess', 'MasterProcess', self.name)

    for _wf_step in self.iterate(self.SETTINGS.user.workflow_order):
      _msg = f"Workflow Step {_wf_step.upper()}"

      _wf_plugins = self.SETTINGS.user.workflow[_wf_step]

      if self.is_iterable(_wf_plugins) and len(_wf_plugins) > 0:
        self.log_info(_msg, hr=True)

      for _plugin_uid in self.iterate(_wf_plugins):
        _plugin_args = {
          "path_base": self.path_base,
          "SETTINGS": self.SETTINGS,
          "current_assignment": _wf_step
        }

        _tmp_plg_ref = self.SETTINGS.PLUGIN_REFS[_plugin_uid]
        if _tmp_plg_ref is None or not callable(_tmp_plg_ref):
          self.log_error(f"MASTER_00: Plugin {_plugin_uid} not found for Workflow Step {_wf_step.upper()}!")
          continue

        _tmp_plg_ref = _tmp_plg_ref(**_plugin_args)

        self.log_debug(f"MASTER_01: Adding Plugin Dependency {_wf_step.upper()}:{_plugin_uid}...")

        self.log_info(f"Delegating Task to Plugin: {_plugin_uid}...")
        _tmp_plg_ref.boot()
        _tmp_plg_ref.run()
        _tmp_plg_ref.shutdown()

    self.TASKS.end_step('MasterProcess', 'MasterProcess', self.name)
