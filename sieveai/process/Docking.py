from .ProcessBase import ProcessBase

class Docking(ProcessBase):
  def __init__(self, *args, **kwargs):
    super(Docking, self).__init__(*args, **kwargs)
    self.utility.update_attributes(self, kwargs)

  def process(self, *args, **kwargs):
    kwargs["utility"] = getattr(self, "utility")
    kwargs_copied = kwargs.copy()
    self.utility.update_attributes(self, kwargs)

    if isinstance(self.path_base, (list, set, tuple)):
      for _base in self.path_base:
        kwargs_copied["path_base"] = _base
        self.vina.add(**kwargs_copied)
        self.vina.run()

    elif isinstance(self.path_base, str):
      self.vina.add(**kwargs_copied)
      self.vina.run()

    else:
      raise Exception("There was some problem accessing base path.")
