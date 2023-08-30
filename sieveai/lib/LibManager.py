from .MolConverter import MolConverter

class LibManager(MolConverter):
  def __init__(self, *args, **kwargs):
    super(LibManager, self).__init__(**kwargs)
    self.__update_attr(*args, **kwargs)

  def __update_attr(self, *args, **kwargs):
    if not hasattr(self, "__defaults"): self.__defaults =  {}

    # Set all defaults
    [setattr(self, _k, self.__defaults[_k]) for _k in self.__defaults.keys() if not hasattr(self, _k)]
    self.__defaults = dict() # Unset defaults to prevent running for second time
    [setattr(self, _k, kwargs[_k]) for _k in kwargs.keys()]
