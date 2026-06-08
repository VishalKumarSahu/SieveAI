"""Configuration loading - SieveAI style with DictConfig."""

from UtilityLib.lib.obj import ObjDict as DictConfig
from UtilityLib.lib.path import EntityPath
import tomllib
from pathlib import Path
from typing import Optional


class ConfigManager:
    """Load and manage SieveAI configuration."""
    
    DEFAULTS = {
        "environment": "local",
        "paths": {
            "base_dir": None,
            "receptors_dir": "receptors",
            "ligands_dir": "ligands",
            "docking_dir": "docking",
            "analysis_dir": "analysis",
            "results_dir": "results",
        },
        "workflow": {
            "order": ["structuresync", "docking", "analysis", "results"],
        },
        "parallel": {
            "backend": "local",  # local, slurm, pbs
            "max_workers": 4,
            "chunk_size": 10,
            "slurm": {
                "array_size": 1000,
                "concurrent": 50,
                "cpus_per_task": 4,
                "mem": "8G",
                "time": "24:00:00",
            },
        },
        "plugins": {},
        "logging": {
            "level": "info",
            "file": None,
        },
        "checkpoint": {
            "enabled": True,
            "file": "sieveai.checkpoint.sqlite",
        },
        "timeouts": {
            "task_timeout": 3600,
            "max_retries": 3,
        },
    }
    
    def __init__(self, config_path: Optional[EntityPath] = None):
        self.config_path = config_path
        self.config = DictConfig()
        self._load_defaults()
    
    def _load_defaults(self):
        """Load default configuration."""
        self.config.update(self.DEFAULTS)
    
    def load(self, path: EntityPath) -> DictConfig:
        """Load configuration from TOML file."""
        if not path.exists():
            return self.config
        
        with open(path, "rb") as f:
            data = tomllib.load(f)
        
        self._merge(data)
        return self.config
    
    def _merge(self, data: dict):
        """Merge loaded config with defaults."""
        for key, value in data.items():
            if isinstance(value, dict) and key in self.config:
                if isinstance(self.config[key], dict):
                    self.config[key].update(value)
                else:
                    self.config[key] = value
            else:
                self.config[key] = value
    
    def save(self, path: EntityPath):
        """Save configuration to TOML file."""
        import tomli_w
        
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "wb") as f:
            tomli_w.dump(dict(self.config), f)
    
    def get_plugin_config(self, plugin_uid: str) -> DictConfig:
        """Get configuration for specific plugin."""
        if plugin_uid not in self.config.plugins:
            self.config.plugins[plugin_uid] = {
                "enabled": True,
                "settings": {},
            }
        return self.config.plugins[plugin_uid]
    
    def resolve_paths(self, base_dir: EntityPath):
        """Resolve all paths relative to base_dir."""
        self.config.paths.base_dir = base_dir
        
        for key, dir_name in self.config.paths.items():
            if key == "base_dir":
                continue
            if isinstance(dir_name, str) and not dir_name.startswith("/"):
                self.config.paths[key] = base_dir / dir_name
