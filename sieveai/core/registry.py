"""Plugin registry with single-file discovery and YAML frontmatter parsing."""

import re
import yaml
import importlib.util
from pathlib import Path
from typing import Dict, Type, Any, Optional
from .base import PluginBase


class PluginRegistry:
    """
    Discovers and registers plugins from a directory.
    
    Plugins are single .py files with YAML frontmatter in docstring:
    
    ```python
    \"\"\"
    ---
    uid: vina
    name: AutoDock VINA
    version: 1.0.0
    assignments: [docking]
    requires: [structuresync]
    provides: [docked_complexes]
    executables: [vina, prepare_receptor]
    python_packages: [biopython, rdkit]
    parallel:
      strategy: compound
      chunk_size: 10
    ---
    \"\"\"
    ```
    """
    
    def __init__(self, plugins_dir: Path):
        self.plugins_dir = Path(plugins_dir)
        self._plugins: Dict[str, Dict[str, Any]] = {}
    
    def discover(self) -> Dict[str, Dict[str, Any]]:
        """
        Find all plugin files and extract metadata.
        
        Returns:
            Dict mapping plugin uid to plugin info:
            {
                "vina": {
                    "class": VinaPlugin,
                    "file": Path(...),
                    "metadata": {...},
                }
            }
        """
        if not self.plugins_dir.exists():
            return {}
        
        for py_file in self.plugins_dir.glob("*.py"):
            if py_file.name == "base.py":
                continue
            
            try:
                metadata = self._parse_metadata(py_file)
                plugin_class = self._import_plugin(py_file)
                
                # Prefer YAML metadata, fallback to class attributes
                uid = metadata.get("uid") or plugin_class.uid
                if not uid:
                    uid = py_file.stem
                
                self._plugins[uid] = {
                    "class": plugin_class,
                    "file": py_file,
                    "metadata": metadata,
                    "name": metadata.get("name") or plugin_class.name or uid,
                    "version": metadata.get("version") or plugin_class.version or "0.0.0",
                    "assignments": metadata.get("assignments") or plugin_class.assignments,
                    "requires": metadata.get("requires") or plugin_class.requires,
                    "provides": metadata.get("provides") or plugin_class.provides,
                    "executables": metadata.get("executables") or plugin_class.executables,
                    "python_packages": metadata.get("python_packages") or plugin_class.python_packages,
                    "parallel": metadata.get("parallel") or {
                        "strategy": plugin_class.parallel_strategy,
                        "chunk_size": plugin_class.chunk_size,
                    },
                }
            except Exception as e:
                print(f"[WARNING] Failed to load plugin {py_file}: {e}")
        
        return self._plugins
    
    def _parse_metadata(self, path: Path) -> Dict[str, Any]:
        """Extract YAML frontmatter from docstring."""
        content = path.read_text()
        
        # Match YAML block in docstring: """---\n...\n---"""
        match = re.search(r'^[\'"]{3}\s*---\n(.*?)\n---\s*[\'"]{3}', content, re.DOTALL)
        if match:
            try:
                return yaml.safe_load(match.group(1)) or {}
            except yaml.YAMLError:
                return {}
        
        return {}
    
    def _import_plugin(self, path: Path) -> Type[PluginBase]:
        """Dynamically import plugin class from file."""
        spec = importlib.util.spec_from_file_location(path.stem, path)
        if not spec or not spec.loader:
            raise ImportError(f"Cannot load spec from {path}")
        
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        
        # Find class that inherits PluginBase
        for name in dir(module):
            obj = getattr(module, name)
            if isinstance(obj, type) and issubclass(obj, PluginBase) and obj is not PluginBase:
                return obj
        
        raise ImportError(f"No PluginBase subclass found in {path}")
    
    def get(self, uid: str) -> Optional[Dict[str, Any]]:
        """Get plugin info by uid."""
        return self._plugins.get(uid)
    
    def list_plugins(self) -> list:
        """List all discovered plugins."""
        return [
            {
                "uid": uid,
                "name": info["name"],
                "version": info["version"],
                "assignments": info["assignments"],
                "requires": info["requires"],
                "file": str(info["file"]),
            }
            for uid, info in self._plugins.items()
        ]
