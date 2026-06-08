"""Abstract base class for all SieveAI plugins."""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Dict, List, Optional
import shutil
import subprocess


class PluginBase(ABC):
    """
    Abstract base class for SieveAI plugins.
    
    Plugin lifecycle:
        1. validate() - Check dependencies, config, inputs
        2. init(**kwargs) - Initialize before processing
        3. run(item) - Process single item (compound, receptor, etc.)
        4. finalize(results) - Aggregate results, write outputs
        5. cleanup() - Post-run cleanup (optional)
    
    Metadata (override in subclass):
        - uid: Unique plugin identifier
        - name: Human-readable name
        - version: Plugin version
        - assignments: Workflow slots this plugin handles
        - requires: Plugins that must run before
        - provides: Outputs this plugin produces
        - executables: Required system executables
        - python_packages: Required Python packages
        - parallel_strategy: none | compound | receptor | batch
        - chunk_size: Items per parallel job
    """
    
    # Required class attributes (override in subclass)
    uid: str = ""
    name: str = ""
    version: str = "0.0.0"
    
    # Workflow
    assignments: List[str] = []
    requires: List[str] = []
    provides: List[str] = []
    
    # Dependencies
    executables: List[str] = []
    python_packages: List[str] = []
    
    # Parallel config
    parallel_strategy: str = "none"
    chunk_size: int = 1
    
    # Runtime (set by engine)
    config: Dict[str, Any] = {}
    workdir: Path = None
    logger: Any = None
    
    def __init__(self, **kwargs):
        """Initialize plugin with runtime config."""
        self.config = kwargs.get("config", {})
        self.workdir = kwargs.get("workdir", Path.cwd())
        self.logger = kwargs.get("logger")
    
    @abstractmethod
    def validate(self) -> bool:
        """
        Check if plugin can run.
        
        Returns:
            True if all dependencies are satisfied, False otherwise.
        """
        pass
    
    @abstractmethod
    def init(self, **kwargs) -> None:
        """
        Initialize before processing.
        
        Setup paths, validate inputs, prepare state.
        """
        pass
    
    @abstractmethod
    def run(self, item: Any) -> Any:
        """
        Process single item.
        
        Args:
            item: Compound dict, receptor path, or other input.
        
        Returns:
            Result dict with status and outputs.
        """
        pass
    
    def finalize(self, results: List[Any]) -> Dict[str, Any]:
        """
        Aggregate results and write outputs.
        
        Args:
            results: List of results from run() calls.
        
        Returns:
            Summary dict with status and output paths.
        """
        return {"status": "ok", "count": len(results)}
    
    def cleanup(self) -> None:
        """Post-run cleanup. Override if needed."""
        pass
    
    # Utility methods
    def _which(self, exe: str) -> Optional[Path]:
        """Find executable in PATH."""
        path = shutil.which(exe)
        return Path(path) if path else None
    
    def check_executable(self, exe: str) -> bool:
        """Check if executable exists."""
        return self._which(exe) is not None
    
    def _import(self, pkg: str) -> bool:
        """Check if Python package is importable."""
        try:
            __import__(pkg)
            return True
        except ImportError:
            return False
    
    def _run_command(
        self,
        cmd: List[str],
        cwd: Optional[Path] = None,
        timeout: Optional[int] = None,
    ) -> subprocess.CompletedProcess:
        """Run shell command."""
        return subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=cwd or self.workdir,
            timeout=timeout,
        )
    
    def _log(self, level: str, msg: str) -> None:
        """Log message if logger available."""
        if self.logger:
            getattr(self.logger, level)(msg)
        else:
            print(f"[{level.upper()}] {msg}")
