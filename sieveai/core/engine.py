"""Main SieveAI engine - orchestrates plugins and workflow execution."""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional
from datetime import datetime

from .config import SieveConfig
from .registry import PluginRegistry
from .environment import EnvironmentDetector, EnvironmentType
from ..plugins.base import PluginBase


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("sieveai.engine")


class SieveEngine:
    """
    Main SieveAI orchestration engine.
    
    Lifecycle:
        1. initialize() - Load config, discover plugins
        2. validate() - Check all plugin dependencies
        3. plan() - Build workflow DAG, determine execution order
        4. execute() - Run workflow with parallel backend
        5. finalize() - Aggregate results, cleanup
    """
    
    def __init__(self, config: Optional[SieveConfig] = None):
        self.config = config or SieveConfig()
        self.registry: Optional[PluginRegistry] = None
        self.plugins: Dict[str, Dict[str, Any]] = {}
        self.workflow_order: List[str] = []
        self.results: Dict[str, Any] = {}
        self.start_time: Optional[datetime] = None
        self.end_time: Optional[datetime] = None
    
    def initialize(self, plugins_dir: Optional[Path] = None) -> None:
        """
        Initialize engine: load config, discover plugins.
        
        Args:
            plugins_dir: Path to plugins directory (default: sieveai/plugins)
        """
        logger.info("Initializing SieveAI engine...")
        
        # Detect environment
        env = EnvironmentDetector.detect()
        logger.info(f"Environment: {env.value}")
        
        if env != EnvironmentType.LOCAL:
            env_info = EnvironmentDetector.get_info()
            logger.info(f"Environment info: {env_info['environment']}")
        
        # Initialize plugin registry
        if plugins_dir is None:
            plugins_dir = Path(__file__).parent.parent / "plugins"
        
        self.registry = PluginRegistry(plugins_dir)
        self.plugins = self.registry.discover()
        
        logger.info(f"Discovered {len(self.plugins)} plugins: {list(self.plugins.keys())}")
        
        # Set default paths
        if self.config.receptors_dir is None:
            self.config.receptors_dir = self.config.base_dir / "receptors"
        if self.config.ligands_dir is None:
            self.config.ligands_dir = self.config.base_dir / "ligands"
        if self.config.docking_dir is None:
            self.config.docking_dir = self.config.base_dir / "docking"
        if self.config.analysis_dir is None:
            self.config.analysis_dir = self.config.base_dir / "analysis"
        if self.config.results_dir is None:
            self.config.results_dir = self.config.base_dir / "results"
        
        # Create directories
        for dir_path in [
            self.config.receptors_dir,
            self.config.ligands_dir,
            self.config.docking_dir,
            self.config.analysis_dir,
            self.config.results_dir,
        ]:
            dir_path.mkdir(parents=True, exist_ok=True)
        
        logger.info("Initialization complete.")
    
    def validate(self) -> Dict[str, Any]:
        """
        Validate all plugins and dependencies.
        
        Returns:
            Dict with validation results per plugin.
        """
        logger.info("Validating plugins...")
        
        results = {}
        for uid, info in self.plugins.items():
            plugin_class = info["class"]
            plugin = plugin_class(
                config=self.config.plugins.get(uid, {}).settings,
                workdir=self.config.base_dir,
                logger=logger,
            )
            
            is_valid = plugin.validate()
            results[uid] = {
                "valid": is_valid,
                "name": info["name"],
                "version": info["version"],
                "executables": info["executables"],
                "python_packages": info["python_packages"],
            }
            
            if is_valid:
                logger.info(f"✓ Plugin {uid} ({info['name']}) validated")
            else:
                logger.warning(f"✗ Plugin {uid} ({info['name']}) validation failed")
        
        return results
    
    def plan(self) -> List[str]:
        """
        Build workflow execution order based on plugin dependencies.
        
        Returns:
            Ordered list of plugin uids to execute.
        """
        from ..workflow.dag import WorkflowDAG
        
        logger.info("Planning workflow...")
        
        dag = WorkflowDAG(self.plugins)
        self.workflow_order = dag.build()
        
        logger.info(f"Workflow order: {' → '.join(self.workflow_order)}")
        return self.workflow_order
    
    def execute(self, inputs: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """
        Execute workflow.
        
        Args:
            inputs: Workflow inputs (compounds, receptors, etc.)
        
        Returns:
            Dict with results from each plugin.
        """
        logger.info("Executing workflow...")
        self.start_time = datetime.now()
        
        inputs = inputs or {}
        self.results = {}
        
        for uid in self.workflow_order:
            if uid not in self.plugins:
                logger.warning(f"Plugin {uid} not found, skipping...")
                continue
            
            info = self.plugins[uid]
            plugin_class = info["class"]
            
            # Check if plugin is enabled
            plugin_config = self.config.plugins.get(uid, {})
            if not plugin_config.enabled:
                logger.info(f"Plugin {uid} disabled, skipping...")
                continue
            
            # Initialize plugin
            plugin = plugin_class(
                config=plugin_config.settings,
                workdir=self.config.base_dir,
                logger=logger,
            )
            
            # Run plugin
            logger.info(f"Running plugin: {uid} ({info['name']})...")
            
            # Initialize with inputs
            plugin.init(**inputs)
            
            # Execute with parallel backend
            if info["parallel"]["strategy"] != "none":
                results = self._execute_parallel(plugin, inputs, info["parallel"])
            else:
                results = plugin.run(inputs)
            
            # Finalize
            if isinstance(results, list):
                summary = plugin.finalize(results)
            else:
                summary = plugin.finalize([results])
            
            self.results[uid] = {
                "status": "ok",
                "summary": summary,
            }
            
            # Pass outputs to next plugin
            inputs.update(summary)
            
            # Cleanup
            plugin.cleanup()
        
        self.end_time = datetime.now()
        duration = (self.end_time - self.start_time).total_seconds()
        logger.info(f"Workflow completed in {duration:.2f}s")
        
        return self.results
    
    def _execute_parallel(
        self,
        plugin: PluginBase,
        inputs: Dict[str, Any],
        parallel_config: Dict[str, Any],
    ) -> List[Any]:
        """Execute plugin with parallel backend."""
        from ..parallel.executor import ParallelExecutor
        
        backend = self.config.parallel.backend
        executor = ParallelExecutor(backend, self.config.parallel)
        
        # Get items to process
        items = inputs.get("compounds", [])
        if not items:
            return [plugin.run(inputs)]
        
        # Execute in parallel
        chunk_size = parallel_config.get("chunk_size", self.config.parallel.chunk_size)
        results = executor.map(plugin.run, items, chunk_size=chunk_size)
        
        return results
    
    def finalize(self) -> Dict[str, Any]:
        """
        Finalize workflow: aggregate results, generate reports.
        
        Returns:
            Final summary dict.
        """
        logger.info("Finalizing workflow...")
        
        summary = {
            "status": "ok",
            "start_time": self.start_time.isoformat() if self.start_time else None,
            "end_time": self.end_time.isoformat() if self.end_time else None,
            "duration_seconds": (
                (self.end_time - self.start_time).total_seconds()
                if self.start_time and self.end_time
                else None
            ),
            "plugins_executed": list(self.results.keys()),
            "results": self.results,
        }
        
        # Write summary
        summary_path = self.config.results_dir / "workflow_summary.json"
        import json
        with open(summary_path, "w") as f:
            json.dump(summary, f, indent=2, default=str)
        
        logger.info(f"Summary written to {summary_path}")
        return summary
