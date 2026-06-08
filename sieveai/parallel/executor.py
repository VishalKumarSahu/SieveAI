"""Parallel execution executor - local and HPC backends."""

from typing import List, Any, Callable
from multiprocessing import Pool, cpu_count
import logging

logger = logging.getLogger("sieveai.parallel")


class ParallelExecutor:
    """
    Unified parallel execution interface.
    
    Backends:
        - local: multiprocessing.Pool
        - slurm: SLURM job arrays
        - pbs: PBS/Torque job arrays
        - dask: Dask distributed (optional)
    
    Usage:
        executor = ParallelExecutor("local", config)
        results = executor.map(plugin.run, compounds, chunk_size=10)
    """
    
    def __init__(self, backend: str = "local", config: dict = None):
        self.backend = backend
        self.config = config or {}
    
    def map(
        self,
        func: Callable,
        items: List[Any],
        chunk_size: int = 10,
    ) -> List[Any]:
        """
        Execute function over items in parallel.
        
        Args:
            func: Function to call on each item.
            items: List of inputs.
            chunk_size: Items per job (for batching).
        
        Returns:
            List of results.
        """
        if self.backend == "local":
            return self._local_map(func, items, chunk_size)
        elif self.backend == "slurm":
            return self._slurm_map(func, items, chunk_size)
        elif self.backend == "pbs":
            return self._pbs_map(func, items, chunk_size)
        else:
            logger.warning(f"Unknown backend '{self.backend}', using local")
            return self._local_map(func, items, chunk_size)
    
    def _local_map(
        self,
        func: Callable,
        items: List[Any],
        chunk_size: int,
    ) -> List[Any]:
        """Execute using multiprocessing.Pool."""
        max_workers = self.config.get("max_workers", cpu_count())
        
        logger.info(
            f"Running {len(items)} tasks on {max_workers} workers "
            f"(chunk_size={chunk_size})"
        )
        
        with Pool(processes=max_workers) as pool:
            results = pool.map(func, items, chunksize=chunk_size)
        
        return results
    
    def _slurm_map(
        self,
        func: Callable,
        items: List[Any],
        chunk_size: int,
    ) -> List[Any]:
        """
        Execute using SLURM job array.
        
        Strategy:
            1. Write items to file (one per line)
            2. Generate SLURM array script
            3. Submit job array
            4. Wait for completion
            5. Collect results from output files
        """
        from ..deploy.hpc_slurm import SlurmExecutor
        
        executor = SlurmExecutor(self.config.get("slurm", {}))
        return executor.run(func, items, chunk_size)
    
    def _pbs_map(
        self,
        func: Callable,
        items: List[Any],
        chunk_size: int,
    ) -> List[Any]:
        """Execute using PBS job array."""
        from ..deploy.hpc_pbs import PbsExecutor
        
        executor = PbsExecutor(self.config.get("pbs", {}))
        return executor.run(func, items, chunk_size)
