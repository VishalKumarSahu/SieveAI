"""HPC SLURM executor."""

import subprocess
import tempfile
import os
from pathlib import Path
from typing import List, Any, Callable
import logging

logger = logging.getLogger("sieveai.hpc.slurm")


class SlurmExecutor:
    """Execute tasks using SLURM job arrays."""
    
    def __init__(self, config: dict = None):
        self.config = config or {}
        self.array_size = config.get("array_size", 1000)
        self.concurrent = config.get("concurrent", 50)
        self.cpus_per_task = config.get("cpus_per_task", 4)
        self.mem = config.get("mem", "8G")
        self.time = config.get("time", "24:00:00")
    
    def run(
        self,
        func: Callable,
        items: List[Any],
        chunk_size: int,
    ) -> List[Any]:
        """
        Submit SLURM job array.
        
        Strategy:
            1. Write items to input file
            2. Generate job script
            3. Submit array job
            4. Wait for completion
            5. Collect results
        
        Note: For now, this is a placeholder. Full implementation
        requires serialization of func and result collection.
        """
        logger.info(f"SLURM executor: {len(items)} items, chunk_size={chunk_size}")
        
        # Placeholder - for now fall back to local
        # Full implementation would:
        # - Serialize items to file
        # - Generate .sh script with array
        # - sbatch --array=0-N%concurrent
        # - Wait for completion
        # - Parse output files
        
        logger.warning("SLURM executor not fully implemented, using local fallback")
        
        from multiprocessing import Pool
        with Pool() as pool:
            results = pool.map(func, items, chunksize=chunk_size)
        
        return results
    
    def generate_script(self, items_file: Path, output_dir: Path) -> str:
        """Generate SLURM job array script."""
        script = f"""#!/bin/bash
#SBATCH --job-name=sieveai
#SBATCH --array=0-{len(items_file)}%{self.concurrent}
#SBATCH --cpus-per-task={self.cpus_per_task}
#SBATCH --mem={self.mem}
#SBATCH --time={self.time}
#SBATCH --output={output_dir}/slurm-%A_%a.out
#SBATCH --error={output_dir}/slurm-%A_%a.err

# Get task from array index
ITEM=$(sed -n "$(($SLURM_ARRAY_TASK_ID + 1))p" {items_file})

# Run task
python -c "from sieveai.cli.main import run_item; run_item('$ITEM')"
"""
        return script
    
    def submit(self, script: str) -> str:
        """Submit SLURM job, return job ID."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sh', delete=False) as f:
            f.write(script)
            script_path = f.name
        
        os.chmod(script_path, 0o755)
        
        result = subprocess.run(
            ["sbatch", script_path],
            capture_output=True,
            text=True,
        )
        
        # Parse job ID from output
        # "Submitted batch job 12345"
        job_id = result.stdout.strip().split()[-1]
        
        os.unlink(script_path)
        return job_id
