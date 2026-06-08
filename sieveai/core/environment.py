"""Environment detection: local, HPC (SLURM/PBS), or cloud."""

import os
from pathlib import Path
from typing import Literal
from enum import Enum


class EnvironmentType(str, Enum):
    """Supported deployment environments."""
    LOCAL = "local"
    HPC_SLURM = "hpc-slurm"
    HPC_PBS = "hpc-pbs"
    AWS_BATCH = "aws-batch"
    DOCKER = "docker"
    SINGULARITY = "singularity"


class EnvironmentDetector:
    """Detects the current runtime environment."""
    
    @staticmethod
    def detect() -> EnvironmentType:
        """
        Detect runtime environment by checking system markers.
        
        Priority:
            1. Container (Singularity > Docker)
            2. HPC (SLURM > PBS)
            3. Cloud (AWS Batch)
            4. Local (default)
        """
        # Check for container environments
        if Path("/.singularity.d").exists():
            return EnvironmentType.SINGULARITY
        
        if Path("/.dockerenv").exists() or os.path.exists("/etc/containers"):
            return EnvironmentType.DOCKER
        
        # Check for HPC schedulers
        if os.path.exists("/etc/slurm_version") or os.environ.get("SLURM_JOB_ID"):
            return EnvironmentType.HPC_SLURM
        
        if os.path.exists("/var/spool/pbs") or os.environ.get("PBS_JOBID"):
            return EnvironmentType.HPC_PBS
        
        # Check for cloud batch systems
        if os.environ.get("AWS_BATCH_JOB_ID"):
            return EnvironmentType.AWS_BATCH
        
        # Default to local
        return EnvironmentType.LOCAL
    
    @staticmethod
    def get_info() -> dict:
        """Get detailed environment information."""
        env_type = EnvironmentDetector.detect()
        
        info = {
            "type": env_type.value,
            "hostname": os.uname().nodename,
            "cpu_count": os.cpu_count() or 1,
            "environment": {},
        }
        
        if env_type == EnvironmentType.HPC_SLURM:
            info["environment"] = {
                "job_id": os.environ.get("SLURM_JOB_ID"),
                "array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
                "cpus_per_task": os.environ.get("SLURM_CPUS_PER_TASK"),
                "mem_per_task": os.environ.get("SLURM_MEM_PER_CPU"),
                "nodes": os.environ.get("SLURM_JOB_NUM_NODES"),
            }
        
        elif env_type == EnvironmentType.HPC_PBS:
            info["environment"] = {
                "job_id": os.environ.get("PBS_JOBID"),
                "queue": os.environ.get("PBS_QUEUE"),
                "nodes": os.environ.get("PBS_NUM_NODES"),
                "ppn": os.environ.get("PBS_NUM_PPN"),
            }
        
        elif env_type == EnvironmentType.AWS_BATCH:
            info["environment"] = {
                "job_id": os.environ.get("AWS_BATCH_JOB_ID"),
                "job_definition": os.environ.get("AWS_BATCH_JOB_DEFINITION"),
            }
        
        return info
    
    @staticmethod
    def is_hpc() -> bool:
        """Check if running on HPC cluster."""
        return EnvironmentDetector.detect() in (
            EnvironmentType.HPC_SLURM,
            EnvironmentType.HPC_PBS,
        )
    
    @staticmethod
    def is_container() -> bool:
        """Check if running in container."""
        return EnvironmentDetector.detect() in (
            EnvironmentType.DOCKER,
            EnvironmentType.SINGULARITY,
        )
