"""
---
uid: vina
name: AutoDock VINA
version: 1.0.0
assignments: [docking]
requires: [structuresync]
provides: [docked_complexes]
executables: [vina, prepare_receptor, prepare_ligand]
python_packages: [biopython, rdkit]
parallel:
  strategy: compound
  chunk_size: 10
---
"""

from ..plugins.base import PluginBase
from UtilityLib.lib.path import EntityPath
from typing import Any, Dict, List
import shutil


class VinaPlugin(PluginBase):
    """AutoDock VINA docking plugin."""
    
    uid = "vina"
    name = "AutoDock VINA"
    version = "1.0.0"
    assignments = ["docking"]
    requires = ["structuresync"]
    provides = ["docked_complexes"]
    executables = ["vina", "prepare_receptor", "prepare_ligand"]
    python_packages = ["biopython", "rdkit"]
    parallel_strategy = "compound"
    chunk_size = 10
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.compounds = []
        self.receptor = None
        self.results = []
    
    def validate(self) -> bool:
        """Check VINA dependencies."""
        valid = True
        
        for exe in self.executables:
            if not self.check_executable(exe):
                self._log("error", f"Missing executable: {exe}")
                valid = False
        
        for pkg in self.python_packages:
            if not self._import(pkg):
                self._log("error", f"Missing Python package: {pkg}")
                valid = False
        
        return valid
    
    def init(self, compounds: List[dict] = None, receptor: EntityPath = None, **kwargs):
        """Initialize docking run."""
        self.compounds = compounds or []
        self.receptor = receptor
        self.docking_dir = (self.workdir / "docking" / "vina").validate()
        self.results = []
        
        self._log("info", f"Initialized VINA: {len(self.compounds)} compounds")
    
    def run(self, compound: Dict[str, Any]) -> Dict[str, Any]:
        """
        Dock single compound.
        
        Args:
            compound: Dict with id, path_pdbqt, etc.
        
        Returns:
            Result dict with status, score, output paths.
        """
        cid = compound.get("id", "unknown")
        ligand_pdbqt = compound.get("path_pdbqt")
        
        if not ligand_pdbqt:
            return {
                "compound_id": cid,
                "status": "failed",
                "error": "No ligand PDBQT file",
            }
        
        out_path = self.docking_dir / f"{cid}.docked.pdbqt"
        log_path = self.docking_dir / f"{cid}.log"
        
        # Build VINA command
        cmd = [
            self._which("vina"),
            "--receptor", str(self.receptor),
            "--ligand", str(ligand_pdbqt),
            "--out", str(out_path),
            "--log", str(log_path),
            "--exhaustiveness", str(self.config.get("exhaustiveness", 16)),
            "--num_modes", str(self.config.get("num_modes", 10)),
        ]
        
        # Add center/size if provided
        if "center" in self.config:
            c = self.config["center"]
            cmd.extend(["--center_x", str(c.get("x", 0))])
            cmd.extend(["--center_y", str(c.get("y", 0))])
            cmd.extend(["--center_z", str(c.get("z", 0))])
        
        if "size" in self.config:
            s = self.config["size"]
            cmd.extend(["--size_x", str(s.get("x", 20))])
            cmd.extend(["--size_y", str(s.get("y", 20))])
            cmd.extend(["--size_z", str(s.get("z", 20))])
        
        # Run VINA
        result = self._run_command(cmd, timeout=3600)
        
        if result.returncode == 0:
            score = self._parse_score(log_path)
            return {
                "compound_id": cid,
                "status": "ok",
                "output": str(out_path),
                "log": str(log_path),
                "score": score,
                "modes": self._parse_modes(log_path),
            }
        else:
            return {
                "compound_id": cid,
                "status": "failed",
                "error": result.stderr,
                "log": str(log_path),
            }
    
    def finalize(self, results: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Aggregate docking results."""
        import pandas as pd
        
        df = pd.DataFrame(results)
        
        # Write Excel
        excel_path = self.workdir / "Results.vina.xlsx"
        df.to_excel(excel_path, index=False)
        
        # Write summary
        summary = {
            "status": "ok",
            "count": len(results),
            "success": sum(1 for r in results if r.get("status") == "ok"),
            "failed": sum(1 for r in results if r.get("status") == "failed"),
            "output_excel": str(excel_path),
            "best_score": min(
                (r.get("score", 999) for r in results if r.get("status") == "ok"),
                default=None,
            ),
        }
        
        self._log("info", f"VINA completed: {summary['success']}/{summary['count']} successful")
        
        return summary
    
    def _parse_score(self, log_path: EntityPath) -> float:
        """Parse VINA log for best affinity."""
        if not log_path.exists():
            return None
        
        content = log_path.read_text()
        for line in content.split("\n"):
            if "Affinity:" in line and line.strip().startswith("#"):
                try:
                    return float(line.split(":")[1].strip().split()[0])
                except (ValueError, IndexError):
                    continue
        
        return None
    
    def _parse_modes(self, log_path: EntityPath) -> int:
        """Parse number of modes from log."""
        if not log_path.exists():
            return 0
        
        content = log_path.read_text()
        for line in content.split("\n"):
            if "modes:" in line.lower():
                try:
                    return int(line.split(":")[1].strip().split()[0])
                except (ValueError, IndexError):
                    continue
        
        return 0
