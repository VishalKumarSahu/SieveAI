"""SieveAI CLI - typer-based command line interface."""

import typer
from pathlib import Path
from typing import Optional

app = typer.Typer(help="SieveAI - Automated Drug Discovery Pipeline")


@app.command()
def dock(
    config: Optional[Path] = typer.Option(
        None, "-c", "--config",
        help="Configuration file (default: sieveai.config.toml)"
    ),
    receptors: Optional[Path] = typer.Option(
        None, "-r", "--receptors",
        help="Receptors directory"
    ),
    ligands: Optional[Path] = typer.Option(
        None, "-l", "--ligands",
        help="Ligands directory"
    ),
    output: Optional[Path] = typer.Option(
        None, "-o", "--output",
        help="Output directory"
    ),
    workers: int = typer.Option(
        4, "-w", "--workers",
        help="Number of parallel workers"
    ),
    dry_run: bool = typer.Option(
        False, "--dry-run",
        help="Preview workflow without executing"
    ),
):
    """Run docking pipeline."""
    from sieveai.core.engine import SieveEngine
    from sieveai.core.config import ConfigManager
    from UtilityLib.lib.path import EntityPath
    
    base_dir = EntityPath.cwd()
    
    # Load config
    cfg_mgr = ConfigManager()
    if config:
        cfg_mgr.load(EntityPath(config))
    
    cfg_mgr.resolve_paths(base_dir)
    cfg = cfg_mgr.config
    
    # Override from CLI
    if receptors:
        cfg.paths.receptors_dir = EntityPath(receptors)
    if ligands:
        cfg.paths.ligands_dir = EntityPath(ligands)
    if output:
        cfg.paths.results_dir = EntityPath(output)
    cfg.parallel.max_workers = workers
    
    # Initialize engine
    engine = SieveEngine(cfg)
    engine.initialize()
    
    # Validate
    validation = engine.validate()
    failed = [uid for uid, res in validation.items() if not res["valid"]]
    if failed:
        typer.echo(f"✗ Validation failed for: {', '.join(failed)}")
        raise typer.Exit(1)
    
    # Plan workflow
    workflow_order = engine.plan()
    
    if dry_run:
        typer.echo("Dry run - workflow order:")
        for uid in workflow_order:
            typer.echo(f"  → {uid}")
        return
    
    # Execute
    typer.echo("Starting workflow...")
    results = engine.execute()
    summary = engine.finalize()
    
    typer.echo(f"✓ Completed in {summary['duration_seconds']:.2f}s")
    typer.echo(f"Results: {summary['results']}")


@app.command()
def status(
    checkpoint: Optional[Path] = typer.Option(
        None, "-c", "--checkpoint",
        help="Checkpoint file"
    ),
):
    """Show workflow progress."""
    typer.echo("Workflow status: [placeholder - checkpoint system TBD]")


@app.command()
def resume(
    checkpoint: Optional[Path] = typer.Option(
        None, "-c", "--checkpoint",
        help="Checkpoint file"
    ),
):
    """Resume interrupted workflow."""
    typer.echo("Resume: [placeholder - checkpoint system TBD]")


@app.command()
def doctor():
    """Check system health and dependencies."""
    import shutil
    import sys
    
    typer.echo("SieveAI System Check\n")
    
    # Python
    typer.echo(f"✓ Python: {sys.version.split()[0]}")
    
    # Executables
    executables = ["vina", "prepare_receptor", "prepare_ligand", "chimerax", "obabel"]
    for exe in executables:
        path = shutil.which(exe)
        if path:
            typer.echo(f"✓ {exe}: {path}")
        else:
            typer.echo(f"✗ {exe}: NOT FOUND")
    
    # Python packages
    packages = ["biopython", "rdkit", "pandas", "numpy"]
    for pkg in packages:
        try:
            __import__(pkg.replace("-", "_"))
            typer.echo(f"✓ {pkg}: installed")
        except ImportError:
            typer.echo(f"✗ {pkg}: NOT INSTALLED")
    
    # Environment
    from sieveai.core.environment import EnvironmentDetector
    env = EnvironmentDetector.detect()
    typer.echo(f"\nEnvironment: {env.value}")
    
    if env.is_hpc():
        info = EnvironmentDetector.get_info()
        typer.echo(f"Job ID: {info['environment'].get('job_id', 'N/A')}")


@app.command()
def install(
    package: str = typer.Argument(..., help="Package to install"),
):
    """Install dependencies (guided)."""
    from sieveai.deploy.installers import install_package
    
    typer.echo(f"Installing {package}...")
    success = install_package(package)
    
    if success:
        typer.echo(f"✓ {package} installed successfully")
    else:
        typer.echo(f"✗ Failed to install {package}")
        typer.echo("Manual installation instructions:")
        typer.echo("  - Ubuntu: sudo apt install <package>")
        typer.echo("  - macOS: brew install <package>")
        typer.echo("  - conda: conda install -c conda-forge <package>")


@app.command()
def plugins():
    """List available plugins."""
    from sieveai.core.registry import PluginRegistry
    from UtilityLib.lib.path import EntityPath
    
    plugins_dir = EntityPath(__file__).parent.parent / "plugins"
    registry = PluginRegistry(plugins_dir)
    plugins = registry.discover()
    
    typer.echo("Available plugins:\n")
    for uid, info in plugins.items():
        typer.echo(f"  {uid} v{info['version']}")
        typer.echo(f"    {info['name']}")
        typer.echo(f"    Assignments: {', '.join(info['assignments'])}")
        typer.echo(f"    Requires: {', '.join(info['requires']) or 'none'}")
        typer.echo()


if __name__ == "__main__":
    app()
