# SieveAI v1.0 - Automated Drug Discovery Pipeline

Production-ready, modular drug discovery pipeline with plugin architecture.

## Features

- **Single-file plugins** - Drop-in/drop-out, YAML frontmatter metadata
- **Config-driven workflow** - TOML configuration, auto dependency resolution
- **Parallel execution** - Local multiprocessing, SLURM/PBS HPC support
- **Checkpoint/resume** - SQLite-based progress tracking
- **Production-grade** - 2000+ dockings/day on i5/12GB (tested)

## Installation

```bash
pip install -e .
```

## Quick Start

```bash
# Check system
sieveai doctor

# Run docking
sieveai dock -r receptors/ -l ligands/ -w 8

# List plugins
sieveai plugins
```

## Configuration

Create `sieveai.config.toml`:

```toml
[environment]
type = "local"

[paths]
receptors_dir = "receptors"
ligands_dir = "ligands"

[parallel]
backend = "local"
max_workers = 8
chunk_size = 10

[plugins.vina]
enabled = true
exhaustiveness = 16
num_modes = 10
```

## Plugin Development

Create `plugins/my_plugin.py`:

```python
"""
---
uid: my_plugin
name: My Plugin
version: 1.0.0
assignments: [docking]
requires: [structuresync]
provides: [results]
executables: [my_exe]
python_packages: [mypackage]
parallel:
  strategy: compound
  chunk_size: 10
---
"""

from sieveai.plugins.base import PluginBase

class MyPlugin(PluginBase):
    def validate(self):
        return True
    
    def init(self, **kwargs):
        pass
    
    def run(self, item):
        return {"status": "ok"}
```

## HPC Usage

```bash
# SLURM
sieveai dock --config sieveai-slurm.config.toml

# PBS
sieveai dock --config sieveai-pbs.config.toml
```

## License

MIT License - See LICENSE file.

## Citation

Basu, S., Sahu, V.K., Ranjan, A., et al. SieveAI: Development of an Automated extensible and customisable drug discovery pipeline. bioRxiv (2023).
