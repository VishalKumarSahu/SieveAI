# SieveAI
An automated drug discovery pipeline

## Features
* Automated blind docking and rescoring
* Automated result parsing and aggregation
* Basic version for AutoDock VINA 1.2.3
* Customisable to include local, api, or web based tools

## Quick installation and usage
* `pip install SieveAI`
* `sieveai -h`

## Installing basic requirements for basic docking

### Requirements for AutoDock VINA
1. Python
    - [Download](https://www.python.org/downloads/) and install python
    - Preferred versions 3.8 - 3.10

2. SieveAI > 0.4
`pip install sieveai`

3. AutoDock VINA 1.2.3
  * Step 1: Navigate to [vina binary on GitHub](https://github.com/ccsb-scripps/AutoDock-Vina/releases/tag/v1.2.3)
  * Ubuntu 20 LTS
    - Step 2: Download `vina_1.2.3_linux_x86_64` (for Linux)
    - Step 3: Rename downloaded file to `vina` and move to `/opt/AutoDock/` in Ubuntu
  * Windows
    - Step 2: Download `vina_1.2.3_windows_x86_64.exe` (for windows)
    - Step 3: Rename downloaded file to `vina` and move it program directory. Add directory to environmental variables.

4. AutoDock Tools/MGL Tools
  - Download [ADFR Suit for Linux/Windows](https://ccsb.scripps.edu/adfr/downloads/) and install as per given instructions
  - Add executable path to bashrc in Ubuntu or add directory to environmental variables in Windows

> Make sure that `vina`, `prepare_ligand`, `mk_prepare_ligand` and `prepare_receptor` commands available through commandline

5. ChimeraX
  - Step 1: Download Chimerax from [download page](https://www.rbvi.ucsf.edu/chimerax/download.html)
  - Step 2: `sudo apt-get install ./ucsf-chimerax_1.3ubuntu20.04_amd64.deb`
  - Confirm ChimeraX availble through commandline `chimerax`

6. Openbabel (Optional)
  - Ubuntu 20 LTS
    * STEP 1: `sudo apt-get install openbabel`
    * Confirm OpenBabel availble through commandline `obabel`
  - Windows
    * Get OpenBabel download from [here](https://openbabel.org/docs/dev/Installation/install.html)
    * Add `obabel` path to environment PATH variables

7. FreeSASA (Optional)
  * Follow instructions at https://freesasa.github.io/
  * Download https://freesasa.github.io/freesasa-2.0.3.tar.gz and extract
  * ```sh
    ./configure
    make
    sudo make install
    freesasa -h```

## 0.5.20230606
* Bug Fixes

## 0.1.20210630
* Version 0.1

## Copyright
&copy; Vishal Kumar Sahu, 2023
