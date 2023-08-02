# water-box

Trajectory files and scripts for analyzing a GROMACS simulation of SPC/E water.

## Compilation

* Install [libgmxfort](https://github.com/formularin/libgmxfort)
* Run the command: `make`

## Usage

Scripts output to a directory called "results" by default, and assume it already exists.

### Fortran programs

Options are listed in source files.

- Radial distribution function: `bin/rdf`
- Diffusion constant: `bin/diffusion`
- P(N) in a spherical volume: `bin/count_waters`

### Python programs

Install dependencies using this command:

```bash
pip install -r requirements.txt
```

- Radial distribution function: `python3 scripts/rdf.py`
- Making plots for P(N) in a spherical volume: `python3 scripts/count_waters.py`
