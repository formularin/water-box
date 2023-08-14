# water-box

Fortran and Python scripts that tell you everything there is to know about water (not even close).

## Compilation

* Install [libgmxfort](https://github.com/formularin/libgmxfort) (depends on [libxdrfile](https://github.com/wesbarnett/libxdrfile))
* If on macOS, you may have to install pkg-config (most easily done from [homebrew](https://formulae.brew.sh/formula/pkg-config))
* Run the command: `make` from the project home directory

## Usage

Options are listed in the source file for each program.

Notes about default options:
 - Output goes to a directory called "results." It is assumed that this already exists.
 - Scripts use the GROMACS output from `simulations/water/`.

### Fortran programs

- Radial distribution function: `bin/rdf`
- Diffusion constant: `bin/diffusion`
- P(N) in a spherical volume: `bin/count_waters`
- Number of hydrogen bonds: `bin/hbonds`

### Python programs

Install dependencies using this command:

```bash
pip install -r requirements.txt
```

- Radial distribution function: `python3 scripts/rdf.py`
- Making plots for P(N) in a spherical volume: `python3 scripts/count_waters.py`

### Running with other GROMACS simulations

- All fortran scripts take the options `-xtc` and `-ndx` to specify GROMACS trajectory and index files.
- Index files must contain "OW", "HW1", and "HW2" groups.
    - Make an index file using [`gmx make_ndx`](https://manual.gromacs.org/current/onlinehelp/gmx-make_ndx.html), passing in the .gro file of the production run.
    - Use the commands `a OW`, `a HW1`, and `a HW2` to create the groups.
<br />
- This will work with any three-point water model, like SPC or TIP3P
- `simulations/water/water_box.sh`, as well as `simulations/water/mdp/` are adapted from the [GASERI water tutorial](https://gaseri.org/en/tutorials/gromacs/1-tip4pew-water/)