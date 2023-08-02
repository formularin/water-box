"""Computes the compressibilty of water given energy data from a simulation.
Written by Arin Khare

Usage: python3 compressibility.py -f <energy_file> -o <output_file>

- energy_file must be an XVG output from the `gmx energy` command, with data for pressure and volume.
- the output file will be an XVG containing data sorted by volume, with columns for log volume and pressure.
"""

import sys

import numpy as np


CHUNK_SIZE = 100
TEMP_INTERVAL = 20


if __name__ == "__main__":
    energy_file = "simulations/water/prd_energy.xvg"
    output_file = "results/compressibility.xvg"

    for i in range(1, len(sys.argv)):
        if sys.argv[i] == "-f":
            energy_file = sys.argv[i + 1]
        elif sys.argv[i] == "-o":
            output_file = sys.argv[i + 1]

    time, pres, vol = np.loadtxt(energy_file, comments="#", unpack=True)

    log_vol = np.log(vol)
    A = np.vstack([log_vol, np.ones(len(log_vol))]).T
    slope, y_int = np.linalg.lstsq(A, pres, rcond=None)[0]
    print(f"{slope=}\t{y_int=}")
    print(f"compressibility={-1/slope}")

    np.savetxt(output_file, np.array([log_vol, pres]).T, delimiter="   ", fmt="%20.7f")