"""
Computes the O-O radial distribution function of water
written by Arin Khare

usage: python3 rdf.py -xtc -gro -o

- Outputs an XVG file containing these columns:
    - The right edge of each bin, in nm
    - The frequency of distances for each bin (averaged over all frames)
    - The total volume for each bin
    - The number density corresponding to each bin
    - The normalized RDF value for each bin
"""

import itertools
import math
import sys
import time

from fast_histogram import histogram1d
import mdtraj as md
import numpy as np

import matplotlib.pyplot as plt

NUM_BINS = 1000


def get_rdf(traj):
    print("Calculating volumes... ", end="", flush=True)
    oxygens = traj.topology.select("symbol == O")
    max_dist = np.min(traj.unitcell_volumes) ** (1 / 3) / 2
    dr = max_dist / NUM_BINS
    # The value to divide by is the number of oxygens times the volume of the spherical shell,
    # since we are averaging over all oxygens. This will ensure we are actually calculating
    # number density, not something proportional to it.
    total_volumes = np.zeros(NUM_BINS)
    prev_sphere_vol = 0
    for i in range(NUM_BINS):
        sphere_vol = (4 / 3) * math.pi * (dr * (i + 1)) ** 3
        volume = sphere_vol - prev_sphere_vol
        prev_sphere_vol = sphere_vol
        total_volumes[i] = volume * len(oxygens)
    print("done")

    print("Calculating distances... ", end="", flush=True)
    o_pairs = np.array(list(itertools.combinations(oxygens, 2)))
    distances = md.compute_distances(traj, o_pairs)  # Follows minimum image convention
    distances = np.ravel(distances)
    print("done")

    print("Calculating frequencies... ", end="", flush=True)
    frequencies = histogram1d(distances, NUM_BINS, range=(0, max_dist + dr / 2))
    print("done")

    print("Normalizing distribution... ", end="", flush=True)
    avg_num_density = len(oxygens) / np.mean(traj.unitcell_volumes)
    # Number of frames for each radius that have unitcells large enough to contain that radius.
    num_frames = np.zeros(NUM_BINS)
    for x, y, z in traj.unitcell_lengths:
        diagonal = 0.5 * math.sqrt((x) ** 2 + (y) ** 2 + (z) ** 2)
        for i in range(NUM_BINS):
            r = dr * i
            if diagonal >= r:
                num_frames[i] += 1
            else:
                break
    # Multiply by 2 first because each distance needs to be counted twice.
    frequency_rates = (2 * frequencies) / num_frames
    rdf = frequency_rates / total_volumes
    normalized_rdf = rdf / avg_num_density
    print("done")

    return np.transpose(
        [np.arange(dr, max_dist + dr / 2, dr),
         frequencies,
         total_volumes,
         rdf,
         normalized_rdf]) 


if __name__ == "__main__":
    start_time = time.time()
    out_file = "results/rdf.xvg"
    traj_file = "simulations/water/spce.xtc"
    ndx_file = "simulations/water/spce.gro"

    for i in range(1, len(sys.argv) - 1):
        if sys.argv[i] == "-xtc":
            traj_file = sys.argv[i + 1]
        if sys.argv[i] == "-ndx":
            ndx_file = sys.argv[i + 1]
        if sys.argv[i] == "-o":
            out_file = sys.argv[i + 1]
    print("Reading trajectory file... ", end="", flush=True)
    traj = md.load(traj_file, top=ndx_file)
    print("done")
    rdf = get_rdf(traj)
    print(f"Writing to {out_file}... ", end="", flush=True)
    np.savetxt(out_file, rdf, delimiter="  ")
    print(f"done. RDF data saved to {out_file}")
    print(f"Program finished in {time.time() - start_time} seconds")