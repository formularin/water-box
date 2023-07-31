"""Runs the count_waters fortran script for spherical volumes of multiple radii and makes plots.
Written by Arin Khare
"""

import math
import os
import shlex
import subprocess

import numpy as np
import matplotlib.pyplot as plt


def main():
    try:
        os.mkdir("results/count_waters")
    except FileExistsError:
        pass

    # Plot all histograms
    all_data = []
    for i in range(10, 100, 5):
        r = str(round(i / 100, 2))
        if i % 10 == 0:
            r = r + "0"

        subprocess.run(shlex.split(f"bin/count_waters -r {r} -o results/count_waters/{r}.xvg"), capture_output=True)
        
        with open(f"count_waters/{r}.xvg", "r") as f:
            lines = f.read().strip().split("\n")
        avg_cell_vol = float(lines[0].split()[-1])
        num_waters = int(lines[1].split()[-1])
        data = [[float(val) for val in line.split()] for line in lines[2:]]
        data = np.transpose(data)
        all_data.append(data)

        # Remove beginning zeros for plotting
        for j in range(len(data[1])):
            if data[1][j] != 0:
                break
        data = data[:,j:]
        plt.plot(data[0], data[1], label=f"r={r}nm")
        plt.xlabel("Number of waters")
        plt.ylabel("P(N)")

    plt.legend(ncol=2)
    plt.yticks(np.arange(0, 0.9, 0.1), minor=True)
    plt.xticks(np.arange(0, 130, 10), minor=True)
    plt.tick_params(which='minor', length=0)
    plt.grid()
    plt.grid(which='minor', alpha=0.3)
    plt.gcf().set_size_inches(15, 5)
    plt.savefig(f"results/count_waters/all_hists.png", dpi=100)
    plt.clf()

    # Reconstruct raw counts for boxplot
    raw_counts = []
    for data in all_data:
        counts = []
        for n, pn in zip(data[0], data[1]):
            for i in range(int(pn * num_waters)):
                counts.append(n)
        raw_counts.append(counts)

    # Compare observed means with theoretical means assuming uniform waters
    theoretical = np.array([((4/3) * math.pi * (r/100)**3 / (avg_cell_vol)) * num_waters for r in range(10, 100, 5)])
    plt.plot(np.arange(0.1, 1.0, 0.05), theoretical, label="Volume fraction * total number of waters")
    plt.boxplot(raw_counts, positions=np.arange(0.1, 1.0, 0.05), manage_ticks=False, widths=0.03)
    plt.xlabel("Radius (nm)")
    plt.ylabel("Number of waters")
    plt.legend()
    plt.xticks(np.arange(0, 1, 0.05), minor=True)
    plt.yticks(np.arange(0, 130, 10), minor=True)
    plt.tick_params(which='minor', length=0)
    plt.grid()
    plt.grid(which='minor', alpha=0.3)
    plt.gcf().set_size_inches(10, 10)
    plt.savefig("results/count_waters/boxplot.png", dpi=100)


if __name__ == "__main__":
    main()