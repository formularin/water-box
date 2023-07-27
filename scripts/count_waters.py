"""Runs the count_waters fortran script for spherical volumes of multiple radii and conducts various analyses.
Written by Arin Khare
"""

import math
import shlex
import subprocess

import numpy as np
import matplotlib.pyplot as plt


def main():
    # Plot all histograms
    all_data = []
    for i in range(10, 100, 5):
        r = str(round(i / 100, 2))
        if i % 10 == 0:
            r = r + "0"

        subprocess.run(shlex.split(f"bin/count_waters -r {r} -o count_waters/{r}.xvg"), capture_output=True)
        
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
    plt.gcf().set_size_inches(15, 5)
    plt.savefig(f"count_waters/all_hists.png", dpi=100)
    plt.clf()

    # Add zeros to the end of all_data so it can become a matrix
    # max_len = max([len(data[1]) for data in all_data])
    # for i, data in enumerate(all_data):
    #     new_data = np.zeros(shape=(2, max_len))
    #     new_data[1][:len(data[1])] = data[1]
    #     new_data[0] = np.arange(max_len)
    #     all_data[i] = new_data

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
    plt.gcf().set_size_inches(10, 10)
    plt.savefig("count_waters/boxplot.png", dpi=100)


if __name__ == "__main__":
    main()