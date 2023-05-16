import math
import os
import sys
from math import cos, pi, sin
from operator import itemgetter
import numpy as np

import matplotlib.pyplot as plt
from graphviz import Graph
import numpy as np

# inp = "../../Conferences and Papers/2023 CIBCB/AAMatcher/AAMOut/"
outp = "./AAMTestFigs/"
finame = "./data.txt"
samps = 200 # 2 for each 100 tests
precision = 6
col_width = 8 + precision
popsize = 50

def get_data(filename: str):
    vals = [[[] for _ in range(popsize)] for _ in range(popsize)]
    fits = []
    with open(filename) as f:
        lines = f.readlines()
        holder = []
        mom = -1
        dad = -1
        for line in lines:
            if line.__contains__("Fitness"):
                started = False
                pass
            elif line.__contains__("Idxs"):
                line = line.rstrip()
                line = line.split(": ")[1].split("\t")
                mom = int(line[0])
                dad = int(line[1])
                pass
            elif line.__contains__("Fits"):
                pass
            elif line == "\n":
                vals[mom][dad] = holder
                vals[dad][mom] = holder
                holder = []
                pass
            elif line.__contains__("Crossover"):
                started = True
                pass
            elif started:
                line = line.rstrip()
                holder.append(int(line))
                pass
            else:
                line = line.rstrip()
                fits.append(int(line))
                pass
            pass
        pass

    return vals, fits


def make_boxplot(data: [], fits: [], parent_diff: [], first_parent: int):
    plt.style.use("seaborn-v0_8")
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)

    f = plt.figure()
    f.set_figheight(4.5)
    f.set_figwidth(10)
    plot = f.add_subplot(111)

    if first_parent > -1:
        plot.scatter([i + 1 for i in range(50)], fits, marker='x', color='r', zorder=10)
        bp = plot.boxplot(data, patch_artist=True,  zorder=5)
        plot.bar([i + 1 for i in range(50)], parent_diff, color='lime', zorder=2)

        f.suptitle("Boxplots of 100 Crossovers for Parent " + str(first_parent + 1), fontsize=14)
        plot.set_xlabel("with Parent", fontsize=12)
        plot.set_ylabel("Fitness of Children", fontsize=12)
        pass
    else:
        plot.scatter([i + 1 for i in range(50)], fits, marker='x', color='r', zorder=10)
        bp = plot.boxplot(data, patch_artist=True, zorder=5)

        f.suptitle("Boxplots of All Children", fontsize=14)
        plot.set_xlabel("For Parent", fontsize=12)
        plot.set_ylabel("Fitness of Children", fontsize=12)
        pass

    f.tight_layout()
    f.savefig(outp + "Crossover_" + str(first_parent + 1) + "_boxplot.png", dpi=300)
    plt.close()
    pass


def main():
    data, fits = get_data(finame)
    all_for = []
    all_parent_means = [[] for _ in range(50)]

    for mom in range(50):
        one_parent_means = []
        for dad in range(50):
            if mom == dad:
                one_parent_means.append(0)
            else:
                parent_mean = (fits[mom] + fits[dad])/2
                child_mean = np.mean(data[mom][dad])
                one_parent_means.append(child_mean - parent_mean)
                pass
            pass
        all_parent_means[mom] = one_parent_means
    pass

    for idx, dat in enumerate(data):
        make_boxplot(dat, fits, all_parent_means[idx], idx)
        holder = []
        for vals in dat:
            holder.extend(vals)
            pass
        all_for.append(holder)
        pass

    make_boxplot(all_for, fits, [0 for _ in range(50)], -1)

    # print(data)
    print("DONE")
    pass


main()