import math
import os
import sys
from math import cos, pi, sin
from operator import itemgetter
import numpy as np

import copy
import matplotlib.pyplot as plt
from graphviz import Graph
import numpy as np
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple
from matplotlib.patches import Patch
import matplotlib.patches as mpatches

# inp = "../../Conferences and Papers/2023 CIBCB/AAMatcher/AAMOut/"
inp = "./AAMOut/"
outp = "./AAMTestFigs/"
finame1 = "crossover01_start.dat"
finame2 = "crossover01_end.dat"
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


def make_boxplot(data: [], fits: [], parent_diff: [], first_parent: int, out_path, run, start):
    data = copy.deepcopy(data)
    plt.style.use("seaborn-v0_8")
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)

    f = plt.figure()
    f.set_figheight(4.5)
    f.set_figwidth(10)
    plot = f.add_subplot(111)

    if first_parent > -1:
        out_path = out_path + "Run" + str(run).zfill(2) + "_Start" if start else "_End"
        if not os.path.exists(out_path):
            os.makedirs(out_path)
            pass
        out_path += "/"
        sp = plot.scatter([i + 1 for i in range(50)], fits, marker='2', color='#DB57B2', zorder=10, s=100, linewidth=1)

        xs = [i+1 for i in range(popsize)]
        for idx in range(popsize):
            if data[idx] == []:
                to_del = idx
                pass
            pass
        del data[to_del]
        del xs[to_del]
        vp = plot.violinplot(data, xs, showmedians=True, widths=0.85)
        for pc in vp["bodies"]:
            pc.set_facecolor("#5770DB")
            pc.set_linewidth(1)
            pc.set_edgecolor("black")
            pass
        for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
            vc = vp[partname]
            vc.set_linewidth(1)
            vc.set_alpha(1)
            vc.set_color("#5770DB")
            pass

        pos_bars, neg_bars, pos_xs, neg_xs = [], [], [], []
        for idx in range(popsize):
            if parent_diff[idx] >= 0:
                pos_xs.append(idx + 1)
                pos_bars.append(parent_diff[idx])
            else:
                neg_xs.append(idx + 1)
                neg_bars.append(parent_diff[idx])
                pass
        posb = plot.bar(pos_xs, pos_bars, color="#57DB80", zorder=0.9)
        negb = plot.bar(neg_xs, neg_bars, color="#DB5F57", zorder=0.9)

        labels = ["Parent's Fitness", "Children's Fitness", "Mean of Children's Fitness \u2212 Mean of Parents' Fitness"]
        patches = []
        patches.append(sp)
        patches.append(Patch(color="#5770DB", label="Children's Fitness"))
        patches.append([mpatches.Patch(color="#57DB80", label=labels[2]), mpatches.Patch(color="#DB5F57", label=labels[2])])
        plot.legend(handles=patches, labels=labels, facecolor='white', frameon='true', fontsize=12, framealpha=0.75, loc='upper center', ncol=3, borderaxespad=0.1, handler_map = {list: HandlerTuple(None)})
        plot.grid(visible="True", axis="y", which='minor', color="white", linewidth=0.5)
        plt.minorticks_on()

        if start:
            f.suptitle("Fitness of Children from 100 Two-Point Crossover Operations using Parent " + str(first_parent + 1) + " Before Evolution", fontsize=14)
        else:
            f.suptitle("Fitness of Children from 100 Two-Point Crossover Operations using Parent " + str(first_parent + 1) + " After Evolution", fontsize=14)
            pass
        plot.set_xlabel("And Parent", fontsize=12)

        plt.ylim(-5, 45)
        plt.xlim(0, 51)
        plt.xticks([i + 1 for i in range(popsize)])
        f.subplots_adjust(bottom=0.1, top=0.93, left=0.03, right=0.99)
        f.savefig(out_path + "Crossover_Run" + str(run).zfill(2) + "_P" + str(first_parent + 1) + ".png", dpi=500)
        plt.close()
        pass
    else:
        sp = plot.scatter([i + 1 for i in range(50)], fits, marker='2', color='#DB57B2', zorder=10, s=100, linewidth=1)

        xs = [i + 1 for i in range(popsize)]
        vp = plot.violinplot(data, xs, showmedians=True, widths=0.85)
        for pc in vp["bodies"]:
            pc.set_facecolor("#5770DB")
            pc.set_linewidth(1)
            pc.set_edgecolor("black")
            pass
        for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
            vc = vp[partname]
            vc.set_linewidth(1)
            vc.set_alpha(1)
            vc.set_color("#5770DB")
            pass

        pos_bars, neg_bars, pos_xs, neg_xs = [], [], [], []
        for idx in range(popsize):
            if parent_diff[idx] >= 0:
                pos_xs.append(idx + 1)
                pos_bars.append(parent_diff[idx])
            else:
                neg_xs.append(idx + 1)
                neg_bars.append(parent_diff[idx])
                pass
        posb = plot.bar(pos_xs, pos_bars, color="#57DB80", zorder=0.9)
        negb = plot.bar(neg_xs, neg_bars, color="#DB5F57", zorder=0.9)

        f.suptitle("Boxplots of All Children", fontsize=14)
        plot.set_xlabel("For Parent", fontsize=12)
        plot.set_ylabel("Fitness of Children", fontsize=12)
        labels = ["Parent's Fitness", "Children's Fitness",
                  "Mean of Children's Fitness \u2212 Parent's Fitness"]
        patches = []
        patches.append(sp)
        patches.append(Patch(color="#5770DB", label="Children's Fitness"))
        patches.append([mpatches.Patch(color="#57DB80", label=labels[2]), mpatches.Patch(color="#DB5F57", label=labels[2])])
        plot.legend(handles=patches, labels=labels, facecolor='white', frameon='true', fontsize=12, framealpha=0.75,
                    loc='upper center', ncol=3, borderaxespad=0.1, handler_map={list: HandlerTuple(None)})
        plot.grid(visible="True", axis="y", which='minor', color="white", linewidth=0.5)
        plt.minorticks_on()

        plt.ylim(-5, 45)
        plt.xlim(0, 51)
        plt.xticks([i + 1 for i in range(popsize)])
        f.subplots_adjust(bottom=0.1, top=0.93, left=0.03, right=0.99)
        if start:
            f.savefig(out_path + "Crossover_Run" + str(run) + "_StartSummary.png", dpi=500)
        else:
            f.savefig(out_path + "Crossover_Run" + str(run) + "_EndSummary.png", dpi=500)
            pass
        plt.close()
        pass


    pass


def main():
    folder_names = os.listdir(inp)

    for fold in folder_names:
        out_path = outp + fold
        if not os.path.exists(out_path):
            os.makedirs(out_path)
            pass
        out_path += "/"

        data, fits = get_data(inp + fold + "/" + finame1)
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
            make_boxplot(dat, fits, all_parent_means[idx], idx, out_path, 1, True)
            holder = []
            for vals in dat:
                holder.extend(vals)
                pass
            all_for.append(holder)
            pass

        make_boxplot(all_for, fits, [np.mean(all_for[i]) - fits[i] for i in range(popsize)], -1, out_path, 1, True)
        pass

    print("DONE")
    pass


main()