import math
import os
import sys
from math import cos, pi, sin
from operator import itemgetter

import matplotlib.pyplot as plt
from graphviz import Graph
import numpy as np

inp = "../../Conferences and Papers/2023 CIBCB/AAMatcher/AAMOut/"
outp = "./AAMFigs/"
finame = "exp.dat"
samps = 50
precision = 6
col_width = 8 + precision


def writeStat(data: [], out, lower_better):
    mean = float(np.mean(data))
    mean = round(mean, precision)
    std = float(np.std(data, ddof=0))
    std = round(std, precision)  # Population standard deviation
    diff = 1.96 * std / math.sqrt(30)  # 95% CI
    diff = round(diff, precision)
    if lower_better:
        maxima = float(min(data))
        pass
    else:
        maxima = float(max(data))
        pass
    maxima = round(maxima, precision)
    out.write(str(maxima).ljust(col_width))
    out.write(str(mean).ljust(col_width))
    out.write(str(std).ljust(col_width))
    out.write(str(diff).ljust(col_width))
    return mean, maxima


def make_table(many_data: [], best_run, exp_info: [], fname: str, minimizing: bool):
    with open(fname, "w") as f:
        f.write("EXP".ljust(col_width))
        f.write("Parameters".ljust(3 * col_width))
        f.write("Best Run".ljust(col_width))
        f.write("Best Fit".ljust(col_width))
        f.write("Mean".ljust(col_width))
        f.write("SD".ljust(col_width))
        f.write("95% CI".ljust(col_width))
        f.write("\n")
        for di, data in enumerate(many_data):
            f.write(str("EXP" + str(di + 1)).ljust(col_width))
            f.write(exp_info[di].ljust(3 * col_width))
            f.write(str("Run " + str(best_run[di])).ljust(col_width))
            # f.write(str("Run " + str(data[0][0])).ljust(col_width))
            writeStat(data, f, minimizing)
            f.write("\n")
            pass
        pass
    pass


def box_plot(bp, num_splits: int, split_info: []):
    for whisker in bp['whiskers']:
        whisker.set(color='#8B008B', linewidth=1)
        pass

    for cap in bp['caps']:
        cap.set(color='#8B008B', linewidth=1)
        pass

    for median in bp['medians']:
        median.set(color='#FF7F0E', linewidth=1)
        pass

    for flier in bp['fliers']:
        flier.set(marker='.', color='#e7298a', alpha=0.5, markersize=5)
        pass

    for info in split_info:
        for idx in info[0]:
            bp['boxes'][idx].set(facecolor=info[1][idx % len(info[1])])
            pass
        pass
    pass


def calc(data):
    mean = float(np.mean(data))
    std = float(np.std(data, ddof=0))
    diff = 1.96 * std / math.sqrt(30)  # 95% CI
    print(str(mean) + "+-" + str(diff))
    pass


def combine(dir: str, start: str, end: str):
    if os.path.exists(dir + start + end):
        os.remove(dir + start + end)
        pass
    out_file = open(dir + start + end, "a")
    for idx in range(1, 31):
        with open(dir + start + str(idx).zfill(2) + end, "r") as f:
            out_file.write(f.read())
            pass
        pass
    out_file.close()
    pass


def get_data(dir_path: str):
    fits = []
    SDAs = []
    seqs = []
    with open(dir_path + finame) as f:
        lines = f.readlines()
        next_SDA = False
        SDA = []
        for line in lines:
            if line.__contains__("best fitness"):
                line = line.rstrip()
                line = line.split(" ")
                fits.append(int(line[-1]))
                pass
            elif line.__contains__(str("Best Match")):
                line = line.rstrip()
                line = line.split(" ")
                seqs.append(line[-1])
                pass
            elif line.__contains__("SDA"):
                next_SDA = True
                pass
            elif line.__contains__("Fitness Values"):
                next_SDA = False
                SDAs.append(SDA)
                pass
            elif next_SDA:
                SDA.append(line)
                pass
            pass
        pass

    # run number, fitness, profileS, dnaS, edge count
    if len(fits) != samps:
        print("ERROR in fits: " + dir_path)
        pass
    if len(SDAs) != samps:
        print("ERROR in SDAs: " + dir_path)
        pass
    if len(seqs) != samps:
        print("ERROR in seqs: " + dir_path)
        pass

    data = [[i + 1, fits[i], SDAs[i], seqs[i]] for i in range(samps)]
    data.sort(key=itemgetter(1))  # Ascending
    data.reverse()
    return data


def str_to_list(orig: str):
    rtn = []
    for c in orig:
        rtn.append(int(c))
        pass
    return rtn

def print_best_info(path: str, info: str, dat: [], true_seq: []):
    with open(path, "w") as f:
        f.write(info + "\n")
        f.write("Run number: " + str(dat[0]) + "\n")
        f.write("With Fitness: " + str(dat[1]) + "\n")
        test_seq = str_to_list(dat[3])
        similar_str = ''
        for idx in range(len(true_seq)):
            if test_seq[idx] == true_seq[idx]:
                similar_str += 'X'
            else:
                similar_str += '-'
                pass
            pass
        f.write(int_to_DNA(test_seq) + '\n')
        f.write(similar_str + '\n')
        f.write(int_to_DNA(true_seq) + '\n')
        f.writelines(dat[2])
        pass
    pass


def gen_sequences(path: str):
    seqs = []
    with open(path, "r") as f:
        lines = f.readlines()
        till_seq = 2
        for line in lines:
            line = line.rstrip()
            if line.__contains__('>'):
                till_seq = 2
                pass
            elif till_seq == 1:
                seqs.append(line)
                till_seq = 0
                pass
            else:
                till_seq = 1
                pass
            pass
        pass
    seqs.append("ATGGGACGCAAGGACGAGCAGAAGCAAACGAGCGCCACAAGCACGCCGGGGCAGGGG")
    return seqs


def int_to_DNA(vals: []):
    rtn = ''
    for val in vals:
        if val == 0:
            rtn += 'G'
        elif val == 1:
            rtn += 'C'
        elif val == 2:
            rtn += 'A'
        elif val == 3:
            rtn += 'T'
            pass
        pass
    return rtn


def DNA_to_int(seq: str):
    rtn = []
    for c in seq:
        if c == 'g' or c == 'G':
            rtn.append(0)
        elif c == 'c' or c == 'C':
            rtn.append(1)
        elif c == 'a' or c == 'A':
            rtn.append(2)
        elif c == 't' or c == 'T':
            rtn.append(3)
            pass
        pass
    return rtn


def main():
    print("START")
    folder_names = os.listdir(inp)
    seq_idxs = [0, 1, 2, 3, 4, 5]
    seq_idxs = [1, 2, 3, 4]
    popsizes = [500, 5000]
    num_states = [5, 20]
    max_muts = [1, 3, 10]
    tourn_sizes = [5, 15]
    sequences = gen_sequences("./Sequences.dat")
    for idx in range(len(sequences)):
        sequences[idx] = DNA_to_int(sequences[idx])
        print(len(sequences[idx]))
        print(int_to_DNA(sequences[idx]))
        pass
    return 0


    exp_dat = []
    exp_lbls = []
    exp_descriptions = []
    exp_num = 1
    for ps in popsizes:
        for ns in num_states:
            for mm in max_muts:
                for ts in tourn_sizes:
                    exp_dat.append([str(ps).rjust(4, '0') + "Pop", str(ns).rjust(2, '0') + "Sta",
                                    str(mm).rjust(2, '0') + "Mut", str(ts).rjust(2, '0') + "Tsz"])
                    exp_lbls.append(str(exp_num) + "(" + str(ps).rjust(4, '_') + ", " + str(ns).rjust(2, '_') +
                                    ", " + str(mm).rjust(2, '_') +
                                    ", " + str(ts).rjust(2, '_') + ")")
                    exp_descriptions.append(str(ps).rjust(4, '0') + "Pop, " + str(ns).rjust(2, '0') + "Sta, " +
                                            str(mm).rjust(2, '0') + "Mut, " +
                                            str(ts).rjust(2, '0') + "Tsz")
                    exp_num += 1
                    pass
                pass
            pass
        pass

    seq_dirs = []
    for seq in seq_idxs:
        one_dirs = []
        for eidx, dat in enumerate(exp_dat):
            for fld in folder_names:
                if all(fld.__contains__(itm) for itm in dat) and str("Seq" + str(seq)) in fld:
                    one_dirs.append(fld)
                    pass
                pass
            pass
        seq_dirs.append(one_dirs)
        pass

    # mode_data[seq][exp][run][0] = run num
    # mode_data[seq][exp][run][1] = run's fit (sorted based on this)
    # mode_data[seq][exp][run][2][:] = run's SDA
    # mode_data[seq][exp][run][3] = run's sequence
    mode_data = [[] for _ in range(len(seq_idxs))]
    for sidx, seq in enumerate(seq_idxs):
        for fld in seq_dirs[sidx]:
            mode_data[sidx].append(get_data(inp + fld + "/"))
            pass
        pass

    # mode_stats[seq][exp] = [run's fitness vals]
    mode_stats = [[] for _ in range(len(seq_idxs))]
    make_all = False
    make_any = False
    for sidx, seq in enumerate(seq_idxs):
        for expidx, exp in enumerate(mode_data[sidx]):
            exp_fits = []
            for runidx, run in enumerate(exp):
                exp_fits.append(run[1])
                pass
            mode_stats[sidx].append(exp_fits)
            pass
        pass

    title = "Sequence Matching with Sequence "
    # xsp = [[i for i in range(len(all_data[0]))], [i for i in range(len(all_data[1]))]]
    # xpos = [xsp[0], xsp[1], xsp[0], xsp[1], xsp[0], xsp[1], xsp[0], xsp[1]]
    ylb = "Fitness"
    xlb = "Experiment (Popsize, SDA States, Max Muts., Tourn. Size)"

    # lxpos = []
    # for i in range(2, len(all_data[0]) - 3, 3):
    #     lxpos.append(i + 0.5)
    #     pass
    # base_colors = ["#808080", "#696969"]
    colors = ['#0077AA','#0000FF', '#00AA00', '#00FF00','#AAAA00', '#FFFF00']
    for sidx, seq in enumerate(seq_idxs):
        plt.style.use("seaborn-v0_8")
        plt.rc('xtick', labelsize=10)
        plt.rc('ytick', labelsize=10)

        f = plt.figure()
        f.set_figheight(4.5)
        f.set_figwidth(8)
        plot = f.add_subplot(111)

        bp = plot.boxplot(mode_stats[sidx], patch_artist=True)
        box_plot(bp, 1, [[[i for i in range(len(mode_stats[sidx]))], colors]])

        # plot.set_xticks(xpos[idx])
        plot.set_xticklabels(exp_lbls, rotation=90)

        # f.suptitle(title + str(seq), fontsize=14)
        plot.set_xlabel(xlb, fontsize=12)
        plot.set_ylabel(ylb, fontsize=12)

        plot.hlines([len(sequences[seq])], 0.5, len(mode_stats[sidx]) + 0.5, colors="#0000FF", linestyles="dashed", linewidth=1)
        # for x in lxpos:
        #     plot.axvline(x=x, color='black', linestyle='--', linewidth=0.75)
        #     pass
        plot.grid(visible="True", axis="y", which='major', color="darkgray", linewidth=0.75)
        f.tight_layout()
        f.savefig(outp + "AAMatchSeq" + str(seq) + "_boxplot.png", dpi=300)
        plt.close()
        pass

    for sidx, seq in enumerate(seq_idxs):
        best_runs = []
        best_of_best_exp = -1
        best_of_best_val = 0
        for didx, dat in enumerate(mode_data[sidx]):
            best_runs.append(dat[0][0])
            if dat[0][1] > best_of_best_val:
                best_of_best_val = dat[0][1]
                best_of_best_exp = didx
            pass
        make_table(mode_stats[sidx], best_runs, exp_descriptions, outp +
                   "Seq" + str(seq) + "table" + ".dat", True)
        info = "Best for Sequence " + str(seq) + " is " + \
               "EXP" + str(best_of_best_exp + 1) + ": " + str(exp_descriptions[best_of_best_exp])
        print_best_info(outp + "Seq" + str(seq) + "_best.dat", info,
                        mode_data[sidx][best_of_best_exp][0], sequences[seq])
        pass

    # mode_data[mode][exp][run][0] = run num
    # mode_data[mode][exp][run][1] = run's fit (sorted based on this)
    # mode_data[mode][exp][run][2][:] = run's SDA
    # mode_data[mode][exp][run][3][:] = run's network
    # mode_data[mode][exp][run][4] = run's network's edges
    # mode_data[mode][exp][run][5] = run's network's weight
    # mode_data[mode][exp][run][6][:] = run's network's weight hist
    # net_stats = [[[[] for _ in range(2)] for _ in range(len(sizes))] for _ in range(len(mode_itms))]
    # net_bests = [[[[] for _ in range(2)] for _ in range(len(sizes))] for _ in range(len(mode_itms))]
    # for midx, itms in enumerate(mode_itms):
    #     for sidx, size in enumerate(sizes):
    #         for expidx, exp in enumerate(mode_data[midx][sidx]):
    #             exp_edges = []
    #             exp_weights = []
    #             for runidx, run in enumerate(exp):
    #                 if runidx == 0:
    #                     net_bests[midx][sidx][0].append(run[4])
    #                     net_bests[midx][sidx][1].append(run[5])
    #                     pass
    #                 exp_edges.append(run[4])
    #                 exp_weights.append(run[5])
    #                 pass
    #             net_stats[midx][sidx][0].append(exp_edges)
    #             net_stats[midx][sidx][1].append(exp_weights)
    #             pass
    #         pass
    #     pass

    # titles = ["Epidemic Length", "Epidemic Profile Matching P1", "Epidemic Profile Matching P7",
    #           "Epidemic Profile Matching Dublin"]
    # names = ["EL_netstats", "PM1_netstats", "PM7_netstats", "PMDUB_netstats"]
    # ylb = ["Network Edges", "Network Weight"]
    # xlb = ["Experiment (Num. States, Max. Mutations)",
    #        "Experiment (Num. States, Max. Mutations)",
    #        "Experiment (Num. States, Max. Mutations)",
    #        "Experiment (Num. States, Max. Mutations)", ]
    # out_path = outp + "Network Stats Boxplots"
    # if not os.path.exists(out_path):
    #     os.makedirs(out_path)
    #     pass
    # out_path += "/"
    # xsp = [[i for i in range(1, 10)], [i - 0.22 for i in range(1, 10)], [i + 0.22 for i in range(1, 10)]]
    # for idx in range(len(titles)):
    #     for sidx, size in enumerate(sizes):
    #         if idx < 3 or sidx == 0:
    #             if idx >= 3:
    #                 size = 200
    #                 pass
    #             plt.rc('xtick', labelsize=6)
    #             plt.rc('ytick', labelsize=6)
    #
    #             f = plt.figure()
    #             f.set_figheight(5)
    #             f.set_figwidth(8)
    #             plot = f.add_subplot(111)
    #             plot2 = plot.twinx()
    #
    #             # plot2.set_axisbelow(True)
    #             # plot2.grid(visible="True", axis="y", which='major', color="lightgray", linewidth=0.75)
    #
    #             bp = plot.boxplot(net_stats[idx][sidx][0], positions=xsp[1], patch_artist=True, zorder=1, widths=0.4)
    #             box_plot(bp, 1, [[[i for i in range(9)], ["#0000FF"]]])
    #             plot.plot(xsp[1], net_bests[idx][sidx][0], 'x', color="#FF0000")
    #             bp = plot2.boxplot(net_stats[idx][sidx][1], positions=xsp[2], patch_artist=True, zorder=2, widths=0.4)
    #             box_plot(bp, 1, [[[i for i in range(9)], ["#00FF00"]]])
    #             plot2.plot(xsp[2], net_bests[idx][sidx][1], 'x', color="#FF0000")
    #
    #             plot.hlines([size * 1, size * 4], 0.5, 9.5, colors="#0000FF", linestyles="dashed", linewidth=1)
    #             plot2.hlines([size * 4, size * 16], 0.5, 9.5, colors="#00FF00", linestyles="dotted", linewidth=1)
    #
    #             plot.set_xticks(xsp[0])
    #             plot.set_xticklabels(exp_lbls[2:], rotation=90)
    #
    #             if not titles[idx].__contains__("Dublin"):
    #                 f.suptitle(titles[idx] + " w " + str(size) + " Nodes", fontsize=12)
    #             else:
    #                 f.suptitle(titles[idx], fontsize=12)
    #                 pass
    #
    #             plot.set_xlabel(xlb[idx], fontsize=10)
    #             plot.set_ylabel(ylb[0], fontsize=10, color="#0000FF")
    #             plot2.set_ylabel(ylb[1], fontsize=10, rotation=270, labelpad=10, color="#00FF00")
    #
    #             plot.set_axisbelow(True)
    #             plot.grid(visible="True", axis="y", which='major', color="darkgray", linewidth=0.75)
    #             f.tight_layout()
    #
    #             if idx < 3:
    #                 f.savefig(out_path + str(size) + names[idx] + ".png", dpi=450)
    #             else:
    #                 f.savefig(out_path + names[idx] + ".png", dpi=450)
    #                 pass
    #             plt.close()
    #         pass
    #     pass
    print("END")
    pass


main()
