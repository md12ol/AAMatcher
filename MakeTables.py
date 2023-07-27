from random import randint

exec_path = "./cmake-build-release---remote/SDATester"
popsizes = [50]
num_chars = 4
num_states = [6, 12]
num_runs = 50
gens = 1000000
num_muts = [2, 4]
seq_nums = [1]
tourn_sizes = [5]
crossOp = [0]
crossRate = [1]
mutateRate = [1, 0.5]
cullingRate = [1.0]
randomCulling = [0]
dynamic_muts = [0]


def main():
    with open("./table.dat", "w") as f:
        for sn in seq_nums:
            for ps in popsizes:
                for ns in num_states:
                    for mm in num_muts:
                        for dm in dynamic_muts:
                            for ts in tourn_sizes:
                                for cop in crossOp:
                                    for cra in crossRate:
                                        for mr in mutateRate:
                                            for cur in cullingRate:
                                                for rc in randomCulling:
                                                    if cur == 1 and randomCulling == 1:
                                                        pass
                                                    else:
                                                        rnd = randint(1000, 9999)
                                                        line = exec_path + " " + str(ps) + " " + str(
                                                            num_chars) + " " + str(ns) + " " + str(rnd) + " " + str(
                                                            num_runs) + " " + str(gens) + " " + str(mm) + " " + str(
                                                            mm) + " " + str(dm) + " " + str(sn) + " " + str(
                                                            ts) + " " + str(cop) + " " + "%0.2f" % cra + " " + \
                                                               "%0.2f" % mr + " " + "%0.2f" % cur + " " + str(rc)
                                                f.write(line + "\n")
                                                pass
                                            pass
                                        pass
                                    pass
                                pass
                            pass
                        pass
                    pass
                pass
            pass
        pass
    pass
    # with open("./table.dat", "w") as f:
    #     for sn in seq_nums:
    #         for ps in popsizes:
    #             for ns in num_states:
    #                 for mm in max_muts:
    #                     for ts in tourn_sizes:
    #                         for cop in crossOp:
    #                             for cra in crossRate:
    #                                 for mr in mutateRate:
    #                                     for cur in cullingRate:
    #                                         for rc in randomCulling:
    #                                             if cur == 1 and randomCulling == 1:
    #                                                 pass
    #                                             else:
    #                                                 rnd = randint(1000, 9999)
    #                                                 line = exec_path + " " + str(ps) + " " + str(num_chars) + " " + \
    #                                                        str(ns) + " " + str(rnd) + " " + str(num_runs) + " " + \
    #                                                        str(gens) + " " + str(mm) + " " + str(sn) + " " + \
    #                                                        str(ts) + " " + str(cop) + " " + "%0.2f" % cra + " " + \
    #                                                        "%0.2f" % mr + " " + "%0.2f" % cur + " " + \
    #                                                        str(rc)
    #                                             f.write(line + "\n")
    #                                             pass
    #                                         pass
    #                                     pass
    #                                 pass
    #                             pass
    #                         pass
    #                     pass
    #                 pass
    #             pass
    #         pass
    #     pass


print("DONE!")
pass

main()
