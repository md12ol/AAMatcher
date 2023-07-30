from random import randint

exec_path = "./cmake-build-release---graham/SDATester"
popsizes = [50, 500]
num_chars = 4
num_states = [5]
num_runs = 50
gens = 10000000
num_muts = [2, 4]
seq_nums = [0, 1, 2, 3, 4, 5]
tourn_sizes = [7]
crossOp = [0]
crossRate = [0.5]
mutateRate = [1]
cullingRate = [0.25]
randomCulling = [1]
dynamic_muts = [0, 1, 2, 3, 4]
culling_every = [1, 5]


def main():
    with open("./table2Graham.dat", "w") as f:
        for sn in seq_nums:
            for ps in popsizes:
                for ns in num_states:
                    for mm in num_muts:
                        for dm in num_muts:
                            for dm2 in dynamic_muts:
                                for ts in tourn_sizes:
                                    for cop in crossOp:
                                        for cra in crossRate:
                                            for mr in mutateRate:
                                                for cur in cullingRate:
                                                    for rc in randomCulling:
                                                        for ce in culling_every:
                                                            if cur == 1 and randomCulling == 1:
                                                                pass
                                                            else:
                                                                rnd = randint(1000, 9999)
                                                                line = exec_path + " " + str(ps).zfill(3) + " " + str(
                                                                    num_chars) + " " + str(ns) + " " + str(rnd) + " " + str(
                                                                    num_runs) + " " + str(gens) + " " + str(mm) + " " + str(
                                                                    dm) + " " + str(dm2) + " " + " " + str(sn) + " " + str(
                                                                    ts) + " " + str(cop) + " " + "%0.2f" % cra + " " + \
                                                                       "%0.2f" % mr + " " + "%0.2f" % cur + " " + str(
                                                                    rc) + " " + str(ce)
                                                                pass
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
