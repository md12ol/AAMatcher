from random import randint

exec_path = "./cmake-build-release---remote/AAMatcher"
popsizes = [500, 5000]
num_chars = 4
num_states = [5, 20, 50]
num_runs = 50
gens = 10000000
max_muts = [1, 3, 10]
seq_nums = [5]
tourn_sizes = [5, 15]


def main():
    with open("./table.dat", "w") as f:
        for sn in seq_nums:
            for ps in popsizes:
                for ns in num_states:
                    for mm in max_muts:
                        for ts in tourn_sizes:
                            rnd = randint(1000, 9999)
                            line = exec_path + " " + str(ps) + " " + str(num_chars) + " " + str(ns) + " " + \
                                   str(rnd) + " " + str(num_runs) + " " + str(gens) + " " + str(mm) + " " + \
                                   str(sn) + " " + str(ts)
                            f.write(line + "\n")
                            pass
                        pass
                    pass
                pass
            pass
        pass
    pass


main()
