import os
import random

prog_name = "./cmake-build-release---graham/AAMatcher"


def main():
    # popsize = ["50", "250", "500"]
    popsize = ["500"]
    first = "4 20"
    # second = "50 10000000 2 2 0 10"
    second = "5 10000000 2 2 0 10"
    sequence = ["0", "1", "2", "3", "4", "5"]
    tourn_size = ["7", "15"]
    # third = "0 0.5 1.0 0.25 1 5 1"
    third = "0 0.5 1.0 0.25 1 5"
    start_run = ["1", "6", "11", "16", "21", "26", "31", "36", "41", "46"]

    with open("./table2.dat", "w") as f:
        for ps in popsize:
            for seq in sequence:
                for ts in tourn_size:
                    for sr in start_run:
                        f.writelines(prog_name + " " + ps + " " + first + " " + str(random.randint(1000, 9999)) + " "
                                    + second + " " + seq + " " + ts + " " + third + " " + sr + "\n")
                    pass
                pass
            pass
        pass
    print("DONE")
    pass


main()
