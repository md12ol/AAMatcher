#include <iostream>
#include "main.h"

/**
 * This method...
 *
 * Input variables:
 * 1.   Population size
 * 2.   Number of characters (in SDA)
 * 3.   Number of states
 * 4.   Seed
 * 5.   Number of runs
 * 6.   Maximum number of mating events
 * 7.   Default Number of Transition Mutations
 * 8.   Default Number of Response Mutations
 * 9.   Dynamic Mutation Operator? (0 -> static, >0 -> dynamic implementation to use)
 * 10.  Upper Bound on Mutations
 * 11.  Sequence number
 * 12.  Tournament size
 * 13.  Crossover operator
 * 14.  Crossover Rate
 * 15.  Mutation Rate
 * 16.  Culling Rate
 * 17.  Random Culling
 * 18.  Culling Every _ Reporting Intervals
 *
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char *argv[]) {
    getArgs(argv);
    string pathToSeqs = "./Sequences.dat";
    string filename;
    ofstream runStats, expStats, readMe;
    vector<double> bests;
    bests.reserve(runs);
    double expBestFit = (BIGGER_BETTER ? 0 : MAXFLOAT);
    SDA expBestSDA = SDA(sdaStates, numChars, 2, seqLen);

    initAlg(pathToSeqs);
    cmdLineIntro(cout);
    char dynamicMessage[20];
    sprintf(dynamicMessage, "%s%d", (dynamicMutOperator == 0 ? "Static" : "Dynamic"), dynamicMutOperator);
    sprintf(pathToOut, "./AAMOut/AAMatch on Seq%d with %.1fmilMMEs, %04dPS, %02dSt, %02dNTM, %02dNRM, %s,"
                       " %dTS, %sCO, %03d%%CrR, %03d%%MR, %03d%%CuR, %sCu/", seqNum, (double) maxGens / 1000000,
            popsize, sdaStates, initNumTransMuts, initNumRespMuts, dynamicMessage, tournSize,
            (crossoverOp == 0 ? "2Pt" : "1St"), (int) (crossoverRate * 100), (int) (mutationRate * 100),
            (int) (cullingRate * 100), (randomCulling ? "Rand" : "Worst"));
    mkdir(pathToOut, 0777);
    expStats.open(string(pathToOut) + "./exp.dat", ios::out);
    readMe.open(string(pathToOut) + "./read.me", ios::out);
    makeReadMe(readMe);
    readMe.close();

    int tmp;
    for (int run = 1; run < runs + 1; ++run) {
        curNumTransMuts = initNumTransMuts;
        curNumRespMuts = initNumRespMuts;
        initPop(run);
        char runNumStr[10];
        sprintf(runNumStr, "%02d", run);
        filename = string(pathToOut) + "run" + string(runNumStr) + ".dat";
        runStats.open(filename, ios::out);
        printExpStatsHeader(cout);
        printExpStatsHeader(runStats);
        report(runStats, run, 0, BIGGER_BETTER);

        int gen = 1;
        int stallCount = 0;
        double best = (BIGGER_BETTER ? 0 : MAXFLOAT);
        while (gen <= maxGens && stallCount < TERM_CRIT) {
            matingEvent(BIGGER_BETTER);

            if (gen % REPORT_EVERY == 0) {
                tmp = report(runStats, run, (int) gen / (REPORT_EVERY), BIGGER_BETTER);
                if ((BIGGER_BETTER && tmp > best) || (!BIGGER_BETTER && tmp < best)) {
                    best = tmp;
                    stallCount = 0;
                } else {
                    stallCount++;
                }
            }

            if (gen % (int) (CULLING_EVERY * REPORT_EVERY) == 0 && stallCount < TERM_CRIT) {
                culling(cullingRate, randomCulling, BIGGER_BETTER);
            }
            gen++;
        }

        tmp = runReport(expStats, BIGGER_BETTER);
        if ((BIGGER_BETTER && fits[tmp] > expBestFit) || (!BIGGER_BETTER && fits[tmp] < expBestFit)) {
            expBestFit = fits[tmp];
            expBestSDA.copy(pop[tmp]);
        }
        bests.push_back(fits[tmp]);
        runStats.close();
    }

    ofstream best;
    best.open(string(pathToOut) + "./best.dat", ios::out);
    expReport(best, bests, expBestSDA, BIGGER_BETTER);
    best.close();
    delete[] pop;
    cout << "Program Completed Successfully!" << endl;
    return 0;
}





