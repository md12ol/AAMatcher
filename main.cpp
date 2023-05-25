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
 * 7.   Maximum number of mutations
 * 8.   Sequence number
 * 9.   Tournament size
 * 10.  Crossover operator
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
    double expBestFit = (BIGGERBETTER ? 0 : MAXFLOAT);
    SDA expBestSDA = SDA(sdaStates, numChars, 2, seqLen);

    initAlg(pathToSeqs);
    cmdLineIntro(cout);
    sprintf(pathToOut, "./AAMOut/AAMatch on Seq%d with %04dPop, %02dSta, %02dMut, %02dTsz, %dCross/", seqNum, popsize,
            sdaStates, maxMuts, tournSize, crossOp);
    mkdir(pathToOut, 0777);
    expStats.open(string(pathToOut) + "./exp.dat", ios::out);
    readMe.open(string(pathToOut) + "./read.me", ios::out);
    makeReadMe(readMe);
    readMe.close();

    int tmp;
    for (int run = 1; run < runs + 1; ++run) {
        initPop(run);
        filename = string(pathToOut) + "run" + to_string(run) + ".dat";
        runStats.open(filename, ios::out);
        printExpStatsHeader(cout);
        printExpStatsHeader(runStats);
        report(runStats, run, 0, BIGGERBETTER);

        int gen = 1;
        int stallCount = 0;
        double best = (BIGGERBETTER ? 0 : MAXFLOAT);
        while (gen <= maxGens && stallCount < TERM_CRIT) {
            matingEvent(BIGGERBETTER);

            if (gen % REPORT_EVERY == 0) {
                tmp = report(runStats, run, (int) gen / (REPORT_EVERY), BIGGERBETTER);
                if ((BIGGERBETTER && tmp > best) || (!BIGGERBETTER && tmp < best)) {
                    best = tmp;
                    stallCount = 0;
                } else {
                    stallCount++;
                }
            }

            if (gen % (int) (CULLING_EVERY * REPORT_EVERY) == 0) {
                culling(CULLING_ODDS, RANDOM_CULLING, BIGGERBETTER);
            }
            gen++;
        }

        tmp = runReport(expStats, BIGGERBETTER);
        if ((BIGGERBETTER && fits[tmp] > expBestFit) || (!BIGGERBETTER && fits[tmp] < expBestFit)) {
            expBestFit = fits[tmp];
            expBestSDA.copy(pop[tmp]);
        }
        bests.push_back(fits[tmp]);
        runStats.close();
    }

    ofstream best;
    best.open(string(pathToOut) + "./best.dat", ios::out);
    expReport(best, bests, expBestSDA, BIGGERBETTER);
    best.close();
    delete[] pop;
    cout << "Program Completed Successfully!" << endl;
    return 0;
}





