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
 * 17.  Random Culling (1 -> Random, 0 -> Worst %)
 * 18.  Culling Every _ Reporting Intervals
 * 19.  First Run
 *
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char *argv[]) {
    getArgs(argc, argv);
    string pathToSeqs = "./Sequences.dat";
    string filename;
    ofstream runStats, expStats, readMe, sdaFile;

    vector<pair<int, int>> bests;
    bests.reserve(runs);
    pair<int, int> expBestFit = make_pair((BIGGER_BETTER ? 0.0 : MAXFLOAT), (false ? INT16_MAX : 0.0));
    SDA expBestSDA = SDA(sdaStates, numChars, 2, seqLen);

    initAlg(pathToSeqs);
    cmdLineIntro(cout);
    char dynamicMessage[20];
    char startRun[20];
    // TODO: Add dynamic mutation update.
    sprintf(dynamicMessage, "%s%d", (dynamicMutOperator == 0 ? "Static" : "Dynamic"), dynamicMutOperator);
    sprintf(startRun, "%02d", initRun);
    sprintf(pathToOut, "./AAMOut/AAMatch on Seq%d with %.1fmilMMEs, %04dPS, %02dSt, %02dNTM, %02dNRM, %s,"
                       " %dTS, %sCO, %03d%%CrR, %03d%%MR, %03d%%CuR, %sCu, %dCE, V2/", seqNum, (double) maxGens / 1000000,
            popsize, sdaStates, initNumTransMuts, initNumRespMuts, dynamicMessage, tournSize,
            (crossoverOp == 0 ? "2Pt" : "1St"), (int) (crossoverRate * 100), (int) (mutationRate * 100),
            (int) (cullingRate * 100), (randomCulling ? "Rand" : "Worst"), CULLING_EVERY);
    mkdir(pathToOut, 0777);
    expStats.open(string(pathToOut) + "./exp" + string(startRun) + ".dat", ios::out);
    readMe.open(string(pathToOut) + "./read.me", ios::out);
    makeReadMe(readMe);
    readMe.close();

    pair<double, double> rptVals;
    for (int run = initRun; run < initRun + runs + 1; ++run) {
        curNumTransMuts = initNumTransMuts;
        curNumRespMuts = initNumRespMuts;
        initPop(run);
        char runNumStr[20];
        sprintf(runNumStr, "%02d", run);
        filename = string(pathToOut) + "run" + string(runNumStr) + ".dat";
        runStats.open(filename, ios::out);
        filename = string(pathToOut) + "SDAs" + string(runNumStr) + ".dat";
        sdaFile.open(filename, ios::out);
        printExpStatsHeader(cout);
        printExpStatsHeader(runStats);
        report(runStats, run, 0, BIGGER_BETTER);
        sdaCheck(sdaFile, 0);

        int gen = 1;
        int stallCount = 0;
        pair<int, int> runBestFit = make_pair((BIGGER_BETTER ? 0.0 : MAXFLOAT), (false ? INT16_MAX : 0.0));
        while (gen <= maxGens && (stallCount < TERM_CRIT || gen <= MIN_GEN_RATIO * maxGens)) {
            matingEvent(BIGGER_BETTER);
            if (gen % REPORT_EVERY == 0) {
                rptVals = report(runStats, run, (int) gen / (REPORT_EVERY), BIGGER_BETTER);
//                printPopFits(cout);
//                cout << "[" << populationBestFit.first << ", " << populationBestFit.second << "]" << endl;

//                vector<int> sortedIdxs = tournSelect(popsize, BIGGER_BETTER);
//                cout << "Fitness Values: ";
//
//                bool first = true;
//                for (int idx: sortedIdxs) {
//                    if (!first) {
//                        cout << ", ";
//                    }
//                    cout << "[" << matchFits[idx] << ", " << noveltyFits[idx] << "]";
//                    first = false;
//                }
//                cout << "[" << populationBestFit.first << ", " << populationBestFit.second << "]" << endl;

//                cout << "Match Fits: ";
//                vector<double> thing = matchFits;
//                sort(thing.begin(), thing.end());
//                reverse(thing.begin(), thing.end());
//                for (int i = 0; i < popsize; i++) {
//                    cout << setw(2) << thing[i] << ", ";
//                }
//                cout << endl;
//
//                cout << "Nove. Fits: ";
//                for (int i = 0; i < popsize; i++) {
//                    cout << setw(2) << noveltyFits[i] << ", ";
//                }
//                cout << endl;

                if ((BIGGER_BETTER && rptVals.first > runBestFit.first) ||
                    (!BIGGER_BETTER && rptVals.first < runBestFit.first)) {
                    runBestFit = rptVals;
                    stallCount = 0;
                } else {
                    runBestFit.second = rptVals.second;
                    stallCount++;
                }

                if (gen == 2500000) {
                    sdaCheck(sdaFile, gen);
                }

                // TODO: Add dynamic mutation update.
            }

            if (gen < maxGens && gen % (int) (CULLING_EVERY * REPORT_EVERY) == 0 &&
                    (stallCount < TERM_CRIT || gen < MIN_GEN_RATIO * maxGens)) {
                culling(cullingRate, randomCulling, BIGGER_BETTER);
            }
            gen++;
        }
        runReport(expStats, BIGGER_BETTER, make_pair(matchFits[popBestIdx], noveltyFits[popBestIdx]));
        if ((BIGGER_BETTER && runBestFit.first > expBestFit.first) ||
            (!BIGGER_BETTER && runBestFit.first < expBestFit.first)) {
            expBestFit = runBestFit;
            expBestSDA.copy(pop[popBestIdx]);
        } else if (runBestFit.first == expBestFit.first) {
            if (runBestFit.second < expBestFit.second) {
                expBestFit = runBestFit;
                expBestSDA.copy(pop[popBestIdx]);
            }
        }
        bests.push_back(rptVals);
        sdaCheck(sdaFile, gen);
        runStats.close();
        sdaFile.close();
    }

    ofstream best;
    best.open(string(pathToOut) + "./best" + string(startRun) + ".dat", ios::out);
    expReport(best, bests, expBestSDA, BIGGER_BETTER);
    best.close();
    delete[] pop;
    cout << "Program Completed Successfully!" << endl;
    return 0;
}
