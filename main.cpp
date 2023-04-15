#include <iostream>
#include "main.h"

int main(int argc, char *argv[]) {
    double bestOfBest = (BIGGERBETTER ? 0 : MAXFLOAT);
    string filename;
    ofstream runStats, expStats, readMe;

    vector<int> bests;
    bests.reserve(runs);

    getArgs(argv);
    initAlg();
    cmdLineIntro(cout);
    sprintf(pathToOut, "./AAMOut/AAMatch on Seq%d with %04dPop, %02dSta, %02dMut, %02dTsz/",
            seqNum, popsize, sdaStates, maxMuts, tournSize);
    filesystem::create_directory(pathToOut);
    expStats.open(string(pathToOut) + "./exp.dat", ios::out);
    readMe.open(string(pathToOut) + "./read.me", ios::out);
    makeReadMe(readMe);
    readMe.close();

    double tmp;
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

        tmp = finalReport(expStats, true);
        if ((BIGGERBETTER && tmp > bestOfBest) || (!BIGGERBETTER && tmp < bestOfBest)) {
            bestOfBest = tmp;
        }
        bests.push_back(tmp);
        runStats.close();
    }

    vector<double> stats = calcStats2(bests, BIGGERBETTER);
    cout << "BEST OF BEST: " << bestOfBest << " of " << seqLen << endl;
    cout << "Fitness for Exp: " << stats[0] << " +- " << stats[2] << endl;
    delete[] pop;
    cout << "Program Completed Successfully!" << endl;
    return 0;
}

int culling(double percentage, bool rndPick, bool biggerBetter) {
    int numKillings = (int) (popsize * percentage);
    vector<int> contestants;
    vector<int> winners;
    winners.reserve(numKillings);

    if (rndPick) {
        contestants = tournSelect(numKillings, !biggerBetter);
    } else {
        contestants = tournSelect(popsize, !biggerBetter);
    }

    for (int idx = 0; idx < numKillings; idx++) {
        winners.push_back(contestants[idx]);
    }

    for (int idx: winners) {
        pop[idx].randomize();
        fits[idx] = fitness(pop[idx]);
    }

    return 0;
}

double fitness(Bitsprayer &sda) {
    double val = 0.0;
    vector<int> result = sda.getBitsVec(seqLen);
    for (int i = 0; i < seqLen; ++i) {
        if (result[i] == testSeq[i]) {
            val += 1;
        }
    }
    return val;
}

bool compareFitness(int popIdx1, int popIdx2) {
    return fits[popIdx1] < fits[popIdx2];
}

vector<int> tournSelect(int size, bool decreasing) {
    vector<int> tournIdxs;
    int idxToAdd;

    tournIdxs.reserve(size);
    if (size == popsize) {
        for (int idx = 0; idx < size; idx++) {
            tournIdxs.push_back(idx);
        }
    } else {
        do {
            idxToAdd = (int) lrand48() % popsize;
            if (count(tournIdxs.begin(), tournIdxs.end(), idxToAdd) == 0) {
                tournIdxs.push_back(idxToAdd);
            }
        } while (tournIdxs.size() < tournSize);
    }

    sort(tournIdxs.begin(), tournIdxs.end(), compareFitness);
    if (decreasing) {
        reverse(tournIdxs.begin(), tournIdxs.end());
    }
    return tournIdxs;
}

int matingEvent(bool biggerBetter) {
    int numMuts;
    Bitsprayer child1, child2;

    vector<int> idxs = tournSelect(tournSize, biggerBetter);
//    printIdxsOfVector<double>(fits, idxs, "Tournament Fitness Values: ", "\t", true);

    child1.copy(pop[idxs[0]]);
    child2.copy(pop[idxs[1]]);
    child1.twoPtCrossover(child2);

    numMuts = (int) lrand48() % maxMuts + 1;
    child1.mutate(numMuts);
    numMuts = (int) lrand48() % maxMuts + 1;
    child2.mutate(numMuts);

    pop[idxs.end()[-1]] = child1;
    pop[idxs.end()[-2]] = child2;

    fits[idxs.end()[-1]] = fitness(child1);
    fits[idxs.end()[-2]] = fitness(child2);
    return 0;
}

vector<double> calcStats(bool biggerBetter) {
    double sum = 0.0;
    double bestVal = (biggerBetter ? 0.0 : MAXFLOAT);

    for (double fit: fits) {
        sum += fit;
        if ((biggerBetter && fit > bestVal) || (!biggerBetter && fit < bestVal)) {
            bestVal = fit;
        }
    }
    double mean = sum / popsize;
    double stdDevSum = 0.0;
    for (double fit: fits) {
        stdDevSum += pow(fit - mean, 2);
    }
    double stdDev = sqrt(stdDevSum / (popsize - 1));
    double CI95 = 1.96 * (stdDev / sqrt(popsize));

    return {mean, stdDev, CI95, (double) bestVal}; // {mean, stdDev, 95CI, best}
}

vector<double> calcStats2(vector<int> vals, bool biggerBetter) {
    double sum = 0.0;
    double bestVal = (biggerBetter ? 0.0 : MAXFLOAT);

    for (int val: vals) {
        sum += val;
        if ((biggerBetter && val > bestVal) || (!biggerBetter && val < bestVal)) {
            bestVal = val;
        }
    }
    double mean = sum / (double) vals.size();
    double stdDevSum = 0.0;
    for (int val: vals) {
        stdDevSum += pow((double) val - mean, 2);
    }
    double stdDev = sqrt(stdDevSum / ((double) vals.size() - 1.0));
    double CI95 = 1.96 * (stdDev / sqrt(vals.size()));

    return {mean, stdDev, CI95, bestVal}; // {mean, stdDev, 95CI, best}
}

double report(ofstream &outp, int run, int rptNum, bool biggerBetter) {
    vector<double> stats = calcStats(biggerBetter); // {mean, stdDev, 95CI, best}
    multiStream printAndSave(cout, outp);

    printAndSave << left << setw(5) << run;
    printAndSave << left << setw(4) << rptNum;
    printAndSave << left << setw(10) << stats[0];
    printAndSave << left << setw(12) << stats[2];
    printAndSave << left << setw(10) << stats[1];
    printAndSave << left << setw(8) << stats[3];
    printAndSave << left << setw(8) << (stats[3] / seqLen) * 100 << "%";
    printAndSave << "\n";
    return stats[3];
}

double finalReport(ostream &outp, bool biggerBetter) {
    auto maxIterator = minmax_element(fits.begin(), fits.end());
    int bestIdx;
    if (biggerBetter) {
        bestIdx = (int) distance(fits.begin(), maxIterator.second);
    } else {
        bestIdx = (int) distance(fits.begin(), maxIterator.first);
    }

    multiStream printAndSave(cout, outp);
    printAndSave << "The best fitness is " << fits[bestIdx] << "\n";
    vector<int> bestSeq = pop[bestIdx].getBitsVec(seqLen);
    printAndSave << left << setw(20) << "Best Match: ";
    printVector<multiStream, int>(printAndSave, bestSeq, "", "", false);
    printAndSave << "\n";
    printAndSave << left << setw(20) << "Matches: ";
    for (int idx = 0; idx < seqLen; idx++) {
        if (testSeq[idx] == bestSeq[idx]) {
            printAndSave << "X";
        } else {
            printAndSave << " ";
        }
    }
    printAndSave << "\n";
    printAndSave << left << setw(20) << "Desired Sequence: ";
    printVector<multiStream, int>(printAndSave, testSeq, "", "", false);
    printAndSave << "\n";

    outp << "SDA" << endl;
    pop[bestIdx].print(outp);

    sort(fits.begin(), fits.end());
    if (biggerBetter) {
        reverse(fits.begin(), fits.end());
    }
    printVector<multiStream, double>(printAndSave, fits, "Fitness Values: ", " ", true);
    cout << endl;
    return fits[0];
}





