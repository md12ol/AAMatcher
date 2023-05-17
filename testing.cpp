#include "testing.h"

vector<int> runMultiMating(int numEvents, SDA mom, SDA dad) {
    SDA child1, child2;
    vector<int> fitVals;

    fitVals.reserve(popsize);

    for (int i = 0; i < numEvents; ++i) {
        child1.copy(mom);
        child2.copy(dad);
        if (crossOp == 0) child1.twoPtCrossover(child2);
        else if (crossOp == 1) child1.oneStateCrossover(child2);
        fitVals.push_back((int) fitness(child1));
        fitVals.push_back((int) fitness(child2));
    }
    return fitVals;
}

int crossoverCheck(ofstream& outp){
    outp << "Population Fitness Values:" << endl;

    for (int i = 0; i < popsize; ++i) {
        outp << fits[i] << endl;
    }

    outp << "Population Crossover Checks:" << endl;

    for (int mom = 0; mom < popsize; mom++) {
        for (int dad = mom + 1; dad < popsize; dad++) {
            outp << "Parent Idxs: " << mom << "\t" << dad << endl;
            outp << "Parent Fits: " << fits[mom] << "\t" << fits[dad] << endl;
            printVector<ofstream, int>(outp, runMultiMating(mateTests, pop[mom], pop[dad]), "", "\n", true);
        }
    }
    return 0;
}

int main(int argc, char *argv[]) {
    double expBestFit = (BIGGERBETTER ? 0 : MAXFLOAT);
    SDA expBestSDA = SDA(sdaStates, numChars);

    char filename[100];
    ofstream runStats, expStats, readMe, crossStart, crossEnd;

    vector<double> bests;
    bests.reserve(runs);

    getArgs(argv);
    initAlg();
    cmdLineIntro(cout);
    sprintf(pathToOut, "./AAMOut/AAMatchTest on Seq%d with %04dPop, %02dSta, %02dMut, %02dTsz, %dCross/",
            seqNum, popsize, sdaStates, maxMuts, tournSize, crossOp);
    filesystem::create_directory(pathToOut);
    expStats.open(string(pathToOut) + "./exp.dat", ios::out);
    readMe.open(string(pathToOut) + "./read.me", ios::out);
    makeReadMe(readMe);
    readMe.close();

    int tmp;
    for (int run = 1; run < runs + 1; ++run) {
        initPop(run);
        sprintf(filename, "%scrossover%02d_start.dat", pathToOut, run);
        crossStart.open(filename, ios::out);
        crossoverCheck(crossStart);
        crossStart.close();
        cout<<"Crossover Check Start Complete!"<<endl;

        sprintf(filename, "%srun%02d.dat", pathToOut, run);
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

        sprintf(filename, "%scrossover%02d_end.dat", pathToOut, run);
        crossEnd.open(filename, ios::out);
        crossoverCheck(crossEnd);
        crossEnd.close();
        cout<<"Crossover Check End Complete!"<<endl;
    }

    ofstream best;
    best.open(string(pathToOut) + "./best.dat", ios::out);
    expReport(best, bests, expBestSDA, BIGGERBETTER);
    best.close();
    delete[] pop;
    cout << "Program Completed Successfully!" << endl;
    return 0;
}