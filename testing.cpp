#include "testing.h"

/**
 * This method mirrors the main method in main.cpp.  For command line argument descriptions see that method.
 *
 * @param argc
 * @param argv
 * @return
 */

int main(int argc, char *argv[]) {
    getArgs(argv);
    string pathToSeqs = "./Sequences.dat";
    char filename[200];
    ofstream runStats, expStats, readMe, crossFile, mutateFile, sdaFile, runGains, runPopulation;

    vector<double> bests;
    bests.reserve(runs);
    double expBestFit = (BIGGER_BETTER ? 0 : MAXFLOAT);
    SDA expBestSDA = SDA(sdaStates, numChars, 2, seqLen);

    initAlg(pathToSeqs);
    cmdLineIntro(cout);
    char dynamicMessage[20];
    sprintf(dynamicMessage, "%s%d", (dynamicMutOperator == 0 ? "Static" : "Dynamic"), dynamicMutOperator);
    sprintf(pathToOut, "./AAMTestOut/AAMatch on Seq%d with %.1fmilMMEs, %04dPS, %02dSt, %02dNTM, %02dNRM, %s,"
                       " %dTS, %sCO, %03d%%CrR, %03d%%MR, %03d%%CuR, %sCu/", seqNum, (double) maxGens / 1000000,
            popsize, sdaStates, numTransMuts, numRespMuts, dynamicMessage, tournSize,
            (crossoverOp == 0 ? "2Pt" : "1St"), (int) (crossoverRate * 100), (int) (mutationRate * 100),
            (int) (cullingRate * 100), (randomCulling ? "Rand" : "Worst"));
    mkdir(pathToOut, 0777);
    sprintf(filename, "%sCrossover Checks/", pathToOut);
    mkdir(filename, 0777);
    sprintf(filename, "%sMutate Checks/", pathToOut);
    mkdir(filename, 0777);
    sprintf(filename, "%sSDA Checks/", pathToOut);
    mkdir(filename, 0777);
    expStats.open(string(pathToOut) + "./exp.dat", ios::out);
    readMe.open(string(pathToOut) + "./read.me", ios::out);
    makeReadMe(readMe);
    readMe.close();

    int tmp;
    for (int run = 1; run < runs + 1; ++run) {
        initPop(run);
        sprintf(filename, "%s/Crossover Checks/crossover%02d_%05dk.dat", pathToOut, run, 0);
        crossFile.open(filename, ios::out);
        crossoverCheck(crossFile);
        crossFile.close();
        sprintf(filename, "%s/Mutate Checks/mutate%02d_%05dk.dat", pathToOut, run, 0);
        mutateFile.open(filename, ios::out);
        mutateCheck(mutateFile);
        mutateFile.close();

        cout << "Initial Checks Complete!" << endl;
        sprintf(filename, "%s/SDA Checks/pop%02d.dat", pathToOut, run);
        sdaFile.open(filename, ios::out);
        sprintf(filename, "%srun%02d.dat", pathToOut, run);
        runStats.open(filename, ios::out);
        sprintf(filename, "%sgains%02d.dat", pathToOut, run);
        runGains.open(filename, ios::out);

        printExpStatsHeader(cout);
        printExpStatsHeader(runStats);
        sdaCheck(sdaFile, 0);
        report(runStats, run, 0, BIGGER_BETTER);
        int gen = 1;
        int stallCount = 0;
        double best = (BIGGER_BETTER ? 0 : MAXFLOAT);
        while (gen <= maxGens && stallCount < TERM_CRIT) {
            matingEvent(BIGGER_BETTER, gen, runGains);

            if (gen % REPORT_EVERY == 0) {
                tmp = report(runStats, run, (int) gen / (REPORT_EVERY), BIGGER_BETTER);
                if ((BIGGER_BETTER && tmp > best) || (!BIGGER_BETTER && tmp < best)) {
                    best = tmp;
                    stallCount = 0;
                } else {
                    stallCount++;0
                }
            }

            if (gen % (int) (CULLING_EVERY * REPORT_EVERY) == 0 && stallCount < TERM_CRIT) {
                sprintf(filename, "%s/Crossover Checks/crossover%02d_%05dk.dat", pathToOut, run, gen / 1000);
                crossFile.open(filename, ios::out);
                crossoverCheck(crossFile);
                crossFile.close();
                sprintf(filename, "%s/Mutate Checks/mutate%02d_%05dk.dat", pathToOut, run, gen / 1000);
                mutateFile.open(filename, ios::out);
                mutateCheck(mutateFile);
                mutateFile.close();
                sdaCheck(sdaFile, gen);
                if (gen != maxGens){
                    culling(cullingRate, randomCulling, BIGGER_BETTER);
                }
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
        runGains.close();

        sprintf(filename, "%s/Crossover Checks/crossover%02d_%05dk.dat", pathToOut, run, (gen - 1) / 1000);
        crossFile.open(filename, ios::out);
        crossoverCheck(crossFile);
        crossFile.close();
        sprintf(filename, "%s/Mutate Checks/mutate%02d_%05dk.dat", pathToOut, run, (gen - 1) / 1000);
        mutateFile.open(filename, ios::out);
        mutateCheck(mutateFile);
        mutateFile.close();
        sdaCheck(sdaFile, gen);
        sdaFile.close();
        cout << "Final Checks Complete!" << endl;
    }

    ofstream best;
    best.open(string(pathToOut) + "./best.dat", ios::out);
    expReport(best, bests, expBestSDA, BIGGER_BETTER);
    best.close();
    delete[] pop;
    cout << "Program Completed Successfully!" << endl;
    return 0;
}