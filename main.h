#include <iomanip>
#include <algorithm>
#include <vector>
#include <cmath>
#include <fstream>
#include <sys/stat.h>
#include "SDA.h"

using namespace std;

// System Parameters
#define REPORT_EVERY (int)10000
#define TERM_CRIT 50
#define CULLING_EVERY 1
#define BIGGER_BETTER (bool)true

// Experiment Parameters
int popsize;
int numChars;
int sdaStates;
double *sdaProbs;
int seed;
int runs;
int maxGens;
int numTransMuts;
int numRespMuts;
int dynamicMutOperator;
int tournSize;
int seqNum;
int crossoverOp;
double crossoverRate;
double mutationRate;
double cullingRate;
bool randomCulling;
int populationBestIdx;
double populationBestFit;

vector<int> goalSeq;
vector<int> testSeq;
vector<char> charSeq;
int seqLen;

SDA *pop;
vector<double> fits;

char pathToOut[150];

// Program Method Declarations:
int getArgs(char *arguments[]);
int initAlg(const string &pathToSeqs);
vector<vector<int>> getSequences(const string &pathToSeqs);
int cmdLineIntro(ostream &outp);
int makeReadMe(ostream &outp);
int initPop(int run);
double fitness(SDA &sda);
int printExpStatsHeader(ostream &outp);
double report(ofstream &outp, int run, int rptNum, bool biggerBetter);
template<class T>
vector<double> calcStats(vector<T> vals, bool biggerBetter);
int matingEvent(bool biggerBetter);
vector<int> tournSelect(int size, bool decreasing);
bool compareFitness(int popIdx1, int popIdx2);
int culling(double percentage, bool rndPick, bool biggerBetter);
int runReport(ostream &outp, bool biggerBetter);
template<class T1, class T2>
int printMatches(T1 &outp, const vector<T2> &test, const vector<T2> &goal, bool newline);
int expReport(ostream &outp, vector<double> bestFits, SDA bestSDA, bool biggerBetter);

// Helper Method Declarations:
vector<int> seqToVector(const string &seq);
int intToChar(const vector<int> &from, vector<char> &to);
template<class T1, class T2>
int printVector(T1 &outp, vector<T2> vec, const string &msg, const string &sep, bool newline);
template<class T1, class T2>
int printIdxsOfVector(T1 &outp, vector<T2> vec, const vector<int> &idxs, const string &msg, const string &sep,
                      bool newline);

// Helper Class
class multiStream : public ostream {
public:
    multiStream(ostream &os1, ostream &os2) : os1(os1), os2(os2) {}

    template<class T>
    multiStream &operator<<(const T &x) {
        os1 << x;
        os2 << x;
        return *this;
    }

private:
    ostream &os1;
    ostream &os2;
};

// Method Definitions
/**
 * This method collects the command line arguments and places them in the respective variable.
 *
 * @param arguments popsize, numChars, sdaStates, seed, runs, maxGens, maxMuts, seqNum, tournSize, crossoverOp,
 *                  crossoverRate, mutationRate, cullingRate, randomCulling
 * @return
 */
int getArgs(char *arguments[]) {
    size_t pos;
    string arg;
    arg = arguments[1]; // popsize
    popsize = stoi(arg, &pos);
    arg = arguments[2]; // numChars
    numChars = stoi(arg, &pos);
    arg = arguments[3]; // sdaStates
    sdaStates = stoi(arg, &pos);
    arg = arguments[4]; // seed
    seed = stoi(arg, &pos);
    arg = arguments[5]; // runs
    runs = stoi(arg, &pos);
    arg = arguments[6]; // maxGens
    maxGens = stoi(arg, &pos);
    arg = arguments[7]; // numTransMuts
    numTransMuts = stoi(arg, &pos);
    arg = arguments[8]; // numRespMuts
    numRespMuts = stoi(arg, &pos);
    arg = arguments[9]; // dynamicMutOperator
    dynamicMutOperator = stoi(arg, &pos);
    arg = arguments[10]; // seqNum
    seqNum = stoi(arg, &pos);
    arg = arguments[11]; // tournSize
    tournSize = stoi(arg, &pos);
    arg = arguments[12]; // crossoverOp
    crossoverOp = stoi(arg, &pos);
    arg = arguments[13]; // crossoverRate
    crossoverRate = stod(arg, &pos);
    arg = arguments[14]; // mutationRate
    mutationRate = stod(arg, &pos);
    arg = arguments[15]; // cullingRate
    cullingRate = stod(arg, &pos);
    arg = arguments[16]; // randomCulling
    randomCulling = stoi(arg, &pos) == 1;
    cout << "Arguments Captured!" << endl;
    return 0;
}

int initAlg(const string &pathToSeqs) {
    srand48(time(nullptr)); // use system time as random number seed
    srand48(seed);           // read the random number seed
    goalSeq = getSequences(pathToSeqs)[seqNum];
    seqLen = (int) goalSeq.size();
    fits.reserve(popsize);

    pop = new SDA[popsize];
    for (int idx = 0; idx < popsize; ++idx) {
        pop[idx] = SDA(sdaStates, numChars, 2, seqLen);
    }

    testSeq.reserve(seqLen);
    charSeq.reserve(seqLen);
    for (int idx = 0; idx < seqLen; ++idx) {
        testSeq.push_back(-1);
        charSeq.push_back('a');
    }

    populationBestIdx = -1;
    populationBestFit = (BIGGER_BETTER ? 0.0 : MAXFLOAT);
    return 0;
}

vector<vector<int>> getSequences(const string &pathToSeqs) {
    string tmp;
    ifstream in(pathToSeqs);
    vector<vector<int>> rtn;
    for (int seqIdx = 0; seqIdx < 6; ++seqIdx) {
        if (in.is_open()) {
            getline(in, tmp);
            getline(in, tmp);
            getline(in, tmp);
            rtn.push_back(seqToVector(tmp));
            getline(in, tmp);
        }
    }
    in.close();
    return rtn;
}

int cmdLineIntro(ostream &outp) {
    outp << "Amino Acid/DNA Sequence Matcher." << endl;
    outp << "Fitness Function: Naive." << endl;
    outp << "Matching Sequence (length " << seqLen << "): ";
    intToChar(goalSeq, charSeq);
    printVector<ostream, char>(outp, charSeq, "", "", true);
    outp << "See read.me for system and experiment parameters." << endl;
    outp << endl;
    return 0;
}

int makeReadMe(ostream &outp) {
    cmdLineIntro(outp);
    outp << "Experiment Parameters" << endl;
    outp << "Population Size: " << popsize << endl;
    outp << "Number of States: " << sdaStates << endl;
    outp << "Alphabet Size: " << numChars << endl;
    outp << "Max Generations: " << maxGens << endl;
    outp << "Report Every: " << REPORT_EVERY << " generations" << endl;
    outp << "Termination Criteria: maximum fitness or no improvement for " << TERM_CRIT * REPORT_EVERY
         << " generations" << endl;
    outp << "Tournament Size: " << tournSize << endl;
    outp << "Crossover Operator: " << (crossoverOp == 0 ? "Two-Point Crossover" : "One-State Crossover") << endl;
    outp << "Crossover Rate: " << (int) (crossoverRate * 100) << "%" << endl;
//    outp << "Maximum Number of Mutations: " << maxMuts << endl; //TODO: Update
    outp << "Mutation Rate: " << (int) (mutationRate * 100) << "%" << endl;
    outp << "Culling Every: " << CULLING_EVERY * REPORT_EVERY << " generations" << endl;
    if (randomCulling) {
        outp << "Culling Winners: Random " << (int) (cullingRate * 100) << "% of the population" << endl;
    } else {
        outp << "Culling Winners: Worst " << (int) (cullingRate * 100) << "% of the population" << endl;
    }
    return 0;
}

int initPop(int run) {
    if (runs > 1) {
        cout << "Beginning Run " << run << " of " << runs << endl;
    }
    fits.clear();
    for (int idx = 0; idx < popsize; ++idx) {
        pop[idx].randomize();
        fits.push_back(fitness(pop[idx]));
    }
    cout << "Population Generated!" << endl;
    return 0;
}

double fitness(SDA &sda) {
    double val = 0.0;
    sda.fillOutput(testSeq);
    for (int i = 0; i < seqLen; ++i) {
        if (testSeq[i] == goalSeq[i]) {
            val += 1;
        }
    }
    return val;
}

int printExpStatsHeader(ostream &outp) {
    outp << left << setw(5) << "Run";
    outp << left << setw(4) << "RI";
    outp << left << setw(10) << "Mean";
    outp << left << setw(12) << "95% CI";
    outp << left << setw(10) << "SD";
    outp << left << setw(8) << "Best";
    outp << left << setw(10) << "% Correct";
    outp << endl;
    return 0;
}

double report(ofstream &outp, int run, int rptNum, bool biggerBetter) {
    vector<double> stats = calcStats<double>(fits, biggerBetter); // {mean, stdDev, 95CI, best}
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

template<class T>
vector<double> calcStats(vector<T> vals, bool biggerBetter) {
    double sum = 0.0;
    double bestVal = (biggerBetter ? 0.0 : MAXFLOAT);

    int val;
    for (int idx = 0; idx < vals.size(); ++idx) {
        val = vals[idx];
        sum += val;
        if ((biggerBetter && val > bestVal) || (!biggerBetter && val < bestVal)) {
            bestVal = val;
            populationBestIdx = idx;
            populationBestFit = bestVal;
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

int matingEvent(bool biggerBetter) {
    int numMuts;
    SDA child1, child2;

    vector<int> idxs = tournSelect(tournSize, biggerBetter);

    child1.copy(pop[idxs[0]]);
    child2.copy(pop[idxs[1]]);
    if (drand48() < crossoverRate) {
        if (crossoverOp == 0) child1.twoPointCrossover(child2);
        else if (crossoverOp == 1) child1.oneStateCrossover(child2);
    }

    if (drand48() < mutationRate) {
        if (dynamicMutOperator == 0){
            child1.mutate(numTransMuts, numRespMuts);
            child2.mutate(numTransMuts, numRespMuts);
        }
    }

    pop[idxs.end()[-1]] = child1;
    pop[idxs.end()[-2]] = child2;

    fits[idxs.end()[-1]] = fitness(child1);
    fits[idxs.end()[-2]] = fitness(child2);
    return 0;
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
        } while (tournIdxs.size() < size);
    }

    sort(tournIdxs.begin(), tournIdxs.end(), compareFitness);
    if (decreasing) {
        reverse(tournIdxs.begin(), tournIdxs.end());
    }
    return tournIdxs;
}

bool compareFitness(int popIdx1, int popIdx2) {
    return fits[popIdx1] < fits[popIdx2];
}

int culling(double percentage, bool rndPick, bool biggerBetter) {
    if (percentage == 1) { // Cull all but the best SDA
        pop[0].copy(pop[populationBestIdx]);
        fits[0] = fitness(pop[0]);
        for (int idx = 1; idx < popsize; ++idx) {
            pop[idx].randomize();
            fits[idx] = fitness(pop[idx]);
        }
        return 0;
    }

    // Otherwise, cull percentage% of the population
    int numKillings = (int) (popsize * percentage);
    vector<int> contestants;
    if (rndPick) { // Cull random percentage% of the population
        contestants.reserve(numKillings);
        contestants = tournSelect(numKillings, !biggerBetter);
    } else { // Cull worst percentage% of the population
        contestants.reserve(popsize);
        contestants = tournSelect(popsize, !biggerBetter);
    }

    int idxToCull;
    for (int cnt = 0; cnt < numKillings; cnt++) {
        idxToCull = contestants[cnt];
        pop[idxToCull].randomize();
        fits[idxToCull] = fitness(pop[idxToCull]);
    }
    return 0;
}

int runReport(ostream &outp, bool biggerBetter) {
    auto maxIterator = minmax_element(fits.begin(), fits.end());
    int bestIdx;
    if (biggerBetter) {
        bestIdx = (int) distance(fits.begin(), maxIterator.second);
    } else {
        bestIdx = (int) distance(fits.begin(), maxIterator.first);
    }

    multiStream printAndSave(cout, outp);
    printAndSave << "The best fitness is " << fits[bestIdx] << "\n";
    pop[bestIdx].fillOutput(testSeq);
    printAndSave << left << setw(20) << "Best Match: ";
    intToChar(testSeq, charSeq);
    printVector<multiStream, char>(printAndSave, charSeq, "", "", true);
    printAndSave << left << setw(20) << "Matches: ";
    printMatches<multiStream, int>(printAndSave, testSeq, goalSeq, true);
    printAndSave << left << setw(20) << "Desired Sequence: ";
    intToChar(goalSeq, charSeq);
    printVector<multiStream, char>(printAndSave, charSeq, "", "", true);

    outp << "SDA" << endl;
    pop[bestIdx].printSDA(outp);

    vector<double> fitsCopy = fits;
    sort(fitsCopy.begin(), fitsCopy.end());
    if (biggerBetter) {
        reverse(fitsCopy.begin(), fitsCopy.end());
    }
    printVector<multiStream, double>(printAndSave, fitsCopy, "Fitness Values: ", " ", true);
    cout << endl;
    return bestIdx;
}

template<class T1, class T2>
int printMatches(T1 &outp, const vector<T2> &test, const vector<T2> &goal, bool newline) {
    int count = 0;
    for (int idx = 0; idx < goal.size(); idx++) {
        if (test[idx] == goal[idx]) {
            outp << "X";
            count++;
        } else {
            outp << " ";
        }
    }
    if (newline) outp << "\n";
    return count;
}

int expReport(ostream &outp, vector<double> bestFits, SDA bestSDA, bool biggerBetter) {
    vector<double> stats = calcStats<double>(bestFits, BIGGER_BETTER);
    multiStream printAndSave(cout, outp);
    printAndSave << "Experiment Report:" << "\n";
    printAndSave << "Best Run: " << populationBestIdx + 1 << "\n";
    printAndSave << "Best Fitness: " << bestFits[populationBestIdx] << " of " << seqLen << "\n";
    printAndSave << "Fitness 95% CI: " << stats[0] << " +- " << stats[2] << "\n";
    printAndSave << "\n";
    printAndSave << left << setw(20) << "Best Match: ";
    bestSDA.fillOutput(testSeq);
    intToChar(testSeq, charSeq);
    printVector<multiStream, char>(printAndSave, charSeq, "", "", true);
    printAndSave << left << setw(20) << "Matches: ";
    printMatches<multiStream, int>(printAndSave, testSeq, goalSeq, true);
    printAndSave << left << setw(20) << "Desired Sequence: ";
    intToChar(goalSeq, charSeq);
    printVector<multiStream, char>(printAndSave, charSeq, "", "", true);
    printAndSave << "\nSDA:\n";
    bestSDA.printSDA(cout);
    bestSDA.printSDA(outp);
    return 0;
}

vector<int> seqToVector(const string &seq) {
    vector<int> sequence;
    for (char c: seq) {
        if (c == 'g' || c == 'G') {
            sequence.push_back(0);
        } else if (c == 'c' || c == 'C') {
            sequence.push_back(1);
        } else if (c == 'a' || c == 'A') {
            sequence.push_back(2);
        } else if (c == 't' || c == 'T') {
            sequence.push_back(3);
        }
    }
    return sequence;
}

int intToChar(const vector<int> &from, vector<char> &to) {
    for (int idx = 0; idx < seqLen; ++idx) {
        switch (from[idx]) {
            case 0:
                to[idx] = 'G';
                break;
            case 1:
                to[idx] = 'C';
                break;
            case 2:
                to[idx] = 'A';
                break;
            case 3:
                to[idx] = 'T';
                break;
        }
    }
    return 0;
}

template<class T1, class T2>
int printVector(T1 &outp, vector<T2> vec, const string &msg, const string &sep, bool newline) {
    outp << msg;
    for (auto val: vec) {
        outp << val << sep;
    }
    if (newline) outp << "\n";
    return 0;
}

template<class T1, class T2>
int printIdxsOfVector(T1 &outp, vector<T2> vec, const vector<int> &idxs, const string &msg, const string &sep,
                      bool newline) {
    outp << msg;
    for (auto idx: idxs) {
        outp << vec[idx] << sep;
    }
    if (newline) outp << "\n";
    return 0;
}


