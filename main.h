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
#define TERM_CRIT 250
int CULLING_EVERY;
#define BIGGER_BETTER (bool)true
double MIN_GEN_RATIO = 0.5;

// Experiment Parameters
int popsize;
int numChars;
int sdaStates;
double *sdaProbs;
int seed;
int runs;
int maxGens;
int initRun;
int numRuns;
int initNumTransMuts;
int initNumRespMuts;
int curNumTransMuts;
int curNumRespMuts;
int upBoundMuts;
int dynamicMutOperator;
int tournSize;
int seqNum;
int crossoverOp;
double crossoverRate;
double mutationRate;
double cullingRate;
bool randomCulling;
int popBestIdx;
//pair<int, int> populationBestFit;// current best population fitness
//double prevBestFitness = 0;
int RICounter;// Report Interval counter

vector<int> goalSeq;
vector<int> testSeq;
vector<char> charSeq;
int seqLen;

SDA *pop;
vector<double> matchFits;
vector<int> noveltyFits;

char pathToOut[150];

// Program Method Declarations:
int getArgs(int numArgs, char *arguments[]);
int initAlg(const string &pathToSeqs);
vector<vector<int>> getSequences(const string &pathToSeqs);
int cmdLineIntro(ostream &outp);
int makeReadMe(ostream &outp);
int initPop(int run);
double calcMatchFit(SDA &sda);
int calcNoveltyFit(int idx);
double updateNoveltyFit(int forIdx, vector<int> cmprStr, int score);
int printExpStatsHeader(ostream &outp);
pair<double, double> report(ofstream &outp, int run, int rptNum, bool biggerBetter);
template<class T>
vector<double> calcStats(vector<T> vals, bool biggerBetter);
int matingEvent(bool biggerBetter);
vector<int> tournSelect(int size, bool decreasing);
bool compareFitness(int popIdx1, int popIdx2);
int culling(double percentage, bool rndPick, bool biggerBetter);
int runReport(ostream &outp, bool biggerBetter, pair<int, int> bestFit);
template<class T1, class T2>
int printMatches(T1 &outp, const vector<T2> &test, const vector<T2> &goal, bool newline);
int expReport(ostream &outp, vector<double> bestFits, SDA bestSDA, bool biggerBetter);
int sdaCheck(ofstream &outStream, int currentGen);

// Helper Method Declarations:
vector<int> seqToVector(const string &seq);
int intToChar(const vector<int> &from, vector<char> &to);
template<class T1, class T2>
int printVector(T1 &outp, vector<T2> vec, const string &msg, const string &sep, bool newline);
template<class T1, class T2>
int printIdxsOfVector(T1 &outp, vector<T2> vec, const vector<int> &idxs, const string &msg, const string &sep,
                      bool newline);
int printPopFits(ostream &outStrm);

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
int getArgs(int numArgs, char *arguments[]) {
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
    arg = arguments[7]; // initNumTransMut
    initNumTransMuts = stoi(arg, &pos);
    arg = arguments[8]; // initNumRespMut
    initNumRespMuts = stoi(arg, &pos);
    arg = arguments[9]; // dynamicMutOperator
    dynamicMutOperator = stoi(arg, &pos);
    arg = arguments[10]; // upBoundMuts
    upBoundMuts = stoi(arg, &pos);
    arg = arguments[11]; // seqNum
    seqNum = stoi(arg, &pos);
    arg = arguments[12]; // tournSize
    tournSize = stoi(arg, &pos);
    arg = arguments[13]; // crossoverOp
    crossoverOp = stoi(arg, &pos);
    arg = arguments[14]; // crossoverRate
    crossoverRate = stod(arg, &pos);
    arg = arguments[15]; // mutationRate
    mutationRate = stod(arg, &pos);
    arg = arguments[16]; // cullingRate
    cullingRate = stod(arg, &pos);
    arg = arguments[17]; // randomCulling
    randomCulling = stoi(arg, &pos) == 1;
    arg = arguments[18];
    CULLING_EVERY = stoi(arg, &pos);
    if (numArgs > 19){
        arg = arguments[19];
        initRun = stoi(arg, &pos);
    }
    cout << "Arguments Captured!" << endl;
    return 0;
}

int initAlg(const string &pathToSeqs) {
    srand48(time(nullptr)); // use system time as random number seed
    srand48(seed);           // read the random number seed
    //vector<vector<int>> sequences = getSequences(pathToSeqs);
    goalSeq = getSequences(pathToSeqs)[seqNum];
    seqLen = (int) goalSeq.size();
    matchFits.reserve(popsize);

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

    popBestIdx = -1;
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
    outp << "Mutation Rate: " << (int) (mutationRate * 100) << "%" << endl;
    outp << "Default Number of Transition Mutations: " << initNumTransMuts << endl;
    outp << "Default Number of Response Mutations: " << initNumRespMuts << endl;
    if (dynamicMutOperator == 0) {
        outp << "Static Number of Mutations" << endl;
    } else {
        outp << "Dynamic Number of Mutations Using Version: " << dynamicMutOperator << endl;
    }
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
        cout << "Beginning Run " << run << " of " << initRun + runs << endl;
    }
    matchFits.clear();
    for (int idx = 0; idx < popsize; ++idx) {
        pop[idx].randomize();
        matchFits.push_back(calcMatchFit(pop[idx]));
        noveltyFits.push_back(-1);
    }

    for (int idx = 0; idx < popsize; idx++) {
        noveltyFits[idx] = calcNoveltyFit(idx);
    }

    cout << "Population Generated!" << endl;
    return 0;
}

double calcMatchFit(SDA &sda) {
    double val = 0.0;
    sda.fillOutput(testSeq);
    for (int i = 0; i < seqLen; ++i) {
        if (testSeq[i] == goalSeq[i]) {
            val += 1;
        }
    }
    return val;
}

int calcNoveltyFit(int idx) {
    int val = 0;

    vector<int> idx_same_match_fit;
    idx_same_match_fit.reserve(popsize);
    for (int i = 0; i < popsize; ++i) {
        if (i != idx && matchFits[i] == matchFits[idx]) {
            idx_same_match_fit.push_back(i);
        }
    }

    vector<vector<int>> strings_same_match_fit;
    strings_same_match_fit.reserve(idx_same_match_fit.size());

    for (int i: idx_same_match_fit) {
        pop[i].fillOutput(testSeq);
        strings_same_match_fit.push_back(testSeq);
    }

    pop[idx].fillOutput(testSeq);

    for (int i = 0; i < seqLen; i++) {
        if (testSeq[i] != goalSeq[i]) {
            for (vector<int> vec: strings_same_match_fit) {
                if (testSeq[i] == vec[i]) {
                    val++;
                }
            }
        }
    }

    return val;
}

double updateNoveltyFit(int forIdx, vector<int> cmprStr, int score) {
    pop[forIdx].fillOutput(testSeq);
    int val = noveltyFits[forIdx];

    for (int i = 0; i < testSeq.size(); i++) {
        if (testSeq[i] != goalSeq[i]) {
            if (testSeq[i] == cmprStr[i]) {
                val += score;
            }
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
    outp << left << setw(12) << "Best Match";
    outp << left << setw(12) << "Best Nov.";
    outp << left << setw(10) << "% Correct";
    outp << endl;
    return 0;
}

pair<double, double> report(ofstream &outp, int run, int rptNum, bool biggerBetter) {
    vector<double> stats = calcStats<double>(matchFits, biggerBetter); // {mean, stdDev, 95CI, bestVal, bestIdx}
    int bestMatchFit = (int) stats[3];
    vector<int> noveltyFitsSameMatchFit;
    vector<int> noveltyIdxs;
    noveltyFitsSameMatchFit.reserve(popsize);
    noveltyIdxs.reserve(popsize);

    for (int idx = 0; idx < popsize; idx++) {
        if (matchFits[idx] == bestMatchFit) {
            noveltyFitsSameMatchFit.push_back(noveltyFits[idx]);
            noveltyIdxs.push_back(idx);
        }
    }

    vector<double> stats2 = calcStats<int>(noveltyFitsSameMatchFit, false); // {mean, stdDev, 95CI, bestVal, bestIdx}

    popBestIdx = noveltyIdxs[(int) stats2[4]];

    multiStream printAndSave(cout, outp);

    printAndSave << left << setw(5) << run;
    printAndSave << left << setw(4) << rptNum;
    printAndSave << left << setw(10) << stats[0];
    printAndSave << left << setw(12) << stats[2];
    printAndSave << left << setw(10) << stats[1];
    printAndSave << left << setw(12) << stats[3];
    printAndSave << left << setw(12) << stats2[3];
    printAndSave << left << setw(7) << (stats[3] / seqLen) * 100 << "%";
    printAndSave << "\n";
    return make_pair(stats[3], stats2[3]);
}

template<class T>
vector<double> calcStats(vector<T> vals, bool biggerBetter) {
    double sum = 0.0;
    double bestVal = (biggerBetter ? 0.0 : MAXFLOAT);
    int bestIdx = 0;

    int val;
    for (int idx = 0; idx < vals.size(); ++idx) {
        val = vals[idx];
        sum += val;
        if ((biggerBetter && val > bestVal) || (!biggerBetter && val < bestVal)) {
            bestVal = val;
            bestIdx = idx;
        }
    }

    if (vals.size() == 1){
        return {bestVal, 0.0, 0.0, bestVal, (double) bestIdx}; // {mean, stdDev, 95CI, bestVal, bestIdx}
    }

    double mean = sum / (double) vals.size();
    double stdDevSum = 0.0;
    for (int val: vals) {
        stdDevSum += pow((double) val - mean, 2);
    }
    double stdDev = sqrt(stdDevSum / ((double) vals.size() - 1.0));
    double CI95 = 1.96 * (stdDev / sqrt(vals.size()));

    return {mean, stdDev, CI95, bestVal, (double) bestIdx}; // {mean, stdDev, 95CI, bestVal, bestIdx}
}

int matingEvent(bool biggerBetter) {
    int p1idx, p2idx;
    int c1idx, c2idx;
    double oldMatchFit1, oldMatchFit2;
    SDA child1, child2;
    SDA old1, old2;

    vector<int> idxs = tournSelect(tournSize, biggerBetter);

//    cout<< "Tourn: ";
//    for (int i = 0; i < tournSize; i++){
//        cout << "[" << matchFits[i] << ", " << noveltyFits[i] << "], ";
//    }
//    cout << endl;

    p1idx = idxs[0];
    p2idx = idxs[1];
    c1idx = idxs.end()[-1];
    c2idx = idxs.end()[-2];

    old1.copy(pop[c1idx]);
    old2.copy(pop[c2idx]);
    oldMatchFit1 = matchFits[c1idx];
    oldMatchFit2 = matchFits[c2idx];

    child1.copy(pop[p1idx]);
    child2.copy(pop[p2idx]);

    if (drand48() < crossoverRate) {
        if (crossoverOp == 0) child1.twoPointCrossover(child2);
        else if (crossoverOp == 1) child1.oneStateCrossover(child2);
    }

    if (drand48() < mutationRate) {
        if (dynamicMutOperator == 0) {
            child1.mutate(curNumTransMuts, curNumRespMuts);
            child2.mutate(curNumTransMuts, curNumRespMuts);
        }
    }

    pop[c1idx] = child1;
    pop[c2idx] = child2;

    matchFits[c1idx] = calcMatchFit(child1);
    matchFits[c2idx] = calcMatchFit(child2);

    noveltyFits[c1idx] = calcNoveltyFit(c1idx);
    noveltyFits[c2idx] = calcNoveltyFit(c2idx);

    // Update novelty for all solutions that matched the parent's matching fitness
    for (int i = 0; i < popsize; i++) {
        if (i != c1idx && i != c2idx) {
            if (matchFits[i] == oldMatchFit1) { // Same as old string, remove score
                old1.fillOutput(testSeq);
                noveltyFits[i] = updateNoveltyFit(i, testSeq, -1);
            }
            if (matchFits[i] == matchFits[c1idx]) { // Same as new string, add score
                child1.fillOutput(testSeq);
                noveltyFits[i] = updateNoveltyFit(i, testSeq, 1);
            }
            if (matchFits[i] == oldMatchFit2) { // Same as old string, remove score
                old2.fillOutput(testSeq);
                noveltyFits[i] = updateNoveltyFit(i, testSeq, -1);
            }
            if (matchFits[i] == matchFits[c2idx]) { // Same as new string, add score
                child2.fillOutput(testSeq);
                noveltyFits[i] = updateNoveltyFit(i, testSeq, 1);
            }
        }
    }

//    for (int i = 0; i < popsize; i++){
//        noveltyFits[i] = calcNoveltyFit(i);
//    }
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
    if (matchFits[popIdx1] < matchFits[popIdx2]) {
        return true;
    }
    if (matchFits[popIdx1] > matchFits[popIdx2]) {
        return false;
    }
    if (matchFits[popIdx1] == matchFits[popIdx2]) {
        if (noveltyFits[popIdx1] > noveltyFits[popIdx2]) {
            return true;
        } else {
            return false;
        }
    }
    cout << "ERROR: compare fitness not working as intended!" << endl;
    return false;
}

int culling(double percentage, bool rndPick, bool biggerBetter) {
    if (percentage == 1) { // Cull all but the best SDA
        pop[0].copy(pop[popBestIdx]);
        matchFits[0] = calcMatchFit(pop[0]);
        for (int idx = 1; idx < popsize; ++idx) {
            pop[idx].randomize();
            matchFits[idx] = calcMatchFit(pop[idx]);
        }

        for (int idx = 0; idx < popsize; idx++) {
            noveltyFits[idx] = calcNoveltyFit(idx);
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
        if (idxToCull != popBestIdx) {
            pop[idxToCull].randomize();
            matchFits[idxToCull] = calcMatchFit(pop[idxToCull]);
        }
    }

    for (int idx = 0; idx < popsize; idx++) {
        noveltyFits[idx] = calcNoveltyFit(idx);
    }
    return 0;
}

int runReport(ostream &outp, bool biggerBetter, pair<int, int> bestFit) {
//    auto maxIterator = minmax_element(matchFits.begin(), matchFits.end());
//    int bestIdx;
//    if (biggerBetter) {
//        bestIdx = (int) distance(matchFits.begin(), maxIterator.second);
//    } else {
//        bestIdx = (int) distance(matchFits.begin(), maxIterator.first);
//    }

    multiStream printAndSave(cout, outp);
//    printAndSave << "The best fitness is " << matchFits[bestIdx] << "\n";
    printAndSave << "The best fitness is " << bestFit.first;
    printAndSave << " matches with novelty " << bestFit.second << "\n";
    pop[popBestIdx].fillOutput(testSeq);
    printAndSave << left << setw(20) << "Best Match: ";
    intToChar(testSeq, charSeq);
    printVector<multiStream, char>(printAndSave, charSeq, "", "", true);
    printAndSave << left << setw(20) << "Matches: ";
    printMatches<multiStream, int>(printAndSave, testSeq, goalSeq, true);
    printAndSave << left << setw(20) << "Desired Sequence: ";
    intToChar(goalSeq, charSeq);
    printVector<multiStream, char>(printAndSave, charSeq, "", "", true);

    outp << "SDA" << endl;
    pop[popBestIdx].printSDA(outp);

    vector<int> sortedIdxs = tournSelect(popsize, BIGGER_BETTER);
    printAndSave << "Fitness Values: ";

    bool first = true;
    for (int idx: sortedIdxs) {
        if (!first) {
            printAndSave << ", ";
        }
        printAndSave << "[" << matchFits[idx] << ", " << noveltyFits[idx] << "]";
        first = false;
    }
    printAndSave << "\n";
    return 0;
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

int expReport(ostream &outp, vector<pair<int, int>> bestFits, SDA bestSDA, bool biggerBetter) {
    vector<int> bestMatchFits, bestNoveltyFits;
    bestMatchFits.reserve(runs);
    bestNoveltyFits.reserve(runs);
    for (pair<int, int> vals: bestFits) {
        bestMatchFits.push_back(vals.first);
        bestNoveltyFits.push_back(vals.second);
    }

    vector<double> stats = calcStats<int>(bestMatchFits, biggerBetter); // {mean, stdDev, 95CI, bestVal, bestIdx}
    int bestMatchFit = (int) stats[3];
    vector<int> noveltyFitsSameMatchFit;
    vector<int> noveltyIdxs;
    noveltyFitsSameMatchFit.reserve(runs);
    noveltyIdxs.reserve(runs);

    for (int idx = 0; idx < runs; idx++) {
        if (bestMatchFits[idx] == bestMatchFit) {
            noveltyFitsSameMatchFit.push_back(bestNoveltyFits[idx]);
            noveltyIdxs.push_back(idx);
        }
    }

    vector<double> stats2 = calcStats<int>(noveltyFitsSameMatchFit, false); // {mean, stdDev, 95CI, bestVal, bestIdx}
    int bestRun = noveltyIdxs[(int) stats2[4]] + 1;

    multiStream printAndSave(cout, outp);
    printAndSave << "Experiment Report:" << "\n";
    printAndSave << "Best Run: " << bestRun << "\n";
    printAndSave << "Best Fitness: " << stats[3] << " of ";
    printAndSave << seqLen << " characters matched and a novelty score of ";
    printAndSave << stats2[3] << "\n";
    printAndSave << "Match Fitness 95% CI: " << stats[0] << " +- " << stats[2] << "\n";
    printAndSave << "Novelty Fitness 95% CI (for best match fitness): ";
    printAndSave << stats2[0] << " +- " << stats2[2] << "\n";
//    printAndSave << "\n";
    printAndSave << left << setw(20) << "Best Match: ";
    bestSDA.fillOutput(testSeq);
    intToChar(testSeq, charSeq);
    printVector<multiStream, char>(printAndSave, charSeq, "", "", true);
    printAndSave << left << setw(20) << "Matches: ";
    printMatches<multiStream, int>(printAndSave, testSeq, goalSeq, true);
    printAndSave << left << setw(20) << "Desired Sequence: ";
    intToChar(goalSeq, charSeq);
    printVector<multiStream, char>(printAndSave, charSeq, "", "", true);
//    printAndSave << "\nSDA:\n";
    bestSDA.printSDA(cout);
    bestSDA.printSDA(outp);
    return 0;
}

int sdaCheck(ofstream &outStream, int currentGen) {
    outStream << "Population After " << currentGen << " Mating Events" << endl;

    for (int sda = 0; sda < popsize; sda++) {
        outStream << "SDA " << sda << endl;
        outStream << "Fitness: " << matchFits[sda] << endl;
        pop[sda].printSDA(outStream);
        pop[sda].fillOutput(testSeq, true, outStream);
    }
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

int printPopFits(ostream &outStrm) {
    outStrm << "Fitness Values [Matching, Novelty]: ";
    vector<int> sortedIdxs = tournSelect(popsize, BIGGER_BETTER);
    bool first = true;
    for (int idx: sortedIdxs) {
        if (!first) {
            outStrm << ", ";
        }
        outStrm << "[" << matchFits[idx] << ", " << noveltyFits[idx] << "]";
        first = false;
    }
    outStrm << "\n";
    return 0;
}