#include <iomanip>
#include <algorithm>
#include <vector>
#include <filesystem>
#include <thread>
#include "SDA.h"
namespace filesystem = std::filesystem;

// System parameters
#define REPORT_EVERY 10000
#define TERM_CRIT 50
#define CULLING_EVERY 4
#define CULLING_ODDS 0.25
#define RANDOM_CULLING false
#define BIGGERBETTER true

// Experiment parameters
int popsize;
int numChars;
int sdaStates;
double *sdaProbs;
int seed;
int runs;
int maxGens;
int maxMuts;
int tournSize;
int seqNum;
int crossOp;

vector<int> goalSeq;
int seqLen;

SDA *pop;
vector<double> fits;

string pathToSeqs = "./Sequences.dat";
char pathToOut[150];

int getArgs(char *args[]);
int initAlg();
int initPop(int run);
double fitness(SDA &sda);
int printExpStatsHeader(ostream &outp);
int matingEvent(bool biggerBetter);
vector<int> tournSelect(int size, bool decreasing);
int cmdLineIntro(ostream &outp);
vector<int> generateTestSeq();
bool compareFitness(int popIdx1, int popIdx2);
template<class T1, class T2>
int printVector(T1 &outp, vector<T2> vec, const string &msg, const string &sep, bool newline);
template<class T1, class T2>
int printIdxsOfVector(T1 &outp, vector<T2> vec, const vector<int> &idxs, const string &msg, const string &sep,
                      bool newline);
template<class T1, class T2>
int printMatches(T1 &outp, const vector<T2> &test, const vector<T2> &goal, bool newline);
int runReport(ostream &outp, bool biggerBetter);
double report(ofstream &outp, int run, int rptNum, bool biggerBetter);
int culling(double percentage, bool rndPick, bool biggerBetter);
vector<vector<int>> getSequences();
vector<int> seqToVector(const string& seq);
int makeReadMe(ostream &outp);
template<class T>
vector<double> calcStats(vector<T> vals, bool biggerBetter);
vector<char> intToChar(const vector<int>& seq);
int expReport(ostream &outp, vector<double> bestFits, SDA bestSDA, bool biggerBetter);

/**
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
 * 11.
 */

using namespace std;

double fitness(SDA &sda) {
    double val = 0.0;
    vector<int> result = sda.getBitsVec(seqLen);
    for (int i = 0; i < seqLen; ++i) {
        if (result[i] == goalSeq[i]) {
            val += 1;
        }
    }
    return val;
}

template<class T>
vector<double> calcStats(vector<T> vals, bool biggerBetter) {
    double sum = 0.0;
    double bestVal = (biggerBetter ? 0.0 : MAXFLOAT);

    for (T val: vals) {
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

int initAlg() {
    srand48(time(nullptr)); // use system time as random number seed
    srand48(seed);           // read the random number seed
    goalSeq = getSequences()[seqNum];
    seqLen = (int)goalSeq.size();
    fits.reserve(popsize);

    pop = new SDA[popsize];
    for (int idx = 0; idx < popsize; ++idx) {
        pop[idx] = SDA(sdaStates, numChars);
    }
    return 0;
}

int initPop(int run) {
    if (runs > 1){
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

int cmdLineIntro(ostream &outp) {
    outp << "Amino Acid Sequence Matcher." << endl;
    outp << "Fitness Function: Naive." << endl;
    outp << "Matching Sequence (length " << seqLen << "): ";
    printVector<ostream, char>(outp, intToChar(goalSeq), "", "", true);
    outp << "See read.me for system and experiment parameters." << endl;
    outp << endl;
    return 0;
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

int makeReadMe(ostream &outp){
    cmdLineIntro(outp);
    outp << "Experiment Parameters" << endl;
    outp << "Population Size: " << popsize << endl;
    outp << "Number of States: " << sdaStates << endl;
    outp << "Alphabet Size: " << numChars << endl;
    outp << "Max Generations: " << maxGens << endl;
    outp << "Report Every: " << REPORT_EVERY << " generations" << endl;
    outp << "Termination Criteria: maximum fitness or no improvement for " << TERM_CRIT * REPORT_EVERY << " maxGens" <<endl;
    outp << "Maximum Number of Mutations: " << maxMuts << endl;
    outp << "Tournament Size: " << tournSize << endl;
    outp << "Culling Every: " << CULLING_EVERY * REPORT_EVERY << " generations" << endl;
    if (RANDOM_CULLING){
        outp << "Culling Winners: random " << (int) (CULLING_ODDS * 100) << "% of the population" << endl;
    } else {
        outp << "Culling Winners: worst " << (int) (CULLING_ODDS * 100) << "% of the population" << endl;
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

vector<int> seqToVector(const string& seq) {
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

vector<char> intToChar(const vector<int>& seq) {
    vector<char> sequence;
    for (int i: seq) {
        if (i == 0) {
            sequence.push_back('G');
        } else if (i == 1) {
            sequence.push_back('C');
        } else if (i == 2) {
            sequence.push_back('A');
        } else if (i == 3) {
            sequence.push_back('T');
        }
    }
    return sequence;
}

template<class T1, class T2>
int printMatches(T1 &outp, const vector<T2> &test, const vector<T2> &goal, bool newline){
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

vector<vector<int>> getSequences() {
    string tmp;
    ifstream in(pathToSeqs);
    vector<vector<int>> rtn;
    for (int seqIdx = 0; seqIdx < 5; ++seqIdx) {
        if (in.is_open()) {
            getline(in, tmp);
            getline(in, tmp);
            getline(in, tmp);
            rtn.push_back(seqToVector(tmp));
            getline(in, tmp);
        }
    }
    in.close();
    rtn.push_back(generateTestSeq());
    return rtn;
}

vector<int> generateTestSeq() {
    vector<int> seq;
    string thing = "ATGGGACGCAAGGACGAGCAGAAGCAAACGAGCGCCACAAGCACGCCGGGGCAGGGG";
    seq = seqToVector(thing);
    return seq;
}

int getArgs(char *args[]) {
    string arg = args[1]; // popsize
    try {
        size_t pos;
        popsize = std::stoi(arg, &pos);
        if (pos < arg.size()) {
            std::cerr << "Trailing characters after number: " << arg << endl;
        }
    } catch (std::invalid_argument const &ex) {
        std::cerr << "Invalid number: " << arg << '\n';
    } catch (std::out_of_range const &ex) {
        std::cerr << "Number out of range: " << arg << '\n';
    }

    arg = args[2]; // numChars
    try {
        size_t pos;
        numChars = std::stoi(arg, &pos);
        if (pos < arg.size()) {
            std::cerr << "Trailing characters after number: " << arg << endl;
        }
    } catch (std::invalid_argument const &ex) {
        std::cerr << "Invalid number: " << arg << '\n';
    } catch (std::out_of_range const &ex) {
        std::cerr << "Number out of range: " << arg << '\n';
    }

    arg = args[3]; // sdaStates
    try {
        size_t pos;
        sdaStates = std::stoi(arg, &pos);
        if (pos < arg.size()) {
            std::cerr << "Trailing characters after number: " << arg << endl;
        }
    } catch (std::invalid_argument const &ex) {
        std::cerr << "Invalid number: " << arg << '\n';
    } catch (std::out_of_range const &ex) {
        std::cerr << "Number out of range: " << arg << '\n';
    }

    arg = args[4]; // seed
    try {
        size_t pos;
        seed = std::stoi(arg, &pos);
        if (pos < arg.size()) {
            std::cerr << "Trailing characters after number: " << arg << endl;
        }
    } catch (std::invalid_argument const &ex) {
        std::cerr << "Invalid number: " << arg << '\n';
    } catch (std::out_of_range const &ex) {
        std::cerr << "Number out of range: " << arg << '\n';
    }

    arg = args[5]; // runs
    try {
        size_t pos;
        runs = std::stoi(arg, &pos);
        if (pos < arg.size()) {
            std::cerr << "Trailing characters after number: " << arg << endl;
        }
    } catch (std::invalid_argument const &ex) {
        std::cerr << "Invalid number: " << arg << '\n';
    } catch (std::out_of_range const &ex) {
        std::cerr << "Number out of range: " << arg << '\n';
    }

    arg = args[6]; // maxGens
    try {
        size_t pos;
        maxGens = std::stoi(arg, &pos);
        if (pos < arg.size()) {
            std::cerr << "Trailing characters after number: " << arg << endl;
        }
    } catch (std::invalid_argument const &ex) {
        std::cerr << "Invalid number: " << arg << '\n';
    } catch (std::out_of_range const &ex) {
        std::cerr << "Number out of range: " << arg << '\n';
    }

    arg = args[7]; // maxMuts
    try {
        size_t pos;
        maxMuts = std::stoi(arg, &pos);
        if (pos < arg.size()) {
            std::cerr << "Trailing characters after number: " << arg << endl;
        }
    } catch (std::invalid_argument const &ex) {
        std::cerr << "Invalid number: " << arg << '\n';
    } catch (std::out_of_range const &ex) {
        std::cerr << "Number out of range: " << arg << '\n';
    }

    arg = args[8]; // seqNum
    try {
        size_t pos;
        seqNum = std::stoi(arg, &pos);
        if (pos < arg.size()) {
            std::cerr << "Trailing characters after number: " << arg << endl;
        }
    } catch (std::invalid_argument const &ex) {
        std::cerr << "Invalid number: " << arg << '\n';
    } catch (std::out_of_range const &ex) {
        std::cerr << "Number out of range: " << arg << '\n';
    }

    arg = args[9]; // tournSize
    try {
        size_t pos;
        tournSize = std::stoi(arg, &pos);
        if (pos < arg.size()) {
            std::cerr << "Trailing characters after number: " << arg << endl;
        }
    } catch (std::invalid_argument const &ex) {
        std::cerr << "Invalid number: " << arg << '\n';
    } catch (std::out_of_range const &ex) {
        std::cerr << "Number out of range: " << arg << '\n';
    }

    arg = args[10]; // crossOp
    try {
        size_t pos;
        crossOp = std::stoi(arg, &pos);
        if (pos < arg.size()) {
            std::cerr << "Trailing characters after number: " << arg << endl;
        }
    } catch (std::invalid_argument const &ex) {
        std::cerr << "Invalid number: " << arg << '\n';
    } catch (std::out_of_range const &ex) {
        std::cerr << "Number out of range: " << arg << '\n';
    }

    cout << "Arguments Captured!" << endl;
    return 0;
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

int matingEvent(bool biggerBetter) {
    int numMuts;
    SDA child1, child2;

    vector<int> idxs = tournSelect(tournSize, biggerBetter);
//    printIdxsOfVector<double>(fits, idxs, "Tournament Fitness Values: ", "\t", true);

    child1.copy(pop[idxs[0]]);
    child2.copy(pop[idxs[1]]);
    if (crossOp == 0) child1.twoPtCrossover(child2);
    else if (crossOp == 1) child1.oneStateCrossover(child2);

    if (maxMuts > 0){
        numMuts = (int) lrand48() % maxMuts + 1;
        child1.mutate(numMuts);
        numMuts = (int) lrand48() % maxMuts + 1;
        child2.mutate(numMuts);
    }

    pop[idxs.end()[-1]] = child1;
    pop[idxs.end()[-2]] = child2;

    fits[idxs.end()[-1]] = fitness(child1);
    fits[idxs.end()[-2]] = fitness(child2);
    return 0;
}

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

double report(ofstream &outp, int run, int rptNum, bool biggerBetter) {
    vector<double> stats = calcStats<double>(fits,biggerBetter); // {mean, stdDev, 95CI, best}
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
    vector<int> bestSeq = pop[bestIdx].getBitsVec(seqLen);
    printAndSave << left << setw(20) << "Best Match: ";
    printVector<multiStream, char>(printAndSave, intToChar(bestSeq), "", "", true);
    printAndSave << left << setw(20) << "Matches: ";
    printMatches<multiStream, char>(printAndSave, intToChar(bestSeq), intToChar(goalSeq), true);
    printAndSave << left << setw(20) << "Desired Sequence: ";
    printVector<multiStream, char>(printAndSave, intToChar(goalSeq), "", "", true);

    outp << "SDA" << endl;
    pop[bestIdx].print(outp);

    vector<double> fitsCopy = fits;
    sort(fitsCopy.begin(), fitsCopy.end());
    if (biggerBetter) {
        reverse(fitsCopy.begin(), fitsCopy.end());
    }
    printVector<multiStream, double>(printAndSave, fitsCopy, "Fitness Values: ", " ", true);
    cout << endl;
    return bestIdx;
}

int expReport(ostream &outp, vector<double> bestFits, SDA bestSDA, bool biggerBetter){
    auto maxIterator = minmax_element(bestFits.begin(), bestFits.end());
    int bestIdx;
    if (biggerBetter) {
        bestIdx = (int) distance(bestFits.begin(), maxIterator.second);
    } else {
        bestIdx = (int) distance(bestFits.begin(), maxIterator.first);
    }

    vector<double> stats = calcStats<double>(bestFits, BIGGERBETTER);
    multiStream printAndSave(cout, outp);
    printAndSave << "Experiment Report:" << "\n";
    printAndSave << "Best Fitness: " << bestFits[bestIdx] << " of " << seqLen << "\n";
    printAndSave << "Fitness 95% CI: " << stats[0] << " +- " << stats[2] << "\n";
    printAndSave << "\n";
    printAndSave << left << setw(20) << "Best Match: ";
    printVector<multiStream, char>(printAndSave, intToChar(bestSDA.getBitsVec(seqLen)), "", "", true);
    printAndSave << left << setw(20) << "Matches: ";
    printMatches<multiStream, char>(printAndSave, intToChar(bestSDA.getBitsVec(seqLen)), intToChar(goalSeq), true);
    printAndSave << left << setw(20) << "Desired Sequence: ";
    printVector<multiStream, char>(printAndSave, intToChar(goalSeq), "", "", true);
    printAndSave << "\nSDA:";
    bestSDA.print(cout);
    bestSDA.print(outp);
    return 0;
}

