#include <iomanip>
#include <algorithm>
#include <vector>
#include "Bitsprayer.h"

// This block of code includes the filesystem package (which depends on OS)
#if defined(__cplusplus) && __cplusplus >= 201703L && defined(__has_include)
#if __has_include(<filesystem>)
#define GHC_USE_STD_FS

#include <filesystem>
#include <thread>

namespace filesystem = std::filesystem;
#endif
#endif
#ifndef GHC_USE_STD_FS
#include "filesystem.hpp"
namespace filesystem = ghc::filesystem;
#endif

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

vector<int> testSeq;
int seqLen;

Bitsprayer *pop;
vector<double> fits;

string pathToSeqs = "./Sequences.dat";
char pathToOut[150];

int getArgs(char *args[]);
int initAlg();
int initPop(int run);
double fitness(Bitsprayer &sda);
vector<double> calcStats(bool biggerBetter);
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
double finalReport(ostream &outp, bool biggerBetter);
double report(ofstream &outp, int run, int rptNum, bool biggerBetter);
int culling(double percentage, bool rndPick, bool biggerBetter);
vector<vector<int>> getSequences();
vector<int> seqToVector(const string& seq);
int makeReadMe(ostream &outp);
vector<double> calcStats2(vector<int> vals,bool biggerBetter);

/**
 * Input variables:
 * 1. Population size
 * 2. Number of characters (in SDA)
 * 3. Number of states
 * 4. Seed
 * 5. Number of runs
 * 6. Number of mating events
 * 7. Maximum number of mutations
 * 8. Sequence number
 * 9. Tournament size
 */

using namespace std;

int initAlg() {
    srand48(time(nullptr)); // read the random number seed
//    generateTestSeq();
    testSeq = getSequences()[seqNum];
    seqLen = (int)testSeq.size();
    fits.reserve(popsize);

    pop = new Bitsprayer[popsize];
    for (int idx = 0; idx < popsize; ++idx) {
        pop[idx] = Bitsprayer(sdaStates, numChars);
    }
    return 0;
}

int initPop(int run) {
    cout << "Beginning Run " << run << " of " << runs << endl;
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
    printVector<ostream, int>(outp, testSeq, "", "", true);
    outp << "See readme.txt for system and experiment parameters." << endl;
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
    char letter;
    string thing = "ATGGGACGCAAGGACGAGCAGAAGCAAACGAGCGCCACAAGCACGCCGGGGCAGGGG";
    for (int i = 0; i < thing.size(); i++) {
        letter = thing[i];
        if (letter == 'G') {
            seq.push_back(0);
        } else if (letter == 'C') {
            seq.push_back(1);
        } else if (letter == 'A') {
            seq.push_back(2);
        } else if (letter == 'T') {
            seq.push_back(3);
        }
    }
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

    cout << "Arguments Captured!" << endl;
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