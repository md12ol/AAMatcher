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
#define REPORTS 100
#define CULLINGS 25

// Experiment parameters
int popsize;
int numChars;
int sdaStates;
double *sdaProbs;
int seed;
int runs;
int generations;
int maxMuts;
int tournSize = 7;

vector<int> testSeq;
int seqLen = 50;

Bitsprayer *pop;
vector<double> fits;

string filename;




int getArgs(char *args[]);
int initAlg();
int initPop(int run);
double fitness(Bitsprayer &sda);
vector<double> calcStats(bool biggerBetter);
int printExpStatsHeader(ostream &outp);
int matingEvent(bool biggerBetter);
vector<int> tournSelect(int size ,bool decreasing);
int cmdLineIntro(ostream &outp);
int generateTestSeq();
bool compareFitness(int popIdx1, int popIdx2);
template<class T1, class T2>
int printVector(T1 &outp, vector<T2> vec, const string &msg, const string &sep, bool newline);
template<class T1, class T2>
int printIdxsOfVector(T1 &outp, vector<T2> vec, const vector<int> &idxs, const string &msg, const string &sep, bool newline);
double finalReport(ostream &outp, bool biggerBetter);
int report(ofstream &outp, int run, int rptNum, bool biggerBetter);
int culling(double percentage);

using namespace std;

int initAlg() {
    srand48(time(nullptr)); // read the random number seed
    generateTestSeq();
    fits.reserve(popsize);

    pop = new Bitsprayer[popsize];
    for (int idx = 0; idx < popsize; ++idx) {
        pop[idx] = *new Bitsprayer(sdaStates, numChars);
    }
    return 0;
}

int cmdLineIntro(ostream &outp) {
    outp << "Amino Acid Sequence Matcher." << endl;
    outp << "Fitness Function: Naive." << endl;
    outp << "Marching Sequence (length " << seqLen << "): ";
    printVector<ostream, int>(cout, testSeq, "", "", true);
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
int printIdxsOfVector(T1 &outp, vector<T2> vec, const vector<int> &idxs, const string &msg, const string &sep, bool newline) {
    outp << msg;
    for (auto idx: idxs) {
        outp << vec[idx] << sep;
    }
    if (newline) outp << "\n";
    return 0;
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

    arg = args[6]; // generations
    try {
        size_t pos;
        generations = std::stoi(arg, &pos);
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