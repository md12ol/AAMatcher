#include <iostream>
#include <vector>

using namespace std;

class SDA {
public:
    SDA();
    explicit SDA(int numStates, int numChars, int maxRespLen, int outputLen, int initState = 0, bool verbose = false);
    SDA(SDA &other);
    ~SDA();

    int setOutputLen(int newLen);
    int randomize();
    int copy(SDA &other);
    int twoPtCrossover(SDA &other);
    int oneStateCrossover(SDA &other);
    int mutate(int numMuts);
    vector<int> rtnOutput(bool printToo = false, ostream &outStrm = cout);
    int fillOutput(vector<int> &output, bool printToo = false, ostream &outStrm = cout);
    int printSDA(ostream &outStrm);

private:
    int create();

    int initChar{};
    int numStates{};
    int initState{};
    int curState{};
    int numChars{};
    int maxRespLen{};
    int outputLen{};
    bool verbose{};

    vector<vector<int> > transitions;
    vector<vector<vector<int> > > responses;
};
