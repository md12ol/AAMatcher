#ifndef BITSPRAYER_H
#define BITSPRAYER_H

#include <cmath>
#include <list>
#include <sstream>
#include <vector>
#include <iostream>
#include <ctime>
#include <fstream>
#include <cstdlib>

using namespace std;

class SDA {
public:
    SDA();           //creates an unallocated bitspray
    explicit SDA(int states, int numChars);      //create a bitspray with buffer S states
    SDA(SDA &other);  //copy constructor
    ~SDA();                //destructor

    int create(int states);
    int randomize();
    int copy(SDA &other);
    int print();
    int print(ostream &aus);
    static int destroy();
    int twoPtCrossover(SDA &other);
    int oneStateCrossover(SDA &other);
    int mutate(int numMuts);
    vector<int> getBitsVec(int len);
    int printBitsVec(int len, ostream &aus);

private:
    int initInput;
    int numStates;
    int initState;
    int curState;
    int numChars;
//    double zeroProb;
    vector<int> buf;
    vector<vector<int> > transitions;
    vector<vector<vector<int> > > responses;
};

#endif //MEDICALPREDICTOR_BITSPRAYER_H
