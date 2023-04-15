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

class Bitsprayer {
public:
    Bitsprayer();           //creates an unallocated bitspray
    explicit Bitsprayer(int states, int numChars);      //create a bitspray with buffer S states
    Bitsprayer(Bitsprayer &other);  //copy constructor
    ~Bitsprayer();                //destructor

    int create(int states);
    int randomize();
    int copy(Bitsprayer &other);
    int print();
    int print(ostream &aus);
    static int destroy();
    int twoPtCrossover(Bitsprayer &other);
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
