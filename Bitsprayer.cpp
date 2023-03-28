#include "Bitsprayer.h"

#define VERBOSEBIT false

Bitsprayer::Bitsprayer() {
    initInput = -1;
    numStates = -1;
    initState = -1;
    curState = -1;
//    transitions = {{}};
//    responses = {{{}}};
//    zeroProb = -1;
    if (VERBOSEBIT) cout << "Bitsprayer Made" << endl;
}

Bitsprayer::Bitsprayer(int states, int numChars) {
    initInput = -1;
    numStates = states;
    initState = -1;
    curState = -1;
    this->numChars = numChars;
    transitions.reserve(states);
//    zeroProb = prob;

    for (vector<int> v: transitions) {
        v.reserve(numChars);
    }
    responses.reserve(states);
    for (vector<vector<int> > v: responses) {
        v.reserve(numChars);
    }

//    transitions = {{}};
//    responses = {{{}}};
    if (VERBOSEBIT) cout << "Bitsprayer Made w " << states << " states" << endl;
}

Bitsprayer::Bitsprayer(Bitsprayer &other) {
    copy(other);
}

Bitsprayer::~Bitsprayer() {
    destroy();
}

int Bitsprayer::create(int states) {
    numStates = states;
    randomize();
    return 0;
}

int Bitsprayer::randomize() {
    initInput = (int) lrand48() % numChars;
    if (initInput < 0) {
        cout << "ERROR! initInput is negative!" << endl;
    }
    initState = 0;
    curState = -1;

    transitions.clear();
    responses.clear();

    vector<int> oneState;
    for (int s = 0; s < numStates; ++s) {
        oneState.clear();
        for (int i = 0; i < numChars; ++i) {
            oneState.push_back((int) lrand48() % numStates);
        }
        transitions.push_back(oneState);
    }

    vector<int> oneResponse;
    vector<vector<int>> oneStateResps;
    int respSize;
    for (int s = 0; s < numStates; ++s) {
        oneStateResps.clear();
        for (int t = 0; t < numChars; ++t) {
            oneResponse.clear();
            respSize = (int) lrand48() % 2 + 1;
            for (int i = 0; i < respSize; ++i) {
//                if (drand48() < zeroProb) {
//                    oneResponse.push_back(0);
//                } else {
                    oneResponse.push_back((int) lrand48() % (numChars - 1) + 1);
//                }
                //TODO: Check!
            }
            oneStateResps.push_back(oneResponse);
        }
        responses.push_back(oneStateResps);
    }
    return 0;
}

int Bitsprayer::copy(Bitsprayer &other) {
    initInput = other.initInput;
    numStates = other.numStates;
    initState = other.initState;
    curState = other.curState;
    this->numChars = other.numChars;
//    zeroProb = other.zeroProb;

    transitions = other.transitions;
    responses = other.responses;
    if (VERBOSEBIT) cout << "Bitsprayer Copied" << endl;
    return 0;
}

int Bitsprayer::print() {
    print(cout);
    return 0;
}

int Bitsprayer::print(ostream &aus) {
    aus << initState << " <- " << initInput << endl;
    for (int s = 0; s < numStates; ++s) {
        if (transitions[s].size() > numChars) {
            aus << "ERROR!  More transitions than MAXVAL!" << endl;
        }
        if (responses[s].size() > numChars) {
            aus << "ERROR!  More responses than MAXVAL!" << endl;
        }
        for (int t = 0; t < numChars; ++t) {
            aus << s << " + " << t << " -> " << transitions.at(s).at(t) << " [";
            if (responses[s][t].size() > 2) {
                aus << "ERROR!  Response length more than 2!" << endl;
            }
            for (int v: responses.at(s).at(t)) {
                aus << " " << v;
            }
            aus << " ]" << endl;
        }
    }
    if (transitions.size() > numStates) {
        aus << "ERROR!  More transitions than the number of states!" << endl;
    }
    if (responses.size() > numStates) {
        aus << "ERROR!  More responses than the number of states!" << endl;
    }
    if (VERBOSEBIT) cout << "Bitsprayer Printed" << endl;
    return 0;
}

int Bitsprayer::destroy() {
    return 0;
}

int Bitsprayer::twoPtCrossover(Bitsprayer &other) {
    int cp1, cp2;
    int swapInt;
    vector<int> swapVec;

    if (numStates != other.numStates) {
        return 1;
    }

    do {
        cp1 = (int) lrand48() % numStates;
        cp2 = (int) lrand48() % numStates;
        if (cp1 > cp2) {
            swapInt = cp1;
            cp1 = cp2;
            cp2 = swapInt;
        }
    } while (cp1 == cp2);

    if (cp1 == 0) {
        swapInt = initInput;
        initInput = other.initInput;
        other.initInput = swapInt;
    }

    for (int s = cp1; s < cp2; s++) {
        swapVec = transitions.at(s);
        transitions.at(s) = other.transitions.at(s);
        other.transitions.at(s) = swapVec;
        swapVec = responses.at(s).at(0);
        responses.at(s).at(0) = other.responses.at(s).at(0);
        other.responses.at(s).at(0) = swapVec;
        swapVec = responses.at(s).at(1);
        responses.at(s).at(1) = other.responses.at(s).at(1);
        other.responses.at(s).at(1) = swapVec;
    }
    return 0;
}

int Bitsprayer::mutate(int numMuts) {
    int mutPt;
    vector<int> oneResponse;
    int respSize;

    for (int m = 0; m < numMuts; ++m) {
        mutPt = (int) lrand48() % (2 * numStates + 1);

        if (mutPt == 0) {
            initInput = (int) lrand48() % numChars;
            return 0;
        }
        mutPt = (mutPt - 1) / 2;
        int transNum = (int) lrand48() % numChars;
        if ((int) lrand48() % 2 == 0) { // Mutate transition
            transitions.at(mutPt).at(transNum) = (int) lrand48() % numStates;
        } else { // Mutate response
            oneResponse.clear();
            respSize = (int) lrand48() % 2 + 1;
            for (int i = 0; i < respSize; ++i) {
                oneResponse.push_back((int) lrand48() % numChars);
            }
            responses.at(mutPt).at(transNum) = oneResponse;
        }
    }
    return 0;
}

vector<int> Bitsprayer::getBitsVec(int len) {
    vector<int> rtn;
    rtn.clear();
    int nextBit;
    int cnt = 0;
    curState = initState;
    buf.clear();
    buf.push_back(initInput);

    while (cnt < len) {
        nextBit = buf.front();
        rtn.push_back(nextBit);
        buf.erase(buf.begin());
        for (int i: responses.at(curState).at(nextBit)) {
            buf.push_back(i);
        }
        curState = transitions.at(curState).at(nextBit);
        cnt++;
    }
    return rtn;
}

int Bitsprayer::printBitsVec(int len, ostream &aus) {
    vector<int> vec = getBitsVec(len);
    for (int i: vec) {
        aus << i << " ";
    }
    aus << endl;
    return 0;
}

//int main() {
//    Bitsprayer b(5);
//    b.create(5);
//    b.randomize();
//    b.print();
//
////    Bitsprayer b2(5);
////    b2.create(5);
////    b2.randomize();
////    b2.print();
////
////    b.twoPtCrossover(b2);
////    b.print();
////    b2.print();
//
//    printVec(b.getBitsVec(20));
//
//    if (VERBOSEBIT) cout << "DONE!" << endl;
//    return 0;
//}