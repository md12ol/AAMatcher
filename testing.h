#include "main.h"

#define mateTests (int) 100

vector<int> runMultiCross(int numEvents, SDA mom, SDA dad) {
    SDA child1, child2;
    vector<int> fitVals;

    fitVals.reserve(numEvents * 2);

    for (int i = 0; i < numEvents; ++i) {
        child1.copy(mom);
        child2.copy(dad);
        if (crossoverOp == 0) child1.twoPtCrossover(child2);
        else if (crossoverOp == 1) child1.oneStateCrossover(child2);
        fitVals.push_back((int) fitness(child1));
        fitVals.push_back((int) fitness(child2));
    }
    return fitVals;
}

vector<int> runMultiMutate(int numEvents, SDA parent) {
    SDA child;
    vector<int> fitVals;
    int numMuts;

    fitVals.reserve(numEvents);
    for (int i = 0; i < numEvents; ++i) {
        child.copy(parent);
        numMuts = (int) lrand48() % maxMuts + 1;
        child.mutate(numMuts);
        fitVals.push_back((int) fitness(child));
    }
    return fitVals;
}

int crossoverCheck(ofstream &outStream) {
    outStream << "Population Fitness Values:" << endl;

    for (int i = 0; i < popsize; ++i) {
        outStream << fits[i] << endl;
    }

    outStream << "Population Crossover Checks:" << endl;

    for (int mom = 0; mom < popsize; mom++) {
        for (int dad = mom + 1; dad < popsize; dad++) {
            outStream << "Parent Idxs: " << mom << "\t" << dad << endl;
            outStream << "Parent Fits: " << fits[mom] << "\t" << fits[dad] << endl;
            printVector<ofstream, int>(outStream, runMultiCross(mateTests, pop[mom], pop[dad]), "", "\n", true);
        }
    }
    return 0;
}

int mutateCheck(ofstream &outStream) {
    outStream << "Population Fitness Values:" << endl;

    for (int i = 0; i < popsize; ++i) {
        outStream << fits[i] << endl;
    }

    outStream << "Population Mutation Checks:" << endl;

    for (int mom = 0; mom < popsize; mom++) {
        outStream << "Parent Idx: " << mom << endl;
        outStream << "Parent Fit: " << fits[mom] << endl;
        printVector<ofstream, int>(outStream, runMultiMutate(mateTests, pop[mom]), "", "\n", true);
    }
    return 0;
}

int sdaCheck(ofstream &outStream) {
    outStream << "Population SDA Check:" << endl;

    for (int sda = 0; sda < popsize; sda++) {
        outStream << "SDA " << sda << endl;
        outStream << fits[sda] << endl;
        pop[sda].printSDA(outStream);
        pop[sda].fillOutput(testSeq, true, outStream);
    }
    return 0;
}