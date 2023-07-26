#include "main.h"

#define mateTests (int) 100

int matingEvent(bool biggerBetter, int currentGen, ostream &outp) {
    int numMuts;
    SDA child1a, child2a, child1b, child2b;
    double fit1a, fit2a, fit1b, fit2b;

    vector<int> idxs = tournSelect(tournSize, biggerBetter);

    child1a.copy(pop[idxs[0]]);
    fit1a = fits[idxs[0]];
    child2a.copy(pop[idxs[1]]);
    fit2a = fits[idxs[1]];
    if (drand48() < crossoverRate) {
        if (crossoverOp == 0) child1a.twoPointCrossover(child2a);
        else if (crossoverOp == 1) child1a.oneStateCrossover(child2a);

        fit1a = fitness(child1a);
        fit2a = fitness(child2a);
        // If there is an improvement...
        if (((fit1a > populationBestFit && BIGGER_BETTER) || (fit2a > populationBestFit && BIGGER_BETTER)) ||
            ((fit1a < populationBestFit && !BIGGER_BETTER) || (fit2a < populationBestFit && !BIGGER_BETTER))) {
            outp << "Crossover Improvement during Mating Event " << currentGen << endl;
            if (BIGGER_BETTER){
                outp << populationBestFit << " -> " << max(fit1a, fit2a) << endl;
            } else {
                outp << populationBestFit << " -> " << min(fit1a, fit2a) << endl;
            }
            outp << "Parent 1: " << fits[idxs[0]] << endl;
            pop[idxs[0]].printSDA(outp);
            outp << "Parent 2: " << fits[idxs[1]] << endl;
            pop[idxs[1]].printSDA(outp);
            outp << "Child 1: " << fit1a << endl;
            child1a.printSDA(outp);
            outp << "Child 2: " << fit2a << endl;
            child2a.printSDA(outp);

            if (BIGGER_BETTER) {
                populationBestFit = max(fit1a, fit2a);
            } else {
                populationBestFit = min(fit1a, fit2a);
            }
        }
    }

    child1b.copy(child1a);
    fit1b = fit1a;
    child2b.copy(child2a);
    fit2b = fit2a;

    if (drand48() < mutationRate && maxMuts > 0) {
        numMuts = (int) lrand48() % maxMuts + 1;
        child1b.mutate(numMuts);
        numMuts = (int) lrand48() % maxMuts + 1;
        child2b.mutate(numMuts);

        fit1b = fitness(child1b);
        fit2b = fitness(child2b);
        // If there is an improvement...
        if (((fit1b > populationBestFit && BIGGER_BETTER) || (fit2b > populationBestFit && BIGGER_BETTER)) ||
            ((fit1b < populationBestFit && !BIGGER_BETTER) || (fit2b < populationBestFit && !BIGGER_BETTER))) {
            outp << "Mutation Improvement during Mating Event " << currentGen << endl;
            if (BIGGER_BETTER){
                outp << populationBestFit << " -> " << max(fit1b, fit2b) << endl;
            } else {
                outp << populationBestFit << " -> " << min(fit1b, fit2b) << endl;
            }
            outp << "Child 1 Before: " << fit1a << endl;
            child1a.printSDA(outp);
            outp << "Child 1 After: " << fit1b << endl;
            child1b.printSDA(outp);
            outp << "Child 2 Before: " << fit2a << endl;
            child2a.printSDA(outp);
            outp << "Child 2 After: " << fit2b << endl;
            child2b.printSDA(outp);

            if (BIGGER_BETTER) {
                populationBestFit = max(fit1b, fit2b);
            } else {
                populationBestFit = min(fit1b, fit2b);
            }
        }
    }

    pop[idxs.end()[-1]] = child1b;
    pop[idxs.end()[-2]] = child2b;

    fits[idxs.end()[-1]] = fit1b;
    fits[idxs.end()[-2]] = fit2b;
    return 0;
}

vector<int> runMultiCross(int numEvents, SDA mom, SDA dad) {
    SDA child1, child2;
    vector<int> fitVals;

    fitVals.reserve(numEvents * 2);

    for (int i = 0; i < numEvents; ++i) {
        child1.copy(mom);
        child2.copy(dad);
        if (crossoverOp == 0) child1.twoPointCrossover(child2);
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

int sdaCheck(ofstream &outStream, int currentGen) {
    outStream << "Population After " << currentGen << " Mating Events"<< endl;

    for (int sda = 0; sda < popsize; sda++) {
        outStream << "SDA " << sda << endl;
        outStream << "Fitness: " << fits[sda] << endl;
        pop[sda].printSDA(outStream);
        pop[sda].fillOutput(testSeq, true, outStream);
    }
    return 0;
}