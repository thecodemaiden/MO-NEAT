//
//  ExampleGoals.cpp
//  NEATBox
//
//  Created by Adeola Bannis on 1/7/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#include "ExampleGoals.h"
#include "ExampleNetwork.h"
#include <array>
#include <cmath>

double dummyEvaluation(ExampleNetwork& individual)
{
    std::pair<std::vector<double>, bool> eval =
    individual.simulateTillEquilibrium(std::vector<double>(), 10);
    
    for (long i=0; i<eval.first.size(); i++) {
        std::cout <<"\t" << i << ": " << eval.first[i] <<"\n";
    }
    std::cout << "\tsteady: " << (eval.second ? "y" : "n") << "\n";
    
    return individual.numberOfEdges();
}

bool goodEnoughDummyFitness(double bestFitness)
{
    return bestFitness >= 15.0;
}

double xorEvaluation(ExampleNetwork& individual)
{
    std::vector<std::array<double, 2> > inVals;
    
    inVals.push_back({0,0});
    inVals.push_back({0,1});
    inVals.push_back({1,0});
    inVals.push_back({1,1});
    
    std::array<double, 4> expectedOutput = {0, 1,1,0};
    std::array<double, 4> actualOutput;
    double diff = 0.0;
    for (int i=0; i<4; i++) {
        std::vector<double> input = std::vector<double>(inVals[i].begin(), inVals[i].end());
        std::pair<std::vector<double>, bool> eval =
        individual.simulateTillEquilibrium(input, 30);
        actualOutput[i] = eval.first.front();
        double d = expectedOutput[i] - eval.first.front();
        diff += fabs(d);
    }
    

    return (4-diff)*(4-diff);
}

bool xorFitnessSatisfied(double bestFitness)
{
    return bestFitness == 16.0;
}