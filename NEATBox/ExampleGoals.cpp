//
//  ExampleGoals.cpp
//  NEATBox
//
//  Created by Adeola Bannis on 1/7/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#include "ExampleGoals.h"
#include "BasicNN.h"
#include <array>
#include <cmath>

double dummyEvaluation(BasicNN *individual)
{
    SimReturn eval =
    individual->simulateTillEquilibrium(std::vector<double>(), 10);
    
    for (long i=0; i<eval.outputs.size(); i++) {
        std::cout <<"\t" << i << ": " << eval.outputs[i] <<"\n";
    }
    std::cout << "\tsteady: " << (eval.steady ? "y" : "n") << "\n";
    
    return individual->numberOfEdges();
}

bool goodEnoughDummyFitness(double bestFitness)
{
    return bestFitness >= 15.0;
}

double xorEvaluation(BasicNN& individual)
{
    std::vector<std::array<double, 2> > inVals;
    
    inVals.push_back({0,0});
    inVals.push_back({0,1});
    inVals.push_back({1,0});
    inVals.push_back({1,1});
    
    int maxSteps = 30;
    
    std::array<double, 4> expectedOutput = {0, 1,1,0};
    std::array<double, 4> actualOutput;
    double diff = 0.0;
    for (int i=0; i<4; i++) {
        std::vector<double> input = std::vector<double>(inVals[i].begin(), inVals[i].end());
        SimReturn eval =
        individual.simulateTillEquilibrium(input, maxSteps);
        actualOutput[i] = eval.outputs.front();
        double d = expectedOutput[i] - eval.outputs.front();
        diff += fabs(d);
        
        // penalize unstable networks
        if (!eval.steady)
            diff += 1.0;
        
       // diff += (double)eval.steps/maxSteps;
    }
    

    return (9-diff)*(9-diff);
}

bool xorFitnessSatisfied(double bestFitness)
{
    return bestFitness > 80.0;
}


const int numCases = 8;
double parityEvaluation(BasicNN *individual)
{
    std::vector<std::array<double, 3> > inVals;
    
    inVals.push_back({0,0,0});
    inVals.push_back({0,0,1});
    inVals.push_back({0,1,0});
    inVals.push_back({0,1,1});
    inVals.push_back({1,0,0});
    inVals.push_back({1,0,1});
    inVals.push_back({1,1,0});
    inVals.push_back({1,1,1});
    
    int maxSteps = 50;
    std::array<double, numCases> expectedOutput = {1, -1, -1 ,1,-1,1,1,-1};
    std::array<double, numCases> actualOutput;
    double diff = 0.0;
    for (int i=0; i<numCases; i++) {
        std::vector<double> input = std::vector<double>(inVals[i].begin(), inVals[i].end());
        SimReturn eval =
        individual->simulateTillEquilibrium(input, maxSteps);
        
        double v = eval.outputs.front();
        
        actualOutput[i] = v;
        double d = expectedOutput[i] - v;
        diff += fabs(d);
    
        
        // penalize unstable networks
        //if (!eval.steady)
        //    diff += 1.0;
        
         diff += (double)eval.steps/maxSteps;
    }
    
    
    return (2*numCases-diff)*(2*numCases-diff);
}

bool parityFitnessSatisfied(double bestFitness)
{
    return bestFitness >= (4*numCases*numCases);
}