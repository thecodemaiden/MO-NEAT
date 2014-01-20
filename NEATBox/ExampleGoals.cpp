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

double xorEvaluation(BasicNN *individual)
{
    std::vector<std::array<double, 2> > inVals;
    
    inVals.push_back({-1,-1});
    inVals.push_back({-1,1});
    inVals.push_back({1,-1});
    inVals.push_back({1,1});
    
    int maxSteps = 30;
    
    std::array<double, 4> expectedOutput = {-1, 1,1,-1};
    std::array<double, 4> actualOutput;
    double diff = 0.0;
    for (int i=0; i<4; i++) {
        std::vector<double> input = std::vector<double>(inVals[i].begin(), inVals[i].end());
        SimReturn eval =
        individual->simulateTillEquilibrium(input, maxSteps);
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


double parityEvaluation(BasicNN *individual)
{
    const int numCases = 8;
    std::vector<std::array<double, 3> > inVals;
    const double yesVal = 1.0;
    const double noVal = -1.0;
    inVals.push_back({noVal,noVal,noVal});
    inVals.push_back({noVal,noVal,yesVal});
    inVals.push_back({noVal,yesVal,noVal});
    inVals.push_back({noVal,yesVal,yesVal});
    inVals.push_back({yesVal,noVal,noVal});
    inVals.push_back({yesVal,noVal,yesVal});
    inVals.push_back({yesVal,yesVal,noVal});
    inVals.push_back({yesVal,yesVal,yesVal});
    
    int maxSteps = 50;
    std::array<double, numCases> expectedOutput = {1, -1, -1 ,1,-1,-1,1,-1};
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
        if (!eval.steady)
            diff += 1.0;
        
        // diff += (double)eval.steps/maxSteps;
    }
    
    
    return (2*numCases-diff)*(2*numCases-diff);
}

bool parityFitnessSatisfied(double bestFitness)
{
    return bestFitness >= (4*8*8);
}



double xorWithIndicatorEvaluation(BasicNN *individual)
{
    const int nCases = 4;
    const double yesVal = 1.0;
    const double noVal = -1.0;
    std::vector<std::array<double, 2> > inVals;
    
    inVals.push_back({noVal,noVal});
    inVals.push_back({noVal,yesVal});
    inVals.push_back({yesVal,noVal});
    inVals.push_back({yesVal,yesVal});
    
    int maxSteps = std::max(5, (int)individual->numberOfEdges());
    
    std::array<double, nCases> expectedOutput = {-1, 1,1,-1};
    std::array<double, nCases> actualOutput;
    double diff = 0.0;
    double indicatorDiff = 0.0;
    double instabilityDiff = 0.0;
    for (int i=0; i<4; i++) {
        std::vector<double> input = std::vector<double>(inVals[i].begin(), inVals[i].end());
        SimReturn eval =
        individual->simulateTillEquilibrium(input, maxSteps);
        actualOutput[i] = eval.outputs.front();
        double d = expectedOutput[i] - eval.outputs.front();
        diff += fabs(d)*fabs(d);
        
        double indicatorScore = (eval.steady ? 1.0 : -1.0) - eval.outputs.back();
        //indicatorScore = indicatorScore*indicatorScore;
        
        indicatorDiff += fabs(indicatorScore);
        
        //instabilityDiff += (double)eval.steps/maxSteps;
        if (!eval.steady)
            instabilityDiff += 1.0;
    }
    
    // so indicatorDiff runs from 0 to 4*n_cases
    // and diff runs from 0 to 2*n_cases;
    // instabilityDiff runs from 0 to n_cases;
    
    double finalScore = diff*1.0 + indicatorDiff*1.0 + instabilityDiff*2.0;
    
    const int factor = 8;
    
    return (factor*nCases - finalScore)*(factor*nCases - finalScore);
}

bool xorWithIndicatorFitnessSatisfied(double bestFitness)
{
    const int factor = 8;

    return false;//bestFitness > (factor*4*factor*3);
}