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

    return individual->numberOfEdges();
}

bool goodEnoughDummyFitness(double bestFitness)
{
    return bestFitness >= 15.0;
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
    
   // int maxSteps = 50;
    std::array<double, numCases> expectedOutput = {1, -1, 1 ,-1,1,-1,1,-1};
    std::array<double, numCases> actualOutput;
    double diff = 0.0;
    for (int i=0; i<numCases; i++) {
        std::vector<double> input = std::vector<double>(inVals[i].begin(), inVals[i].end());
        std::vector<std::vector<double> > seqInputs;
        seqInputs.push_back(input);
        std::vector<std::vector<double> > seqOutputs = individual->simulateSequence(seqInputs, 2);

        double v = seqOutputs.front().front();
        
        actualOutput[i] = v;
        double d = expectedOutput[i] - v;
        diff += d*d;

    }
    
    return (4*numCases-diff)*(4*numCases-diff);
}

bool parityFitnessSatisfied(double bestFitness)
{
    //const double factor = 4.0;
    return bestFitness >= (975.0);
}



double xorEvaluation(BasicNN *individual)
{
    const int nCases = 4;
    const double yesVal = 1.0;
    const double noVal = -1.0;
    std::vector<std::array<double, 2> > inVals;
    
    inVals.push_back({noVal,noVal});
    inVals.push_back({noVal,yesVal});
    inVals.push_back({yesVal,noVal});
    inVals.push_back({yesVal,yesVal});
    
    std::array<double, nCases> expectedOutput = {-1, 1,1,-1};
    std::array<double, nCases> actualOutput;
    double diff = 0.0;

    for (int i=0; i<4; i++) {
        std::vector<double> input = std::vector<double>(inVals[i].begin(), inVals[i].end());
        std::vector<std::vector<double> > seqInputs;
        seqInputs.push_back(input);
        std::vector<std::vector<double> > seqOutputs = individual->simulateSequence(seqInputs, 1);
        std::vector<double> outputs = seqOutputs.front();
        actualOutput[i] = outputs.front();
    
        double d = expectedOutput[i] - outputs.front();
        diff += d*d;
        
    }
    
    // and diff runs from 0 to 4*n_cases;
    
    double finalScore = diff;
    
    const int factor = 4;
    
    return (factor*nCases - finalScore)*(factor*nCases - finalScore);
}

bool xorFitnessSatisfied(double bestFitness)
{
  //  const int factor = 2;

    return false;//bestFitness > (factor*4*factor*3);
}


double trainMult3(BasicNN *individual)
{
    
    /* leave out 6 cases for testing:
            1011 = 11, 1001 = 9, 1111 = 15, 0101 = 5, 0010 = 2, 0111 = 7
     */
    
    const int nCases = 10;
    std::vector<std::array<double, 4> > inVals;
    const double yesVal = 1.0;
    const double noVal = -1.0;
    inVals.push_back({noVal,noVal,noVal,noVal});
    inVals.push_back({noVal,noVal,noVal,yesVal});
    inVals.push_back({noVal,noVal,yesVal,yesVal});
    inVals.push_back({noVal,yesVal,noVal,noVal});
    inVals.push_back({noVal,yesVal,yesVal,noVal});
    inVals.push_back({yesVal,noVal,noVal,noVal});
    inVals.push_back({yesVal,noVal,yesVal,noVal});
    inVals.push_back({yesVal,yesVal,noVal,noVal});
    inVals.push_back({yesVal,yesVal,noVal,yesVal});
    inVals.push_back({yesVal,yesVal,yesVal,noVal});
 
    std::array<double, nCases> expectedOutput = {yesVal, noVal, yesVal, noVal, yesVal, noVal, noVal, yesVal, noVal, noVal};
    std::array<double, nCases> actualOutput;

    double diff = 0;
    for (int i=0; i<nCases; i++) {
        std::vector<double> input = std::vector<double>(inVals[i].begin(), inVals[i].end());
        std::vector<std::vector<double> > seqInputs;
        seqInputs.push_back(input);
        std::vector<std::vector<double> > seqOutputs = individual->simulateSequence(seqInputs, 1);
        std::vector<double> outputs = seqOutputs.front();
        actualOutput[i] = outputs.front();
        
        double d = expectedOutput[i] - outputs.front();
        diff += d*d;
        
    }

    const int factor = 4;
    
    return (factor*nCases - diff)*(factor*nCases - diff);
}

double testMult3(BasicNN *individual)
{
    /* the 6 cases for testing:
     1011 = 11, 1001 = 9, 1111 = 15, 0101 = 5, 0010 = 2, 0111 = 7
     */
    
    const int nCases = 6;
    std::vector<std::array<double, 4> > inVals;
    const double yesVal = 1.0;
    const double noVal = -1.0;
    
    inVals.push_back({noVal,noVal,yesVal,noVal});
    inVals.push_back({noVal,yesVal,yesVal,yesVal});
    inVals.push_back({noVal,yesVal,noVal,yesVal});
    inVals.push_back({yesVal,noVal,noVal,yesVal});
    inVals.push_back({yesVal,noVal,yesVal,yesVal});
    inVals.push_back({yesVal,yesVal,yesVal,yesVal});

    
    std::array<double, nCases> expectedOutput = {noVal, noVal, noVal, yesVal, noVal, yesVal};
    std::array<double, nCases> actualOutput;
    
    double diff = 0;
    for (int i=0; i<nCases; i++) {
        std::vector<double> input = std::vector<double>(inVals[i].begin(), inVals[i].end());
        std::vector<std::vector<double> > seqInputs;
        seqInputs.push_back(input);
        std::vector<std::vector<double> > seqOutputs = individual->simulateSequence(seqInputs, 1);
        std::vector<double> outputs = seqOutputs.front();
        actualOutput[i] = outputs.front();
        
        double d = expectedOutput[i] - outputs.front();
        diff += d*d;
        
    }
    
    const int factor = 4;
    
    return (factor*nCases - diff)*(factor*nCases - diff);
}
