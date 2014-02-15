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

double dummyEvaluation(MNIndividual *individual)
{

    return 1.0/(individual->numberOfEdges());
}


double parityEvaluation(MNIndividual *i)
{
    
    BasicNN *individual = dynamic_cast<BasicNN *>(i);
    if (!individual)
        return -INFINITY;
    
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
    
    return diff;
}




double xorEvaluation(MNIndividual *i)
{
    BasicNN *individual = dynamic_cast<BasicNN *>(i);
    if (!individual)
        return -INFINITY;
    
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
        
    double finalScore = diff;
    
    return finalScore;
}


double trainMult23(MNIndividual *i)
{
    return trainMult2(i)+trainMult3(i);
}


double trainMult2(MNIndividual *i)
{
    
    BasicNN *individual = dynamic_cast<BasicNN *>(i);
    if (!individual)
        return -INFINITY;
    
    const int nCases = 10;
    
    const double yesVal = 1.0;
    const double noVal = -1.0;
    
    std::vector<std::array<double, 4> > inVals;
    std::vector<std::array<double, 2> > expectedOutput;
    
    // 10 training cases:
    inVals.push_back({noVal,noVal,noVal,noVal}); //0
    expectedOutput.push_back({yesVal, yesVal});
    
    inVals.push_back({noVal,noVal,noVal,yesVal}); //1
    expectedOutput.push_back({noVal, noVal});
    
    inVals.push_back({noVal,noVal,yesVal,yesVal}); //3
    expectedOutput.push_back({yesVal, noVal});
    
    inVals.push_back({noVal,yesVal,noVal,noVal}); //4
    expectedOutput.push_back({noVal, yesVal});
    
    inVals.push_back({noVal,yesVal,yesVal,yesVal});//7
    expectedOutput.push_back({noVal, noVal});
    
    inVals.push_back({yesVal,noVal,noVal,noVal});//8
    expectedOutput.push_back({noVal, yesVal});
    
    inVals.push_back({yesVal,noVal,yesVal,noVal});//10
    expectedOutput.push_back({noVal, yesVal});
    
    inVals.push_back({yesVal,yesVal,noVal,noVal});//12
    expectedOutput.push_back({yesVal, yesVal});
    
    inVals.push_back({yesVal,yesVal,noVal,yesVal});//13
    expectedOutput.push_back({noVal, noVal});
    
    inVals.push_back({yesVal,yesVal,yesVal,yesVal});//15
    expectedOutput.push_back({yesVal, noVal});
    
    std::vector<std::array<double, 2> > actualOutput;
    
    //  int nCorrect = 0;
    
    double diff_3 = 0;
    double diff_2 = 0;
    for (int i=0; i<nCases; i++) {
        std::vector<double> input = std::vector<double>(inVals[i].begin(), inVals[i].end());
        std::vector<std::vector<double> > seqInputs;
        seqInputs.push_back(input);
        std::vector<std::vector<double> > seqOutputs = individual->simulateSequence(seqInputs, 1);
        std::vector<double> outputs = seqOutputs.front();
        actualOutput.push_back( {outputs[0], outputs[1]});
        
        double d_3 = expectedOutput[i][0] - outputs[0];
        double d_2 = expectedOutput[i][1] - outputs[1];
        diff_2 += d_2*d_2;
        diff_3 += d_3*d_3;
        
    }
    double diff = (diff_2);
//    const int factor = 4;
    
    return diff;//(factor*nCases - diff)*(factor*nCases - diff);
}

double trainMult3(MNIndividual *i)
{
    BasicNN *individual = dynamic_cast<BasicNN *>(i);
    if (!individual)
        return -INFINITY;
    
    const int nCases = 10;
    
    const double yesVal = 1.0;
    const double noVal = -1.0;
    
    std::vector<std::array<double, 4> > inVals;
    std::vector<std::array<double, 2> > expectedOutput;
    
    // 10 training cases:
    inVals.push_back({noVal,noVal,noVal,noVal}); //0
    expectedOutput.push_back({yesVal, yesVal});
    
    inVals.push_back({noVal,noVal,noVal,yesVal}); //1
    expectedOutput.push_back({noVal, noVal});
    
    inVals.push_back({noVal,noVal,yesVal,yesVal}); //3
    expectedOutput.push_back({yesVal, noVal});
    
    inVals.push_back({noVal,yesVal,noVal,noVal}); //4
    expectedOutput.push_back({noVal, yesVal});
    
    inVals.push_back({noVal,yesVal,yesVal,yesVal});//7
    expectedOutput.push_back({noVal, noVal});
    
    inVals.push_back({yesVal,noVal,noVal,noVal});//8
    expectedOutput.push_back({noVal, yesVal});
    
    inVals.push_back({yesVal,noVal,yesVal,noVal});//10
    expectedOutput.push_back({noVal, yesVal});
    
    inVals.push_back({yesVal,yesVal,noVal,noVal});//12
    expectedOutput.push_back({yesVal, yesVal});
    
    inVals.push_back({yesVal,yesVal,noVal,yesVal});//13
    expectedOutput.push_back({noVal, noVal});
    
    inVals.push_back({yesVal,yesVal,yesVal,yesVal});//15
    expectedOutput.push_back({yesVal, noVal});
    
    
    
    std::vector<std::array<double, 2> > actualOutput;
    
    //  int nCorrect = 0;
    
    double diff_3 = 0;
    double diff_2 = 0;
    for (int i=0; i<nCases; i++) {
        std::vector<double> input = std::vector<double>(inVals[i].begin(), inVals[i].end());
        std::vector<std::vector<double> > seqInputs;
        seqInputs.push_back(input);
        std::vector<std::vector<double> > seqOutputs = individual->simulateSequence(seqInputs, 1);
        std::vector<double> outputs = seqOutputs.front();
        actualOutput.push_back( {outputs[0], outputs[1]});
        
        double d_3 = expectedOutput[i][0] - outputs[0];
        double d_2 = expectedOutput[i][1] - outputs[1];
        diff_2 += d_2*d_2;
        diff_3 += d_3*d_3;
        
    }
    double diff = (diff_3);
   // const int factor = 4;
    
    return diff;
}

double testMult2(BasicNN *individual)
{

    
    const int nCases = 6;
    const double yesVal = 1.0;
    const double noVal = -1.0;
    
    std::vector<std::array<double, 4> > inVals;
    std::vector<std::array<double, 2> > expectedOutput;
    
    // 6 test cases
    inVals.push_back({noVal,noVal,yesVal,noVal}); // 2
    expectedOutput.push_back({noVal, yesVal});
    
    inVals.push_back({noVal,yesVal,noVal,yesVal}); // 5
    expectedOutput.push_back({noVal, noVal});
    
    inVals.push_back({noVal,yesVal,yesVal,noVal}); // 6
    expectedOutput.push_back({yesVal, yesVal});
    
    inVals.push_back({yesVal,noVal,noVal,yesVal}); // 9
    expectedOutput.push_back({yesVal, noVal});
    
    inVals.push_back({yesVal,noVal,yesVal,yesVal}); // 11
    expectedOutput.push_back({noVal,noVal});
    
    inVals.push_back({yesVal,yesVal,yesVal,noVal}); // 14
    expectedOutput.push_back({noVal, yesVal});
    
    
    std::vector<std::array<double, 2> > actualOutput;
    
    double diff_3 = 0;
    double diff_2 = 0;
    for (int i=0; i<nCases; i++) {
        std::vector<double> input = std::vector<double>(inVals[i].begin(), inVals[i].end());
        std::vector<std::vector<double> > seqInputs;
        seqInputs.push_back(input);
        std::vector<std::vector<double> > seqOutputs = individual->simulateSequence(seqInputs, 1);
        std::vector<double> outputs = seqOutputs.front();
        actualOutput.push_back({outputs[0], outputs[1]});
        
        double d_3 = expectedOutput[i][0] - outputs[0];
        double d_2 = expectedOutput[i][1] - outputs[1];
        diff_2 += d_2*d_2;
        diff_3 += d_3*d_3;
        
    }
    double diff = diff_2;
    
   // const int factor = 4;
    
    return diff;
}

double testMult3(BasicNN *individual)
{
    
    const int nCases = 6;
    const double yesVal = 1.0;
    const double noVal = -1.0;
    
    std::vector<std::array<double, 4> > inVals;
    std::vector<std::array<double, 2> > expectedOutput;
    
    // 6 test cases
    inVals.push_back({noVal,noVal,yesVal,noVal}); // 2
    expectedOutput.push_back({noVal, yesVal});
    
    inVals.push_back({noVal,yesVal,noVal,yesVal}); // 5
    expectedOutput.push_back({noVal, noVal});
    
    inVals.push_back({noVal,yesVal,yesVal,noVal}); // 6
    expectedOutput.push_back({yesVal, yesVal});
    
    inVals.push_back({yesVal,noVal,noVal,yesVal}); // 9
    expectedOutput.push_back({yesVal, noVal});
    
    inVals.push_back({yesVal,noVal,yesVal,yesVal}); // 11
    expectedOutput.push_back({noVal,noVal});
    
    inVals.push_back({yesVal,yesVal,yesVal,noVal}); // 14
    expectedOutput.push_back({noVal, yesVal});
    
    
    std::vector<std::array<double, 2> > actualOutput;
    
    double diff_3 = 0;
    double diff_2 = 0;
    for (int i=0; i<nCases; i++) {
        std::vector<double> input = std::vector<double>(inVals[i].begin(), inVals[i].end());
        std::vector<std::vector<double> > seqInputs;
        seqInputs.push_back(input);
        std::vector<std::vector<double> > seqOutputs = individual->simulateSequence(seqInputs, 1);
        std::vector<double> outputs = seqOutputs.front();
        actualOutput.push_back({outputs[0], outputs[1]});
        
        double d_3 = expectedOutput[i][0] - outputs[0];
        double d_2 = expectedOutput[i][1] - outputs[1];
        diff_2 += d_2*d_2;
        diff_3 += d_3*d_3;
        
    }
    double diff = diff_3;
    
   // const int factor = 4;
    
    return diff;
}

