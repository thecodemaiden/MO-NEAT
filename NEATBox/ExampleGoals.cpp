 //
//  ExampleGoals.cpp
//  NEATBox
//
//  Created by Adeola Bannis on 1/7/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#include "ExampleGoals.h"
#include "BasicNN.h"
#include "RecurrentNN.h"
#include <array>
#include <cmath>

static const double yesVal = 1.0;
static const double noVal = -1.0;

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
        
        std::vector<double> outputs = individual->simulateSequence(input);
        double v= outputs.front();
        
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
    
    static double bestFitness = INFINITY;
    
    const int nCases = 4;
 
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
 
        std::vector<double> outputs = individual->simulateSequence(input);
        actualOutput[i] = outputs.front();
    
        double d = expectedOutput[i] - outputs.front();
        diff += d*d;
        
    }
        
    double finalScore = diff;
    
    if (finalScore < bestFitness) {
        bestFitness = finalScore;
        std::cout << "BEST: " << bestFitness <<"\n";
    }
    
    
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
        
        std::vector<double> outputs = individual->simulateSequence(input);

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
        
        std::vector<double> outputs = individual->simulateSequence(input);
        
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
        
        std::vector<double> outputs = individual->simulateSequence(input);
        
        actualOutput.push_back( {outputs[0], outputs[1]});
        
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
        
        std::vector<double> outputs = individual->simulateSequence(input);
        
        actualOutput.push_back( {outputs[0], outputs[1]});
        
        double d_3 = expectedOutput[i][0] - outputs[0];
        double d_2 = expectedOutput[i][1] - outputs[1];
        diff_2 += d_2*d_2;
        diff_3 += d_3*d_3;
        
    }
    double diff = diff_3;
    
   // const int factor = 4;
    
    return diff;
}

#pragma mark - Recurrent

static std::vector<double> bin_rep4(int n)
{
    if (n > 15)
        n = 15;
    if (n < 0)
        n = 0;
    
    std::vector<double> toReturn;
    toReturn.push_back(((n>>3)&1) ? yesVal : noVal);
    toReturn.push_back(((n>>2)&1) ? yesVal : noVal);
    toReturn.push_back(((n>>1)&1) ? yesVal : noVal);
    toReturn.push_back((n&1) ? yesVal : noVal);
    
    return toReturn;
}

double isAllEven(MNIndividual *i)
{

    RecurrentNN *individual = dynamic_cast<RecurrentNN *>(i);
    if (!individual)
        return -INFINITY;
    
    double diff = 0;
    
    int allEvenRun = arc4random_uniform(10)+2;

    for (int i=0; i<10; i++) {
      //  bool allEvens = (uniformProbability() > 0.5);
        std::vector<std::vector<double> >seqInputs;
        bool allowOdd = i%allEvenRun == 0;
        std::vector<double> expectedOutputs;
        bool hadOdd = false;
        // a sequence of 3-5 inputs
        int nInputs = 3 + arc4random_uniform(3);
        for (int j=0; j<nInputs; j++) {
            int n;
            if (allowOdd) {
                n = arc4random_uniform(16);
                if (n %2 == 1) {
                    hadOdd = true;
                }
            } else {
                n = arc4random_uniform(8)*2;
            }
            if (!hadOdd)
                expectedOutputs.push_back(yesVal);
            else
                expectedOutputs.push_back(noVal);
            
            std::vector<double> next = bin_rep4(n);
            seqInputs.push_back(next);
        }
        std::vector<std::vector<double> >seqOutputs = individual->simulateSequence(seqInputs);
        for (int j=0; j<seqOutputs.size(); j++) {
            double retVal = seqOutputs[j].front();
            double d = (retVal - expectedOutputs[j]);
            if (d != 0) {
                d = fabs(d);
            }
            diff += d;
        }
    }
    
    return diff;
}

double trainAllEven(BasicNN * individual)
{
//TBD
    
    return 0.0;
}

double testAllEven(BasicNN * individual)
{
    //TBD
    return 0.0;
}

