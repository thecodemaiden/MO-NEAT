    //
//  ExampleRun.cpp
//  NEATBox
//
//  Created by Adeola Bannis on 1/7/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#include "ExampleRun.h"

#include "BasicNN.h"
#include "RecurrentNN.h"
#include "ExampleGoals.h"

#include "SPaNEAT.h"
#include "MONEAT.h"

#include "assert.h"

#include <vector>

MNIndividual *createXorNetwork()
{
    return new BasicNN(2,1);
}

MNIndividual *createParityNetwork()
{
    return new BasicNN(3,1);
}

MNIndividual *createMult3Network()
{
    return new BasicNN(4,2);
}

MNIndividual *createAllEvenNetwork()
{
    return new RecurrentNN(4,1);
}

void runXorExample(BaseNEAT *algo, long maxGenerations)
{
    algo->createInitialIndividual = &createXorNetwork;
    algo->evaluationFunctions.push_back(&xorEvaluation);
    
    algo->w_disjoint = 4.0;
    algo->w_excess = 3.0;
    algo->w_matching = 4.0;
    algo->w_matching_node = 3.0;
    
    algo->d_threshold = 3.0;
    
    for (long i=0; i<maxGenerations; i++)
        (algo->tick());
    
   BasicNN *winner = dynamic_cast<BasicNN *>(algo->optimalSolutions()[0]->individual);
//
 //   std::cout << winner->dotFormat();

    std::cout << "XOR score: " << xorEvaluation(winner) <<"\n\n";
}

void runParityExample(BaseNEAT *algo, long maxGenerations)
{
    algo->createInitialIndividual = &createParityNetwork;
    algo->evaluationFunctions.push_back(&parityEvaluation);
    
    algo->w_disjoint = 3.0;
    algo->w_excess = 2.0;
    algo->w_matching = 2.0;
    algo->w_matching_node = 3.0;
    
    algo->d_threshold = 3.0;
    
    for (long i=0; i<maxGenerations; i++)
        (algo->tick());
    
    BasicNN *winner = dynamic_cast<BasicNN *>(algo->optimalSolutions()[0]->individual);
    //
    std::cout << "Parity score: " << parityEvaluation(winner) <<"\n\n";

    
}

void runMult23SOTestTrain(BaseNEAT *algo, long maxGenerations)
{
    algo->createInitialIndividual = &createMult3Network;
    algo->evaluationFunctions.push_back(&trainMult23);
    
    algo->w_disjoint = 3.0;
    algo->w_excess = 2.0;
    algo->w_matching = 2.0;
    algo->w_matching_node = 3.0;
    
    algo->d_threshold = 3.0;
    
    for (long i=0; i<maxGenerations; i++)
        (algo->tick());
    
    BasicNN *winner = dynamic_cast<BasicNN *>(algo->optimalSolutions()[0]->individual);
    
    std::cout << "2 score: " << testMult2(winner) <<"\n";
    std::cout << "3 score: " << testMult3(winner) <<"\n\n";
}

void runMult23MOTestTrain(BaseNEAT *algo, long maxGenerations)
{
    algo->createInitialIndividual = &createMult3Network;
    algo->evaluationFunctions.push_back(&trainMult2);
    algo->evaluationFunctions.push_back(&trainMult3);
    
    algo->w_disjoint = 4.0;
    algo->w_excess = 4.0;
    algo->w_matching = 2.0;
    algo->w_matching_node = 2.0;
    
    algo->d_threshold = 4.0;
    
    for (long i=0; i<maxGenerations; i++)
        (algo->tick());
    
    std::vector<SystemInfo *>winners = algo->optimalSolutions();
    assert(winners.size() > 0);
    std::cout << winners.size() << " optimal solutions found.\n";
    
    SystemInfo *bestSystem = NULL;
    double min_sum = INFINITY;
    for (std::vector<SystemInfo *>::iterator it = winners.begin(); it!=winners.end(); it++) {
        double d= fabs((*it)->fitnesses[0] + (*it)->fitnesses[1]);
        
        if (d < min_sum) {
            bestSystem = *it;
            min_sum = d;
        }
    }
    
    if (bestSystem) {
        BasicNN *winner = dynamic_cast<BasicNN *>(bestSystem->individual);
        
        std::cout << "2 score: " << testMult2(winner) <<"\n";
        std::cout << "3 score: " << testMult3(winner) <<"\n\n";
    }
}

void runSequenceTest(BaseNEAT *algo, long maxGenerations)
{
    algo->createInitialIndividual = &createAllEvenNetwork;
    algo->evaluationFunctions.push_back(&isAllEven);
    algo->evaluationFunctions.push_back(&hasFewEdges);
    
    algo->w_disjoint = 4.0;
    algo->w_excess = 4.0;
    algo->w_matching = 2.0;
    algo->w_matching_node = 2.0;
    
    algo->d_threshold = 4.0;
    
    for (long i=0; i<maxGenerations; i++)
        (algo->tick());
    
    std::vector<SystemInfo *>winners = algo->optimalSolutions();
    assert(winners.size() > 0);
    
    
    // show all winners
    for (SystemInfo * i : winners) {
        RecurrentNN *winner = dynamic_cast<RecurrentNN *>(i->individual);
        std::cout << winner->dotFormat();
        std::cout << "Score: " << isAllEven(winner) << std::endl;
    }
    
//    RecurrentNN *winner = dynamic_cast<RecurrentNN *>(algo->optimalSolutions()[0]->individual);
//    
//    std::cout << winner->dotFormat();
//    std::cout << "Final score: " << isAllEven(winner) <<"\n";
//
//    
}

