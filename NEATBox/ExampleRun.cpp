    //
//  ExampleRun.cpp
//  NEATBox
//
//  Created by Adeola Bannis on 1/7/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#include "ExampleRun.h"

#include "BasicNN.h"
#include "ExampleGoals.h"

#include "MONEAT.cpp"

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

void runXorExample()
{
    MONEAT *algo = new MONEAT(150, 200, 50);
    algo->createInitialIndividual = &createXorNetwork;
    algo->evaluationFunctions.push_back(&xorEvaluation);
    
    algo->w_disjoint = 3.0;
    algo->w_excess = 2.0;
    algo->w_matching = 2.0;
    algo->w_matching_node = 3.0;
    
    algo->d_threshold = 3.0;
    
    while (!algo->tick());
    
//    BasicNN *winner = algo->bestIndividual();
//    
//    std::cout << "XOR score: " << xorEvaluation(winner) <<"\n";
    delete algo;

}

void runParityExample()
{
    MONEAT *algo = new MONEAT (150, 200, 50);
    algo->createInitialIndividual = &createParityNetwork;
    algo->evaluationFunctions.push_back(&parityEvaluation);
    
    algo->w_disjoint = 3.0;
    algo->w_excess = 2.0;
    algo->w_matching = 2.0;
    algo->w_matching_node = 3.0;
    
    algo->d_threshold = 3.0;
    
    while (!algo->tick());
    
//    BasicNN *winner = algo->bestIndividual();
//    
//    parityEvaluation(winner);
    delete algo;
    
}

void runMult23SOTestTrain()
{
    MONEAT  *algo = new MONEAT (150, 200, 50);
    algo->createInitialIndividual = &createMult3Network;
    algo->evaluationFunctions.push_back(&trainMult23);
    
    algo->w_disjoint = 3.0;
    algo->w_excess = 2.0;
    algo->w_matching = 2.0;
    algo->w_matching_node = 3.0;
    
    algo->d_threshold = 3.0;
    
    while (!algo->tick());
    
  //  BasicNN *winner = algo->bestIndividual();
    
 //   std::cout << "Test score: " << testMult23(winner) <<"\n";
    delete algo;
}

void runMult23MOTestTrain()
{
    MONEAT  *algo = new MONEAT (25, 200, 50);
    algo->createInitialIndividual = &createMult3Network;
    algo->evaluationFunctions.push_back(&trainMult2);
    algo->evaluationFunctions.push_back(&trainMult3);
    
    algo->w_disjoint = 4.0;
    algo->w_excess = 4.0;
    algo->w_matching = 2.0;
    algo->w_matching_node = 2.0;
    
    algo->d_threshold = 2.0;
    
    while (!algo->tick());
    
    std::vector<SystemInfo *>winners = algo->bestIndividuals;

    std::cout << winners.size() << " optimal solutions found.\n";
    
    SystemInfo *bestSystem = NULL;
    double min_diff = INFINITY;
    for (std::vector<SystemInfo *>::iterator it = winners.begin(); it!=winners.end(); it++) {
        double d= fabs((*it)->fitnesses[0] - (*it)->fitnesses[1]);
        if (d < min_diff) {
            bestSystem = *it;
            min_diff = d;
        }
    }
    
    BasicNN *winner = dynamic_cast<BasicNN *>(bestSystem->individual);
    std::cout << "2 score: " << testMult2(winner) <<"\n";
    std::cout << "3 score: " << testMult3(winner) <<"\n";
    
    delete algo;
}


