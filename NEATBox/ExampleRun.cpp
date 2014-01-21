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

#include "NEATPlusAlgorithm.cpp"
template class NEATPlusAlgorithm<BasicNN,Edge>;

BasicNN *createXorNetwork()
{
    return new BasicNN(2,1);
}

BasicNN *createParityNetwork()
{
    return new BasicNN(3,1);
}

BasicNN *createMult3Network()
{
    return new BasicNN(4,1);
}

void runAlgorithmToEnd(NEATPlusAlgorithm<BasicNN, Edge> *algo) {
    while (!algo->tick());
    BasicNN *winner = algo->bestIndividual();
    std::cout << "Solution found in " << algo->getNumberOfIterations() << " generations.\n\n";
    std::cout << winner->dotFormat();
}


void runXorExample()
{
    NEATPlusAlgorithm<BasicNN, Edge> *algo = new NEATPlusAlgorithm<BasicNN, Edge>(150, 200, 50);
    algo->createInitialIndividual = &createXorNetwork;
    algo->evaluationFunc = &xorEvaluation;
    algo->stopFunc = &xorFitnessSatisfied;
    
    algo->w_disjoint = 3.0;
    algo->w_excess = 2.0;
    algo->w_matching = 2.0;
    algo->w_matching_node = 3.0;
    
    algo->d_threshold = 3.0;
    
    runAlgorithmToEnd(algo);
    
    BasicNN *winner = algo->bestIndividual();
    
    xorEvaluation(winner);
    delete algo;

}

void runParityExample()
{
    NEATPlusAlgorithm<BasicNN, Edge> *algo = new NEATPlusAlgorithm<BasicNN, Edge>(150, 200, 50);
    algo->createInitialIndividual = &createParityNetwork;
    algo->evaluationFunc = &parityEvaluation;
    algo->stopFunc = &parityFitnessSatisfied;
    
    algo->w_disjoint = 3.0;
    algo->w_excess = 2.0;
    algo->w_matching = 2.0;
    algo->w_matching_node = 3.0;
    
    algo->d_threshold = 3.0;
    
    runAlgorithmToEnd(algo);
    
    BasicNN *winner = algo->bestIndividual();
    
    parityEvaluation(winner);
    delete algo;
    
}

void runMult3TestTrain()
{
    NEATPlusAlgorithm<BasicNN, Edge> *algo = new NEATPlusAlgorithm<BasicNN, Edge>(150, 200, 50);
    algo->createInitialIndividual = &createMult3Network;
    algo->evaluationFunc = &trainMult3;
    
    algo->w_disjoint = 3.0;
    algo->w_excess = 2.0;
    algo->w_matching = 2.0;
    algo->w_matching_node = 3.0;
    
    algo->d_threshold = 3.0;
    
    runAlgorithmToEnd(algo);
    
    BasicNN *winner = algo->bestIndividual();
    
    std::cout << "Test score: " << testMult3(winner) <<"\n";
    delete algo;
}

