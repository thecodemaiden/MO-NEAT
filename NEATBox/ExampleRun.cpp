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
template class MONEAT<BasicNN,Edge>;

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
    return new BasicNN(4,2);
}

void runAlgorithmToEnd(MONEAT<BasicNN, Edge> *algo) {
    while (!algo->tick());
    BasicNN *winner = algo->bestIndividual();
    std::cout << "Solution found in " << algo->getNumberOfIterations() << " generations.\n\n";
    std::cout << winner->dotFormat();
}


void runXorExample()
{
    MONEAT<BasicNN, Edge> *algo = new MONEAT<BasicNN, Edge>(150, 200, 50);
    algo->createInitialIndividual = &createXorNetwork;
    algo->evaluationFunctions.push_back(&xorEvaluation);
    
    algo->w_disjoint = 3.0;
    algo->w_excess = 2.0;
    algo->w_matching = 2.0;
    algo->w_matching_node = 3.0;
    
    algo->d_threshold = 3.0;
    
    runAlgorithmToEnd(algo);
    
    BasicNN *winner = algo->bestIndividual();
    
    std::cout << "XOR score: " << xorEvaluation(winner) <<"\n";
    delete algo;

}

void runParityExample()
{
    MONEAT<BasicNN, Edge> *algo = new MONEAT<BasicNN, Edge>(150, 200, 50);
    algo->createInitialIndividual = &createParityNetwork;
    algo->evaluationFunctions.push_back(&parityEvaluation);
    
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

void runMult23SOTestTrain()
{
    MONEAT<BasicNN, Edge> *algo = new MONEAT<BasicNN, Edge>(150, 200, 50);
    algo->createInitialIndividual = &createMult3Network;
    algo->evaluationFunctions.push_back(&trainMult23);
    
    algo->w_disjoint = 3.0;
    algo->w_excess = 2.0;
    algo->w_matching = 2.0;
    algo->w_matching_node = 3.0;
    
    algo->d_threshold = 3.0;
    
    runAlgorithmToEnd(algo);
    
    BasicNN *winner = algo->bestIndividual();
    
    std::cout << "Test score: " << testMult23(winner) <<"\n";
    delete algo;
}

void runMult23MOTestTrain()
{
    MONEAT<BasicNN, Edge> *algo = new MONEAT<BasicNN, Edge>(150, 200, 50);
    algo->createInitialIndividual = &createMult3Network;
    algo->evaluationFunctions.push_back(&trainMult2);
    algo->evaluationFunctions.push_back(&trainMult3);
    
    algo->w_disjoint = 3.0;
    algo->w_excess = 2.0;
    algo->w_matching = 2.0;
    algo->w_matching_node = 3.0;
    
    algo->d_threshold = 3.0;
    
    runAlgorithmToEnd(algo);
    
    BasicNN *winner = algo->bestIndividual();
    
    std::cout << "2 score: " << testMult2(winner) <<"\n";
    std::cout << "3 score: " << testMult3(winner) <<"\n";
    delete algo;
}


