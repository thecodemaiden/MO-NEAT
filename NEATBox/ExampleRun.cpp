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

BasicNN *createXorIndicatorNetwork()
{
    return new BasicNN(2,2);
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
    
//    algocreateInitialIndividual = &createParityNetwork;
//    algo.evaluationFunc = &parityEvaluation;
//    algo.stopFunc = &parityFitnessSatisfied;
    
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

void runXorIndicatorExample()
{
    NEATPlusAlgorithm<BasicNN, Edge> *algo = new NEATPlusAlgorithm<BasicNN, Edge>(150, 200, 50);
    algo->createInitialIndividual = &createXorIndicatorNetwork;
    algo->evaluationFunc = &xorWithIndicatorEvaluation;
    algo->stopFunc = &xorWithIndicatorFitnessSatisfied;
    
    algo->w_disjoint = 3.0;
    algo->w_excess = 2.0;
    algo->w_matching = 2.0;
    algo->w_matching_node = 3.0;
    
    algo->d_threshold = 3.0;
    
    runAlgorithmToEnd(algo);
    
    BasicNN *winner = algo->bestIndividual();
    
    xorWithIndicatorEvaluation(winner);
    delete algo;

}

