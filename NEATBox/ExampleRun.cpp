//
//  ExampleRun.cpp
//  NEATBox
//
//  Created by Adeola Bannis on 1/7/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#include "ExampleRun.h"

#include "ExampleNetwork.h"
#include "ExampleGoals.h"

#include "NEATAlgorithm.cpp"
template class NEATAlgorithm<ExampleNetwork,Edge>;

ExampleNetwork createNetwork()
{
    return ExampleNetwork(2,1);
}

void runExample()
{
    
    NEATAlgorithm<ExampleNetwork, Edge> algo(100, 200, 10);
    algo.evaluationFunc = &xorEvaluation;
    algo.stopFunc = &xorFitnessSatisfied;
    algo.createInitialIndividual = &createNetwork;
    
    while (!algo.tick());
    ExampleNetwork *winner = algo.bestIndividual();
    std::cout << winner->display();
    
   // xorEvaluation(*winner);
    
}