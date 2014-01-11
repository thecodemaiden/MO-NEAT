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

#include "NEATAlgorithm.cpp"
template class NEATAlgorithm<BasicNN,Edge>;

BasicNN createNetwork()
{
    return BasicNN(2,1);
}

void runExample()
{
    
    NEATAlgorithm<BasicNN, Edge> algo(150, 200, 10);
    algo.evaluationFunc = &xorEvaluation;
    algo.stopFunc = &xorFitnessSatisfied;
    algo.createInitialIndividual = &createNetwork;
    
    algo.w_disjoint = 2.0;
    algo.w_excess = 3.0;
    algo.w_matching = 2.0;
    
    while (!algo.tick());
    BasicNN *winner = algo.bestIndividual();
    std::cout << winner->display();
    std::cout << "Solution found in " << algo.getNumberOfIterations() << " generations.\n";
    
}