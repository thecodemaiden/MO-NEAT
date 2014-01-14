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

BasicNN createXorNetwork()
{
    return BasicNN(2,1);
}

BasicNN createParityNetwork()
{
    return BasicNN(3,1);
}

void runExample()
{
    
    NEATPlusAlgorithm<BasicNN, Edge> algo(150, 200, 10);
//    algo.createInitialIndividual = &createXorNetwork;
//    algo.evaluationFunc = &xorEvaluation;
//   algo.stopFunc = &xorFitnessSatisfied;
    
    algo.createInitialIndividual = &createParityNetwork;
    algo.evaluationFunc = &parityEvaluation;
    algo.stopFunc = &parityFitnessSatisfied;
    
    algo.w_disjoint = 3.0;
    algo.w_excess = 2.0;
    algo.w_matching = 2.0;
    algo.w_matching_node = 3.0;
    
    algo.d_threshold = 10.0;
    
    while (!algo.tick());
    BasicNN *winner = algo.bestIndividual();
    std::cout << winner->display();
    std::cout << "Solution found in " << algo.getNumberOfIterations() << " generations.\n";
    
   // xorEvaluation(*winner);
    parityEvaluation(*winner);
}