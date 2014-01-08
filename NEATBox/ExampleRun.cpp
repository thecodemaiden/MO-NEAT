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

#include "NEATAlgorithm.h"

void runExample()
{
    
    NEATAlgorithm<ExampleNetwork, Edge> *algo = new NEATAlgorithm<ExampleNetwork, Edge>(100, 10, 5);
//    algo->evaluationFunc = &dummyEvaluation;
//    algo->stopFunc = &goodEnoughDummyFitness;
//    
//    while (algo->tick());
//    ExampleNetwork winner = algo->bestIndividual();
//    winner.display();

    
    
}