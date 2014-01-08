//
//  ExampleGoals.cpp
//  NEATBox
//
//  Created by Adeola Bannis on 1/7/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#include "ExampleGoals.h"

double dummyEvaluation(ExampleNetwork& individual)
{
    return individual.numberOfEdges();
}

bool goodEnoughDummyFitness(double bestFitness)
{
    return bestFitness >= 15.0;
}