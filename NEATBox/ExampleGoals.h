//
//  ExampleGoals.h
//  NEATBox
//
//  Created by Adeola Bannis on 1/7/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#ifndef __NEATBox__ExampleGoals__
#define __NEATBox__ExampleGoals__

#include <iostream>
#include "BasicNN.h"


double dummyEvaluation(MNIndividual *individual);

double xorEvaluation(MNIndividual *individual);

double parityEvaluation(MNIndividual *individual);

// IS THE INPUT A MULTIPLE OF 3?2? 
double trainMult23(MNIndividual *individual);

double trainMult2(MNIndividual *individual);
double trainMult3(MNIndividual *individual);

double testMult2(BasicNN *individual);
double testMult3(BasicNN *individual);

// for when we want to limt complexity
double hasFewNodes(MNIndividual *individual);
double hasFewEdges(MNIndividual *individual);

// --- Sequence test ---

// have all inputs so far been even numbers?
double isAllEven(MNIndividual *individual);
//double trainAllEven(BasicNN * individual);
//double testAllEven(BasicNN * individual);

#endif /* defined(__NEATBox__ExampleGoals__) */
