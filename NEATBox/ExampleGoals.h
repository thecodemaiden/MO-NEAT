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
bool goodEnoughDummyFitness(double bestFitness);

double xorEvaluation(MNIndividual *individual);
bool xorFitnessSatisfied(double bestFitness);

double parityEvaluation(MNIndividual *individual);
bool parityFitnessSatisfied(double bestFitness);

// IS THE INPUT A MULTIPLE OF 3?2? 
double trainMult23(MNIndividual *individual);
double testMult23(BasicNN *individual);

double trainMult2(MNIndividual *individual);
double trainMult3(MNIndividual *individual);

double testMult2(BasicNN *individual);
double testMult3(BasicNN *individual);

#endif /* defined(__NEATBox__ExampleGoals__) */
