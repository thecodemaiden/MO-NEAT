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


double dummyEvaluation(BasicNN& individual);
bool goodEnoughDummyFitness(double bestFitness);

double xorEvaluation(BasicNN& individual);
bool xorFitnessSatisfied(double bestFitness);

#endif /* defined(__NEATBox__ExampleGoals__) */
