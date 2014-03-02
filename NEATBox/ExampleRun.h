//
//  ExampleRun.h
//  NEATBox
//
//  Created by Adeola Bannis on 1/7/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#ifndef __NEATBox__ExampleRun__
#define __NEATBox__ExampleRun__

#include <iostream>
#include <NEATBox/BaseNeat.h>

void runXorExample(BaseNEAT *algo, long maxGenerations=100);
void runParityExample(BaseNEAT *algo, long maxGenerations=100);

void runMult23SOTestTrain(BaseNEAT *algo, long maxGenerations=100);

void runMult23MOTestTrain(BaseNEAT *algo, long maxGenerations=100);

void runSequenceTest(BaseNEAT *algo, long maxGenerations=100);
#endif /* defined(__NEATBox__ExampleRun__) */
