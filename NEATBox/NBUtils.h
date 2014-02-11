//
//  NBUtils.h
//  NEATBox
//
//  Created by Adeola Bannis on 1/11/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#ifndef __NEATBox__NBUtils__
#define __NEATBox__NBUtils__

#include <iostream>

double normallyDistributed(double mean=0, double sd = 1);

long uniformlyDistributed(long upperBound);
double uniformProbability();
bool isUnreasonable(double d);
#endif /* defined(__NEATBox__NBUtils__) */
