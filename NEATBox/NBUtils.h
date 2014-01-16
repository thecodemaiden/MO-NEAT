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

// if both are given, the first is a lower bound and the second is an upper bound
// else
long uniformlyDistributed(long upperBound);
bool isUnreasonable(double d);
#endif /* defined(__NEATBox__NBUtils__) */
