//
//  NBUtils.cpp
//  NEATBox
//
//  Created by Adeola Bannis on 1/11/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#include "NBUtils.h"
#include <cmath>
double normallyDistributed(double mean, double sd)
{
    static double n2 = 0.0;
    static int n2_cached = 0;
    
    if (!n2_cached) {
        double x, y, r;
        do {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;
            r = x*x + y*y;
        } while (r == 0.0 || r > 1.0);
                
        double d = sqrt(-2.0*log(r)/r);
        double n1 = x*d;
        n2 = y*d;
        double result = n1*sd + mean;
        n2_cached = 1;
        return result;
        
    } else {
        n2_cached = 0;
        return n2*sd + mean;
    }
}