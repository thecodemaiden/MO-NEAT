//
//  NNCommon.cpp
//  NEATBox
//
//  Created by Adeola Bannis on 2/13/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#include "NNCommon.h"
#include <cmath>
std::string activationFuncName(ActivationFunc f)
{
    //  std::ostringstream name;
    std::string name;
    switch (f) {
        case SIN_FUNC:
            name = "SIN";
            break;
        case GAUSSIAN_FUNC:
            name = "GAUSSIAN";
            break;
        case TANH_FUNC:
            name = "TANH";
            break;
        case STEP_FUNC:
            name = "THRESHOLD";
            break;
        default:
            break;
    }
    // return name.str();
    return name;
}

double applyActivationFunc(Node n, double inputSum)
{
    double output;
    
    switch (n.type) {
        case STEP_FUNC:
            if (inputSum > n.param1)
                output = n.activatedVal;
            else
                output = n.deactivatedVal;
            break;
        case TANH_FUNC:
            output = tanh(inputSum/exp(n.param1));
            break;
        case SIN_FUNC:
            output = sin(inputSum*exp(n.param1));
            break;
        case GAUSSIAN_FUNC:
        {
            double x = inputSum - n.param1;
            output = 2*exp(-x*x) - 1;
            break;
        }
        default:
            output = inputSum;
            break;
    }
    
    return output;
}