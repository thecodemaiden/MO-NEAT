//
//  NNCommon.cpp
//  NEATBox
//
//  Created by Adeola Bannis on 2/13/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#include "NNCommon.h"
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