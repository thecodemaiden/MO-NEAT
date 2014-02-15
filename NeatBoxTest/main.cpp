//
//  main.cpp
//  NeatBoxTest
//
//  Created by Adeola Bannis on 1/8/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#include <iostream>
#include "ExampleRun.h"
#include <NEATBox/MONEAT.h>
#include <NEATBox/SPaNEAT.h>

int main(int argc, const char * argv[])
{
    std::cout << "--- NSGA-II ---\n";
    for (int i=0; i<10; i++) {
        MONEAT *m = new MONEAT(100,100,50);
        
        runMult23MOTestTrain(m);
        delete m;
        
    }
    std::cout << "--- SPEA2 ---\n";
    
    for (int i=0; i<10; i++) {
        SPaNEAT *s = new SPaNEAT(100,25,100);
        
        runMult23MOTestTrain(s);
        delete s;
        
    }
    
    return 0;
}

