//
//  SPaNEAT.h
//  NEATBox
//
//  Created by Adeola Bannis on 2/9/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#ifndef __NEATBox__SPaNEAT__
#define __NEATBox__SPaNEAT__

#include <iostream>


//
//  SPaNEAT.h
//  SystemGenerator
//
//  Created by Adeola Bannis on 11/4/13.
//
//

#ifndef __SystemGenerator__SPaNEAT__
#define __SystemGenerator__SPaNEAT__

#include <vector>
#include <fstream>

#include <cmath>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <wordexp.h>


#include "NEATCommon.h"
#include "BaseNEAT.h"
struct SPSystemInfo : public SystemInfo
{
    std::vector<double> distances;
    
    SPSystemInfo(MNIndividual *i):SystemInfo(i) {}
    
    SPSystemInfo(const SPSystemInfo &other)
    :SystemInfo(other),
    distances(other.distances)
    {}
};

class SPaNEAT :public BaseNEAT {
protected:
    std::vector<SPSystemInfo *> population;
    std::vector<SPSystemInfo *> archive;
        
    long archiveSize;
    
    void spawnNextGeneration(); // recombine species to get enough children then mutate each one
    
    // finds the Pareto fronts of individuals
    void rankSystems();
    
    // overriden functions cannot be called in constructors, so this is called on the first tick();
    virtual void prepareInitialPopulation();
    
public:
    SPaNEAT(long populationSize, long archiveSize);
    ~SPaNEAT();
    
    void tick();
    std::vector<SystemInfo *> optimalSolutions();

};


#endif /* defined(__SystemGenerator__SPaNEAT__) */



#endif /* defined(__NEATBox__SPaNEAT__) */
