//
//  NEATCommon.h
//  NEATBox
//
//  Created by Adeola Bannis on 2/9/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#ifndef NEATBox_NEATCommon_h
#define NEATBox_NEATCommon_h

#include <map>
#include "MNIndividual.h"
#include <cmath>

struct SystemInfo {
public:
    MNIndividual *individual;
    std::vector<double> fitnesses;
    double rankFitness;
    int dominationCount;
    
    SystemInfo(MNIndividual *i):individual(i) {rankFitness = -INFINITY; dominationCount=0;}
    SystemInfo(const SystemInfo &other)
    :individual(other.individual->clone()),
    fitnesses(other.fitnesses),
    rankFitness(other.rankFitness),
    dominationCount(0)
    {}
    
    ~SystemInfo() {
        delete individual;
    }
};

struct NEATSpecies {
    // each generation the members are cleared out, repopulated, and the representative is updated
    MNIndividual *representative;
    std::vector<SystemInfo *>members;
    double totalSharedFitness;
    int speciesNumber; // for data collection
};

// SO DIRTY
struct InnovationInfo {
    //private:
    InnovationInfo(const InnovationInfo& other)
    :data(other.data), innovationNumber(other.innovationNumber), wasCloned(false)
    {
        
    }
    //public:
    MNEdge *data;
    bool wasCloned;
    int innovationNumber;
    
    InnovationInfo(MNEdge *e, bool doClone=true)
    :data(doClone ? e->clone(): e),innovationNumber(-1), wasCloned(doClone){}
    
    InnovationInfo *clone()
    {
        InnovationInfo *newOne = new InnovationInfo(data, true);
        newOne->innovationNumber = innovationNumber;
        return newOne;
    }
    
    
    ~InnovationInfo(){
        if (wasCloned)
            delete data;
    }
    
    bool operator==(const InnovationInfo& other)const {
        return *(data) == *(other.data);
    }
    
    bool operator==(const MNEdge& other)const {
        return *data == other;
    }
};

/// OH GOD I HATE C++ FUNCTORSSSS
// http://stackoverflow.com/questions/258871/how-to-use-find-algorithm-with-a-vector-of-pointers-to-objects-in-c
template <typename T>
struct pointer_values_equal
{
    const T* to_find;
    
    bool operator()(const T* other) const
    {
        return *to_find == *other;
    }
};



#endif
