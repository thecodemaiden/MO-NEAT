//
//  NNCommon.h
//  NEATBox
//
//  Created by Adeola Bannis on 2/13/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#ifndef NEATBox_NNCommon_h
#define NEATBox_NNCommon_h

#include <string>
#include <vector>
#include "MNEdge.h"
#include "NBUtils.h"
#include "assert.h"

enum ActivationFunc {

    SIN_FUNC,

  FUNC_SENTINEL,
    GAUSSIAN_FUNC,
    TANH_FUNC,
    STEP_FUNC,
    LIN_FUNC,
};

std::string activationFuncName(ActivationFunc f);

struct Edge : public MNEdge {
    long nodeFrom;
    long nodeTo;
    
    double weight;
    
    bool disabled;  // for NEAT
    
    
    Edge(long from, long to)
    :nodeFrom(from), nodeTo(to),disabled(false), weight(1.0)
    {
    
    }
    
    Edge(const Edge& other)
    :nodeFrom(other.nodeFrom), nodeTo(other.nodeTo), disabled(other.disabled), weight(other.weight)
    {
        
    }
    
    virtual Edge *clone() const {return new Edge(*this);}
    
    // equality
    bool operator==(const Edge& other)const {
        return nodeFrom == other.nodeFrom && nodeTo == other.nodeTo;
    }
    
    virtual bool operator==(const MNEdge& other)const
    {
        const Edge *e = dynamic_cast<const Edge *>(&other);
        return e && (*this == *e);
    }
    
    
    bool operator!=(const Edge& other)const {
        return !(*this==other);
    }
    
    // enforce a sorting order: sort by fromNode, then toNode
    bool operator<(const Edge & other)const {
        if (nodeFrom > other.nodeFrom)
            return  false;
        if (nodeFrom == other.nodeFrom && nodeTo >= other.nodeTo)
            return false;
        return true;
    }
    
};

struct Node {
    ActivationFunc type; // function type: linear, logistic, step
    int indegree; // helps with inserting new connections
    int outdegree;
    
    double bias;
    
    double param1 = 0.0;
    
    double activatedVal = 1.0;
    double deactivatedVal = -1.0;
    Node():indegree(0), outdegree(0), bias(0), param1(0){
        type = (ActivationFunc)arc4random_uniform(FUNC_SENTINEL);
    };
};

double applyActivationFunc(Node n, double inputSum);

// Directed acyclic graph - no cycles
class DAGNN
{
public:
    DAGNN(){}
    virtual ~DAGNN(){}
    // simulate the network for evaluation
    virtual std::vector<double> simulateSequence(const std::vector<double> &inputValues) = 0;
};

// allows cycles, should have delayed edges
// allow time-dependent inputs
class CycledNN
{
public:
    CycledNN(){}
    virtual ~CycledNN(){}
    // we generate recurrent networks that may wish to simulate a sequence
    virtual std::vector<std::vector<double> > simulateSequence(const std::vector<std::vector<double> > &inputValues) = 0;

};

#endif
