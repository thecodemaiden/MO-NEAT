//
//  MNIndividual.h
//  NEATBox
//
//  Created by Adeola Bannis on 2/3/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#ifndef NEATBox_MNIndividual_h
#define NEATBox_MNIndividual_h

#include "MNEdge.h"
#include <vector>

class MNIndividual {
    
public:
    virtual MNIndividual * clone() const = 0;
    virtual ~MNIndividual(){};
    
    virtual std::vector<MNEdge *> connectionGenome() = 0;
    virtual void mutateConnectionWeights(double p_m) = 0;
    virtual void mutateNodes(double p_m) = 0;
    virtual MNEdge *createConnection() = 0;
    virtual std::vector<MNEdge *> createNode() = 0;
    virtual long numberOfNodes() = 0;

    virtual long numberOfEdges() = 0;

    virtual double nodeDifference(MNIndividual *other) = 0;
    
    virtual void addGeneFromParentSystem(MNIndividual *parent, MNEdge *gene) = 0;
    
    virtual double connectionDifference(MNEdge *e1, MNEdge *e2) = 0;
};

#endif
