//
//  BayesNet.h
//  NEATBox
//
//  Created by Adeola Bannis on 2/23/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#ifndef __NEATBox__BayesNet__
#define __NEATBox__BayesNet__

#include <iostream>
#include <vector>
#include <map>
#include <NEATBox/NEATBox.h>

class BayesNode {
    int states; // usually just T/F, but not necessarily
    // should we include state names?
    
    std::string name;
    // the cpt of a node looks like
    // (parent1_state, parent2_state, ...) -> prob
    std::map<std::vector<int>, double> cpt;
    
    // make sure we have the right type of entries in our cpt
    // if we must add something, initialize it to equal prob for all parent states
    void updateCPT();
    
public:
    BayesNode(std::string name, int nStates=2):states(nStates), name(name){};

};

// generated each time we ask for the genome - not the worst thing I guess
struct BayesEdge : MNEdge
{
    long nodeFrom;
    long nodeTo;
    
    BayesEdge(long from, long to):nodeFrom(from), nodeTo(to){};
    
    BayesEdge * clone() const { return new BayesEdge(*this);}
    virtual bool operator==(const MNEdge& other)const
    {
       const BayesEdge *be = dynamic_cast<const BayesEdge *>(&other);
        return be && (*this == *be);
    }
    
    virtual bool operator==(const BayesEdge& other)const {
        return (nodeFrom == other.nodeFrom) && (nodeTo == other.nodeTo);
    }
};

class BayesNet : MNIndividual {
    std::vector<BayesNode> nodes;
    std::vector<BayesEdge> edges;
    
public:
    
    long lastEdgeCount; // for numberOfEdges

    BayesNet(std::vector<std::string>nodeNames, std::map<std::vector<int>, double> jdt);
    
    virtual BayesNet * clone() const { return new BayesNet(*this);}
   
    std::vector<MNEdge *> connectionGenome();
    void mutateConnectionWeights(double p_m);
    void mutateNodes(double p_m);
    MNEdge *createConnection();
    std::vector<MNEdge *> createNode();
    long numberOfNodes();
   
    long numberOfEdges();
   
    double nodeDifference(MNIndividual *other);
   
    void addGeneFromParentSystem(MNIndividual *parent, MNEdge *gene);
   
    double connectionDifference(MNEdge *e1, MNEdge *e2);
    
};

#endif /* defined(__NEATBox__BayesNet__) */
