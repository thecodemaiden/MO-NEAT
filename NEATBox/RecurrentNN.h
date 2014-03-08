//
//  BaseNetwork.h
//  NEATBox
//
//  Created by Adeola Bannis on 1/7/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#ifndef __NEATBox__RecurrentNN__
#define __NEATBox__RecurrentNN__

#include <iostream>
#include <vector>
#include <string>
#include <set>
#include "NBUtils.h"
#include "MNIndividual.h"
#include "NNCommon.h"

struct DelayEdge : Edge {
    int delay;

    DelayEdge(long nodeFrom, long nodeTo):Edge(nodeFrom, nodeTo), delay(0){};
    
    virtual DelayEdge *clone() const {return new DelayEdge(*this);}
};

struct RecurrentNode : Node {
    int minDepth;
    int maxDepth;
    
    RecurrentNode():Node(), minDepth(0), maxDepth(0){};
};

// the internal representation of your graph is up to you
// this is a weighted directed graph, the standard model of a NN
class RecurrentNN : public virtual MNIndividual, public CycledNN
{
protected:
    std::vector<Node> nodes;
    std::vector<DelayEdge> edges;
    
    virtual std::vector<MNEdge *> insertNodeOnEdge(long pos);


public:
#pragma mark - Necessary for NEATPlusAlgorithm to work

    RecurrentNN(int nInputs, int nOutputs);
    
    virtual void addGeneFromParentSystem(MNIndividual *parent, MNEdge *gene); // every time we add an edge, we check for cycles
    virtual std::vector<MNEdge *> connectionGenome();
    
    virtual double nodeDifference(MNIndividual *other);
    
    virtual void mutateConnectionWeights(double p_m);
   virtual  void mutateNodes(double p_m);
   virtual  std::vector<MNEdge *> createNode();
   virtual  DelayEdge *createConnection(); // every time we add an edge, we check for cycles
    
   virtual  double connectionDifference(MNEdge *e1, MNEdge *e2);
    
    virtual RecurrentNN *clone() const {return new RecurrentNN(*this);}

    
#pragma mark - End of necessary methods
    
    // we generate recurrent networks that may wish to simulate a sequence
    std::vector<std::vector<double> > simulateSequence(const std::vector<std::vector<double> > &inputValues);
    
    std::string display();
    std::string dotFormat(std::string graphName="RecurrentNN");
    
    long numberOfNodes();
    long numberOfEdges();
    
    std::vector<long> listInputNodes();
    std::vector<long> listOutputNodes();
    
private:
    std::set<long> inputNodes;
    std::set<long> outputNodes;
    
    long maxDelay; // the longest delay of any edge
    
    // -- NOT USED BUT LEAVING THIS FOR NOW
    std::vector<std::vector<long> > cycleSort(); // use Tarjan's algorithm to find cycles, enforce delays in cycles
    struct CycleNode; // helper struct for cycleSort
    // recursive helper function
    void strongConnect(long v, long index, std::vector<long> &stack,
                                    std::vector<RecurrentNN::CycleNode> &nodes,  std::vector<std::vector<long> > &components);
    // -- END --

    void fixCycles(DelayEdge &newEdge); // can't have cycles with no delayed edges

    // helpers for simulation
    std::vector<double> nodeOutputsForInputs(std::vector<double> inputs, std::vector<std::vector<double> > &memory);

    std::vector<DelayEdge> inputsToNode(long n);
    std::vector<DelayEdge> outputsFromNode(long n);
    
    std::vector<DelayEdge> inputsToNode(long n, std::vector<DelayEdge> edgeSet);
    std::vector<DelayEdge> outputsFromNode(long n, std::vector<DelayEdge> edgeSet);

};

#endif /* defined(__NEATBox__RecurrentNN__) */
