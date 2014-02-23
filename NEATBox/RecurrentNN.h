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

    DelayEdge(int nodeFrom, int nodeTo):Edge(nodeFrom, nodeTo), delay(0){};
    
    virtual DelayEdge *clone() const {return new DelayEdge(*this);}
};

struct RecurrentNode : Node {
    int minDepth;
    int maxDepth;
    
    RecurrentNode():Node(), minDepth(0), maxDepth(0){};
};

// the internal representation of your graph is up to you
// this is a weighted directed graph, the standard model of a NN
class RecurrentNN : public virtual MNIndividual
{
protected:
    std::vector<Node> nodes;
    std::vector<Edge> edges;
    
    virtual std::vector<MNEdge *> insertNodeOnEdge(Edge &e);


public:
#pragma mark - Necessary for NEATPlusAlgorithm to work

    RecurrentNN(int nInputs, int nOutputs);
    
    virtual void addGeneFromParentSystem(MNIndividual *parent, MNEdge *gene);
    virtual std::vector<MNEdge *> connectionGenome();
    
    virtual double nodeDifference(MNIndividual *other);
    
    virtual void mutateConnectionWeights(double p_m);
   virtual  void mutateNodes(double p_m);
   virtual  std::vector<MNEdge *> createNode();
   virtual  Edge *createConnection();
    
   virtual  double connectionDifference(MNEdge *e1, MNEdge *e2);
    
    virtual RecurrentNN *clone() const {return new RecurrentNN(*this);}

    
#pragma mark - End of necessary methods
    
    // simulate the network for evaluation
    // we generate recurrent networks that may wish to simulate a sequence
    std::vector<std::vector<double> > simulateSequence(const std::vector<std::vector<double> > &inputValues, int delay);
    
    std::string display();
    std::string dotFormat(std::string graphName="RecurrentNN");
    
    long numberOfNodes();
    long numberOfEdges();
    
    std::vector<long> listInputNodes();
    std::vector<long> listOutputNodes();
    
private:
    std::set<long> inputNodes;
    std::set<long> outputNodes;
    
    
    // helpers for simulation
    std::vector<double> nodeOutputsForInputs(std::vector<double> inputs, std::vector<double> lastOutputs);
    std::vector<Edge> inputsToNode(long n);
    std::vector<Edge> outputsFromNode(long n);
    double visitNode(long i, std::set<long> &visitedNodes, std::vector<double> &lastOutputs);

    void cleanup(); // if all inputs to a node are disable, disable all its outputs too!
};

#endif /* defined(__NEATBox__RecurrentNN__) */
