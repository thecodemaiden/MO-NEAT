//
//  BaseNetwork.h
//  NEATBox
//
//  Created by Adeola Bannis on 1/7/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#ifndef __NEATBox__BaseNetwork__
#define __NEATBox__BaseNetwork__

#include <iostream>
#include <vector>
#include <string>
#include <set>
#include "NBUtils.h"
enum ActivationFunc {
    GAUSSIAN_FUNC,
    TANH_FUNC,
    FUNC_SENTINEL,
    STEP_FUNC,
   SIN_FUNC,

};

std::string activationFuncName(ActivationFunc f);

struct Edge {
    int nodeFrom;
    int nodeTo;
    
    double weight;
    
    int innovationNumber; // needed for NEAT
    bool disabled;  // for NEAT
    
    
    Edge(int from, int to)
    :nodeFrom(from), nodeTo(to), innovationNumber(0),disabled(false)
    {
        // random weight
        weight = normallyDistributed();
    }
    
    // equality
    bool operator==(const Edge& other)const {
        return nodeFrom == other.nodeFrom && nodeTo == other.nodeTo;
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

//struct SimReturn {
//    std::vector<std::vector<double> > outputs;
//    SimReturn(std::vector<std::vector<double> > o, long steps):outputs(o){}
//};


class BasicNN
{
    std::vector<Node> nodes;
    std::vector<Edge> edges;
    
    // the internal representation of your graph is up to you
    // this is a weighted directed graph, the standard model of a NN
public:
#pragma mark - Necessary for NEATPlusAlgorithm to work
    static double connectionDifference(const Edge &mine, const Edge &other);

    BasicNN(int nInputs, int nOutputs);
    
    void addGeneFromParentSystem(BasicNN parent, Edge gene);
    std::vector<Edge> connectionGenome();
    
    double nodeDifference(BasicNN other);
    
    void updateInnovationNumber(const Edge &info);

    void mutateConnectionWeight();
    void mutateNode(long n);
    std::vector<Edge> insertNode();
    Edge createConnection();
                            // allows for decay
    
#pragma mark - End of necessary methods
    
    // simulate the network for evaluation
    // we generate recurrent networks that may wish to simulate a sequence
    std::vector<std::vector<double> > simulateSequence(std::vector<std::vector<double> > inputValues, int delay);
    
    std::string display();
    std::string dotFormat(std::string graphName="BasicNN");
    
    long numberOfNodes();
    long numberOfEdges();
    
    std::vector<long> listInputNodes();
    std::vector<long> listOutputNodes();
    
private:
    std::set<long> inputNodes;
    std::set<long> outputNodes;
    
    std::vector<Edge> insertNodeOnEdge(Edge &e);
    
    // helpers for simulation
    std::vector<double> nodeOutputsForInputs(std::vector<double> inputs, std::vector<double> lastOutputs);
    std::vector<Edge> inputsToNode(long n);
    std::vector<Edge> outputsFromNode(long n);
    double visitNode(long i, std::set<long> &visitedNodes, std::vector<double> &lastOutputs);

};

#endif /* defined(__NEATBox__BaseNetwork__) */
