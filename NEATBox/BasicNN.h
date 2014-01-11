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


enum ActivationFunc {
    STEP_FUNC,
    TANH_FUNC,
    ABS_SIG_FUNC,
    GAUSSIAN_FUNC,
    FUNC_SENTINEL,
};

std::string activationFuncName(ActivationFunc f);

struct Edge {
    int nodeFrom;
    int nodeTo;
    
    double weight;
    
    int innovationNumber; // needed for NEAT
    bool disabled;
        
    Edge(int from, int to)
    :nodeFrom(from), nodeTo(to), innovationNumber(0),disabled(false)
    {
        // random weight
        weight = ((double)rand()/RAND_MAX - 0.5)*2;
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
    
    double threshold = 0.5;
    double activatedVal = 1.0;
    double deactivatedVal = 0.0;
    Node():type(GAUSSIAN_FUNC), indegree(0), outdegree(0), bias(0){};
};

struct SimReturn {
    std::vector<double>outputs;
    bool steady;
    long steps;
    SimReturn(std::vector<double> o, bool steady, long steps):outputs(o), steady(steady), steps(steps){};
};


class BasicNN
{
    std::vector<Node> nodes;
    std::vector<Edge> edges;
    
    // the internal representation of your graph is up to you
    // this is a weighted directed graph, the standard model of a NN
public:
#pragma mark - Necessary for NEATAlgorithm to work
    static double connectionDifference(const Edge &c1, const Edge &c2);

    BasicNN(int nInputs, int nOutputs);
    
    void addGeneFromParentSystem(BasicNN parent, Edge gene);
    std::vector<Edge> connectionGenome();
    void updateInnovationNumber(const Edge &info);

    void mutateConnectionWeight();
    void mutateNode(long n);
    std::pair<Edge, Edge> insertNode();
    Edge createConnection();
    
#pragma mark - End of necessary methods
    // simulate the network for evaluation
    // we can create recurrent networks, and we usually want them to settle
    SimReturn simulateTillEquilibrium(std::vector<double> inputValues, int maxSteps);
    
    std::string display();
    
    long numberOfNodes();
    long numberOfEdges();
    
    std::vector<long> listInputNodes();
    std::vector<long> listOutputNodes();
    
private:
    std::set<long> inputNodes;
    std::set<long> outputNodes;
    
    std::pair<Edge, Edge> insertNodeOnEdge(Edge &e);
    
    // helpers for simulation
    std::vector<double> nodeOuputsForInputs(std::vector<double> inputs, std::vector<double> lastOutputs);
    std::vector<Edge> inputsToNode(long n);
};

#endif /* defined(__NEATBox__BaseNetwork__) */
