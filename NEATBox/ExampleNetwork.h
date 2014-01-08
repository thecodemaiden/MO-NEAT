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

struct Edge {
    int nodeFrom;
    int nodeTo;
    
    double weight;
    
    int innovationNumber; // needed for NEAT
    bool disabled;
        
    Edge(int from, int to)
    :nodeFrom(from), nodeTo(to)
    {
        weight = 1.0; innovationNumber = 0; disabled = false;
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
    int type; // function type: linear, logistic, step
    int indegree; // helps with inserting new connections
    int outdegree;
    
    Node():type(0), indegree(0), outdegree(0){};
   // int name; // for from and to
    //Node(int n) :name(n) {type = 0;}
};



class ExampleNetwork
{
    std::vector<Node> nodes;
    std::vector<Edge> edges;
    
    // the internal representation of your graph is up to you
    // this is a weighted directed graph, the standard model of a NN
public:
    ExampleNetwork();
    
    void addGeneFromParentSystem(ExampleNetwork parent, Edge gene);
    std::vector<Edge> attachmentGenome();
    void mutateAttachmentWeight();
    std::pair<Edge, Edge> insertNode();
    Edge createConnection();
    
    void updateInnovationNumber(const Edge &info);
    
    std::string display();
 
    static double attachmentDifference(const Edge &c1, const Edge &c2);
    static ExampleNetwork smallestSystem();
    
    long numberOfNodes();
    long numberOfEdges();
};

#endif /* defined(__NEATBox__BaseNetwork__) */
